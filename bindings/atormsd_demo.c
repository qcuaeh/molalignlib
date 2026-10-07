/* atormsd_demo.c
 *
 * Demo of atormsd_calculate: a C version of the atormsd program for XYZ
 * files. It reads two XYZ files, passes their atom data and coordinates to
 * the library, and prints the results the way atormsd does. Its options,
 * defaults and output follow atormsd. Only molalign.h and libmolalign are
 * needed to build it.
 *
 * Compile (libmolalign is a static Fortran library, so it must come before
 * -lgfortran and -lm; replace /usr/local with the install prefix used, e.g.
 * $HOME/.local for a user-local install, see "User-local install" in the
 * top-level README):
 *
 *   gcc atormsd_demo.c -o atormsd_demo \
 *       -I/usr/local/include/molalignlib -L/usr/local/lib \
 *       -lmolalign -lgfortran -lm
 *
 * Run:
 *   ./atormsd_demo file1.xyz file2.xyz [options]
 *
 * Pass -help for the full list of options.
 *
 * Output (as in atormsd): one line per record with the RMSD and, with
 * -mapping and -remap, the 1-based atom permutation as a comma-separated
 * list. With -align -aligned FILE the aligned cluster 2 is written to FILE
 * instead, one XYZ frame per record with the RMSD in its title line.
 *
 * Dummy atoms (symbol "X") are skipped when reading, as atormsd does: the
 * library has no dummy element.
 *
 * Differences from atormsd: only XYZ files are read and written, and -help
 * is an extra option.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <stdbool.h>
#include "molalign.h"

/* ---- Command line: single-dash long options ---- */

typedef struct {
    const char *name;
    int has_arg;   /* 0 = no argument, 1 = required argument */
    int val;
} long_opt_t;

typedef struct {
    const char *desc;
    const char *arg;   /* argument metavar, or NULL */
} opt_info_t;

/* Case-insensitive string equality (atormsd lowercases its options). */
static bool equal_nocase(const char *a, const char *b)
{
    for (; *a && *b; a++, b++)
        if (tolower((unsigned char)*a) != tolower((unsigned char)*b)) return false;
    return *a == *b;
}

/* Parse the next option from argv, recognising options in any position
 * relative to the positional arguments, as atormsd does: every token that
 * starts with '-' is an option, and anything else is a file path.
 *
 * Positional tokens are collected into posargs as they are encountered.
 * *npos counts every one of them (it must be 0 before the first call), but
 * only the first max_pos are stored, so the caller can detect too many
 * without posargs overflowing.
 *
 * Sets *option to the option token as typed, and *arg to the option
 * argument when has_arg=1. Like atormsd, an option argument may not start
 * with '-'.
 * Returns the matched val, 0 on error (unknown option or missing
 * argument), -1 when done. */
static int parse_long_opt(int argc, char **argv, int *argi,
                          const char **option, const char **arg,
                          const long_opt_t *opts,
                          const char **posargs, int max_pos, int *npos)
{
    /* Collect positional tokens */
    while (*argi < argc && argv[*argi][0] != '-') {
        if (*npos < max_pos) posargs[*npos] = argv[*argi];
        (*npos)++;
        (*argi)++;
    }
    if (*argi >= argc) return -1;

    const char *p = argv[(*argi)++];
    *option = p;
    for (int i = 0; opts[i].name; i++) {
        if (!equal_nocase(p + 1, opts[i].name)) continue;   /* exact match only */
        if (opts[i].has_arg) {
            if (*argi >= argc || argv[*argi][0] == '-') {
                fprintf(stderr, "Option %s requires an argument\n", p);
                return 0;
            }
            *arg = argv[(*argi)++];
        }
        return opts[i].val;
    }

    fprintf(stderr, "Unknown option: %s\n", p);
    return 0;
}

/* Integer argument of option, at least min_value. Returns 0 on success. */
static int read_int_optarg(const char *option, const char *arg, int min_value, int *value)
{
    char *end;
    long v;

    errno = 0;
    v = strtol(arg, &end, 10);
    if (end == arg || *end != '\0' || errno == ERANGE || v < -2147483647L || v > 2147483647L) {
        fprintf(stderr, "Option %s requires an integer argument\n", option);
        return 1;
    }
    if (v < min_value) {
        fprintf(stderr, "Option %s requires an integer argument of at least %d\n",
                option, min_value);
        return 1;
    }
    *value = (int)v;
    return 0;
}

/* Real argument of option. Returns 0 on success. */
static int read_real_optarg(const char *option, const char *arg, double *value)
{
    char *end;

    errno = 0;
    *value = strtod(arg, &end);
    if (end == arg || *end != '\0' || errno == ERANGE) {
        fprintf(stderr, "Option %s requires a numeric argument\n", option);
        return 1;
    }
    return 0;
}

/* Print a formatted option list to stderr. */
static void print_options(const long_opt_t *opts, const opt_info_t *info)
{
    char head[32];

    for (int i = 0; opts[i].name; i++) {
        const opt_info_t *oi = &info[opts[i].val];
        if (oi->arg)
            snprintf(head, sizeof(head), "-%s %s", opts[i].name, oi->arg);
        else
            snprintf(head, sizeof(head), "-%s", opts[i].name);
        fprintf(stderr, "  %-18s %s\n", head, oi->desc);
    }
}

/* ---- Element table ---- */

/* Symbols indexed by element number, as in the library's element table:
 * the real elements up to Lr (103), then the Lennard-Jones pseudo-element
 * "LJ" (104). Entry 0 is not an element: it is "X", the symbol written for
 * padding atoms. */
static const char *atomic_symbols[] = {
    "X",                                                                  /* padding */
    "H",  "He",                                                           /* 1-2    */
    "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",                     /* 3-10   */
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",                     /* 11-18  */
    "K",  "Ca",                                                           /* 19-20  */
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",         /* 21-30  */
    "Ga", "Ge", "As", "Se", "Br", "Kr",                                  /* 31-36  */
    "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",         /* 37-46  */
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",                     /* 47-54  */
    "Cs", "Ba",                                                           /* 55-56  */
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",                     /* 57-64  */
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",                           /* 65-71  */
    "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",              /* 72-80  */
    "Tl", "Pb", "Bi", "Po", "At", "Rn",                                  /* 81-86  */
    "Fr", "Ra",                                                           /* 87-88  */
    "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm",                     /* 89-96  */
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",                           /* 97-103 */
    "LJ"                                                                 /* 104    */
};
static const int n_elems = (int)(sizeof(atomic_symbols) / sizeof(atomic_symbols[0])) - 1;

/* Element number (1..n_elems) of an element symbol (case-insensitive), or
 * -1 if unknown. */
static int elnum_lookup(const char *elsym)
{
    for (int z = 1; z <= n_elems; z++)
        if (equal_nocase(elsym, atomic_symbols[z])) return z;
    return -1;
}

/* ---- XYZ files ---- */

/* Split an atom label into element number and group, as atormsd does.
 *
 * The label is matched case-insensitively. Its leading letters are the
 * element symbol and anything after them must be digits, giving the group
 * (0 if absent), e.g. "C", "fe", "H12". Groups only matter with -atomtype.
 * "X" marks a dummy atom, which is not an element: *is_dummy is then set
 * and *elnum must not be used.
 *
 * Returns 0 on success, -1 on invalid label. */
static int parse_label(const char *label, int *elnum, int *group, bool *is_dummy)
{
    char elsym[32];
    int i, m = 0;

    while (isalpha((unsigned char)label[m])) m++;
    for (i = m; label[i]; i++)
        if (!isdigit((unsigned char)label[i])) {
            fprintf(stderr, "Invalid atomic label: %s\n", label);
            return -1;
        }

    *group = atoi(&label[m]);   /* 0 when there is no digit suffix */

    /* Also rejects an empty or over-long symbol: neither is in the table. */
    snprintf(elsym, sizeof(elsym), "%.*s", m, label);
    *is_dummy = equal_nocase(elsym, "X");
    if (*is_dummy) return 0;
    *elnum = elnum_lookup(elsym);
    if (*elnum < 0) {
        fprintf(stderr, "Unknown element symbol: %s\n", elsym);
        return -1;
    }
    return 0;
}

/* Format of a path, taken from its extension as atormsd does: the text
 * after the last '.' of the file name. A file name without a '.' has no
 * base name (e.g. "xyz") and the whole name is the format. */
static const char *path_format(const char *path, bool *has_basename)
{
    const char *filename = strrchr(path, '/');
    const char *dot;

    filename = filename ? filename + 1 : path;
    dot = strrchr(filename, '.');
    *has_basename = (dot != NULL);
    return dot ? dot + 1 : filename;
}

/* Read an XYZ file into packed atom_data and coords arrays. Dummy atoms
 * ("X") are skipped, so *n_atoms_out counts the real atoms only.
 *
 * atom_data_out receives a heap-allocated array of length n_atoms*2:
 *   [elnum0, group0, elnum1, group1, ...]
 * coords_out receives a heap-allocated array of length n_atoms*3.
 * The caller is responsible for freeing both arrays.
 *
 * Returns 0 on success, non-zero on error (message written to stderr). */
static int read_xyz(const char *path,
                    int *n_atoms_out,
                    int **atom_data_out,
                    double **coords_out)
{
    FILE *fp;
    int n_atoms, n_real = 0, elnum, group;
    bool has_basename, is_dummy;
    char label[32], line[256];
    int *atom_data = NULL;
    double *coords = NULL;

    const char *in_format = path_format(path, &has_basename);
    if (!has_basename || strcmp(in_format, "xyz") != 0) {
        fprintf(stderr, "File format \"%s\" is not supported (this demo only reads .xyz files)\n",
                in_format);
        return 1;
    }

    fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "Error opening %s for reading\n", path); return 1; }

    if (fscanf(fp, " %d", &n_atoms) != 1) {
        fprintf(stderr, "Invalid XYZ format in %s\n", path);
        goto fail;
    }
    if (n_atoms <= 0) {
        fprintf(stderr, "File %s contains no atoms\n", path);
        goto fail;
    }

    /* Skip the rest of the count line and the title line. */
    if (!fgets(line, sizeof(line), fp) || !fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "Invalid XYZ format in %s\n", path);
        goto fail;
    }

    atom_data = malloc((size_t)n_atoms * 2 * sizeof(int));
    coords    = malloc((size_t)n_atoms * 3 * sizeof(double));
    if (!atom_data || !coords) {
        fprintf(stderr, "Error: out of memory\n");
        goto fail;
    }

    for (int i = 0; i < n_atoms; i++) {
        double *xyz = &coords[n_real*3];
        if (fscanf(fp, " %31s %lf %lf %lf", label, &xyz[0], &xyz[1], &xyz[2]) != 4) {
            fprintf(stderr, "Invalid XYZ format in %s (atom line %d)\n", path, i+1);
            goto fail;
        }
        if (parse_label(label, &elnum, &group, &is_dummy) != 0) goto fail;
        if (is_dummy) continue;   /* skip dummy atoms */
        atom_data[n_real*2]     = elnum;
        atom_data[n_real*2 + 1] = group;
        n_real++;
    }
    if (n_real == 0) {
        fprintf(stderr, "File %s contains no atoms\n", path);
        goto fail;
    }

    fclose(fp);
    *n_atoms_out   = n_real;
    *atom_data_out = atom_data;
    *coords_out    = coords;
    return 0;

fail:
    free(atom_data);
    free(coords);
    fclose(fp);
    return 1;
}

/* Write one XYZ frame of cluster 2, transformed by transform and with its
 * atoms reordered by full_atomperm1, like atormsd -aligned. Line i holds
 * atom full_atomperm1[i] (0-based) of the padded cluster 2: indices
 * >= n_atoms2 are padding atoms, written as "X" with the origin as their
 * input coordinates, as in the library. */
static void write_xyz_frame(FILE *fp, const char *title, int n_padding,
                            int n_atoms2, const int *atom_data2, const double *coords2,
                            const int *full_atomperm1, const double *transform)
{
    fprintf(fp, "%d\n%s\n", n_padding, title);

    for (int i = 0; i < n_padding; i++) {
        int j = full_atomperm1[i];
        const double origin[3] = {0.0, 0.0, 0.0};
        const double *p = j < n_atoms2 ? &coords2[j*3] : origin;
        int elnum = j < n_atoms2 ? atom_data2[j*2] : 0;
        double q[3];

        /* p_out = R * p_in + t, with the row-major 4x4 transform */
        for (int k = 0; k < 3; k++)
            q[k] = transform[k*4]*p[0] + transform[k*4+1]*p[1]
                 + transform[k*4+2]*p[2] + transform[k*4+3];

        fprintf(fp, "%-2s  %12.6f  %12.6f  %12.6f\n", atomic_symbols[elnum], q[0], q[1], q[2]);
    }
}

/* ---- Demo program ---- */

/* Option identifiers (also indices into opt_info) */
enum {
    OPT_ALIGN = 1, OPT_REMAP, OPT_PRUNETOL, OPT_ATOMTYPE,
    OPT_HEAVY, OPT_MASSWEIGHT, OPT_MIRROR, OPT_ATOFREQ,
    OPT_MAXTRIALS, OPT_MAXRECS, OPT_ALIGNED, OPT_MAPPING,
    OPT_STATS, OPT_RANDOM, OPT_HELP
};

/* Same options, in the same order, as the atormsd program, plus -help */
static const long_opt_t long_options[] = {
    {"align",      0, OPT_ALIGN},
    {"remap",      0, OPT_REMAP},
    {"prunetol",   1, OPT_PRUNETOL},
    {"atomtype",   0, OPT_ATOMTYPE},
    {"heavy",      0, OPT_HEAVY},
    {"massweight", 0, OPT_MASSWEIGHT},
    {"mirror",     0, OPT_MIRROR},
    {"atofreq",    1, OPT_ATOFREQ},
    {"maxtrials",  1, OPT_MAXTRIALS},
    {"maxrecs",    1, OPT_MAXRECS},
    {"aligned",    1, OPT_ALIGNED},
    {"mapping",    0, OPT_MAPPING},
    {"stats",      0, OPT_STATS},
    {"random",     0, OPT_RANDOM},
    {"help",       0, OPT_HELP},
    {NULL, 0, 0}
};

static const opt_info_t opt_info[] = {
    [OPT_ALIGN]      = {"Align atoms to minimise the RMSD",                                NULL  },
    [OPT_REMAP]      = {"Remap atoms to minimise the RMSD",                                NULL  },
    [OPT_PRUNETOL]   = {"Prune atom pairs with tolerance TOL (Angstrom)",                  "TOL" },
    [OPT_ATOMTYPE]   = {"Only match atoms with the same label (e.g. C1, C2)",              NULL  },
    [OPT_HEAVY]      = {"Ignore hydrogen atoms",                                           NULL  },
    [OPT_MASSWEIGHT] = {"Use mass-weighted coordinates",                                   NULL  },
    [OPT_MIRROR]     = {"Reflect cluster 2 (x -> -x) before comparison",                   NULL  },
    [OPT_ATOFREQ]    = {"Stop once the best solution is found N times (default: 10)",      "N"   },
    [OPT_MAXTRIALS]  = {"Stop after at most N random orientations (default: 10000)",       "N"   },
    [OPT_MAXRECS]    = {"Record the N lowest RMSDs found (default: 1)",                    "N"   },
    [OPT_ALIGNED]    = {"Write the aligned cluster 2 to FILE (XYZ; needs -align)",         "FILE"},
    [OPT_MAPPING]    = {"Print the atom permutation (1-based; needs -remap)",              NULL  },
    [OPT_STATS]      = {"Print optimisation statistics (needs -align -remap)",             NULL  },
    [OPT_RANDOM]     = {"Seed the random-number generator from the system clock",          NULL  },
    [OPT_HELP]       = {"Show this help message",                                          NULL  },
};

static void print_usage(const char *prog)
{
    fprintf(stderr, "Usage: %s file1.xyz file2.xyz [options]\n\n"
                    "Calculate the RMSD between two atom clusters (XYZ format).\n\n"
                    "Options:\n", prog);
    print_options(long_options, opt_info);
    fputc('\n', stderr);
}

int main(int argc, char **argv)
{
    const char *arg = NULL;
    const char *option = NULL;
    bool align_flag = false, remap_flag = false, heavy_flag = false, massweight_flag = false;
    bool mirror_flag = false, useatomtype_flag = false;
    bool printstats_flag = false, random_flag = false;
    bool printmapping_flag = false, write_aligned = false;
    const char *aligned_path = NULL;
    bool pruning_flag = false;
    double prune_tol = 0.0;
    int ato_freq = 10, max_trials = 10000;
    int max_records = 1;
    int argi = 1, opt;

    /* Positional arguments (the two file paths) */
    const char *posargs[2] = {NULL, NULL};
    int npos = 0;

    while ((opt = parse_long_opt(argc, argv, &argi, &option, &arg, long_options,
                                 posargs, 2, &npos)) != -1) {
        switch (opt) {
        case OPT_ALIGN:      align_flag     = true; break;
        case OPT_REMAP:      remap_flag     = true; break;
        case OPT_PRUNETOL:
            pruning_flag = true;
            if (read_real_optarg(option, arg, &prune_tol) != 0) return 1;
            break;
        case OPT_ATOMTYPE:   useatomtype_flag = true; break;
        case OPT_HEAVY:      heavy_flag     = true; break;
        case OPT_MASSWEIGHT: massweight_flag      = true; break;
        case OPT_MIRROR:     mirror_flag    = true; break;
        case OPT_ATOFREQ:
            if (read_int_optarg(option, arg, 1, &ato_freq) != 0) return 1;
            break;
        case OPT_MAXTRIALS:
            if (read_int_optarg(option, arg, 1, &max_trials) != 0) return 1;
            break;
        case OPT_MAXRECS:
            if (read_int_optarg(option, arg, 1, &max_records) != 0) return 1;
            break;
        case OPT_ALIGNED:
            write_aligned = true;
            aligned_path = arg;
            break;
        case OPT_MAPPING:    printmapping_flag = true; break;
        case OPT_STATS:      printstats_flag      = true; break;
        case OPT_RANDOM:     random_flag     = true; break;
        case OPT_HELP:       print_usage(argv[0]); return 0;
        default:             print_usage(argv[0]); return 1;
        }
    }

    switch (npos) {
    case 0:  fprintf(stderr, "File paths are missing\n"); print_usage(argv[0]); return 1;
    case 1:  fprintf(stderr, "Too few file paths\n");     print_usage(argv[0]); return 1;
    case 2:  break;
    default: fprintf(stderr, "Too many file paths\n");    print_usage(argv[0]); return 1;
    }

    /* Everything below is released at `done`; free(NULL) is a no-op, so
     * every exit path can jump there regardless of how far it got. */
    int status = 1;
    int n_atoms1 = 0, n_atoms2 = 0, n_padding = 0;
    int *atom_data1 = NULL, *atom_data2 = NULL;
    double *coords1 = NULL, *coords2 = NULL;
    double *rmsd_list = NULL, *transform_list = NULL;
    int *mapping_list = NULL;
    int n_records = 0, error_code = MOLALIGN_SUCCESS;
    FILE *aligned_file = NULL;

    if (read_xyz(posargs[0], &n_atoms1, &atom_data1, &coords1) != 0) goto done;
    if (read_xyz(posargs[1], &n_atoms2, &atom_data2, &coords2) != 0) goto done;

    /* As in atormsd, -aligned only takes effect with -align; it then
     * replaces the RMSD output on stdout. A file name without a base name
     * (e.g. "xyz") writes to stdout. */
    write_aligned = write_aligned && align_flag;
    if (write_aligned) {
        bool has_basename;
        const char *out_format = path_format(aligned_path, &has_basename);
        if (strcmp(out_format, "xyz") != 0) {
            fprintf(stderr, "File format \"%s\" is not supported (this demo only writes .xyz files)\n",
                    out_format);
            goto done;
        }
        if (has_basename) {
            aligned_file = fopen(aligned_path, "w");
            if (!aligned_file) {
                fprintf(stderr, "Can't open %s for writing\n", aligned_path);
                goto done;
            }
        }
    }

    /* Output buffers for up to max_records solutions. Each permutation has
     * n_padding = max(n_atoms1, n_atoms2) entries: the clusters may differ
     * in size with -heavy, and the smaller one is then padded with padding
     * atoms. */
    n_padding = n_atoms1 > n_atoms2 ? n_atoms1 : n_atoms2;
    rmsd_list      = malloc((size_t)max_records * sizeof(double));
    mapping_list   = malloc((size_t)max_records * (size_t)n_padding * sizeof(int));
    transform_list = malloc((size_t)max_records * 16 * sizeof(double));
    if (!rmsd_list || !mapping_list || !transform_list) {
        fprintf(stderr, "Error: out of memory\n");
        goto done;
    }

    /* The library prints the -stats output through Fortran I/O */
    fflush(stdout);

    atormsd_calculate(
        n_atoms1, atom_data1, coords1,
        n_atoms2, atom_data2, coords2,
        align_flag, remap_flag, heavy_flag, massweight_flag,
        mirror_flag, useatomtype_flag,
        printstats_flag, random_flag,
        pruning_flag, prune_tol, ato_freq, max_trials,
        max_records,
        rmsd_list, mapping_list,
        transform_list, &n_records, &error_code);

    if (error_code != MOLALIGN_SUCCESS) {
        switch (error_code) {
        case MOLALIGN_ERROR_INVALID_BOUND:
            fprintf(stderr, "Error: -atofreq, -maxtrials and -maxrecs must be at least 1\n"); break;
        case MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER:
            fprintf(stderr, "Error: atomic number out of range\n"); break;
        case MOLALIGN_ERROR_NOT_ISOMERS:
            fprintf(stderr, "These molecules are not isomers\n"); break;
        case MOLALIGN_ERROR_ATOM_TYPE_MISMATCH:
            fprintf(stderr, "Atoms do not match\n"); break;
        case MOLALIGN_ERROR_ASSIGNMENT_FAILED:
            fprintf(stderr, "Error: Assignment failed (pruning tolerance might be too tight)\n"); break;
        default:
            fprintf(stderr, "Error: error code %d\n", error_code); break;
        }
        status = error_code;
        goto done;
    }

    for (int r = 0; r < n_records; r++) {
        double rmsd = rmsd_list[r];
        const double *transform = &transform_list[(size_t)r * 16];
        /* Entry i (0-based) is the atom of cluster 2 placed on line i of
         * cluster 1. Values >= n_atoms2 are padding atoms of cluster 2,
         * and entries i >= n_atoms1 hold the extra atoms of cluster 2. */
        const int *full_atomperm1 = &mapping_list[(size_t)r * n_padding];

        if (write_aligned) {
            char title[64];
            snprintf(title, sizeof(title), "rmsd=%.6f", rmsd);
            write_xyz_frame(aligned_file ? aligned_file : stdout, title, n_padding,
                            n_atoms2, atom_data2, coords2, full_atomperm1, transform);
        } else {
            printf("%.6f", rmsd);
            if (printmapping_flag && remap_flag) {
                /* 1-based, comma-separated, as printed by atormsd */
                putchar(' ');
                for (int i = 0; i < n_padding; i++)
                    printf(i ? ",%d" : "%d", full_atomperm1[i] + 1);
            }
            putchar('\n');
        }
    }
    status = 0;

done:
    if (aligned_file) fclose(aligned_file);
    free(atom_data1); free(coords1);
    free(atom_data2); free(coords2);
    free(rmsd_list); free(mapping_list); free(transform_list);
    return status;
}
