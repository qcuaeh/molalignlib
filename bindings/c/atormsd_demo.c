/* atormsd_demo.c
 *
 * Demo of ato_rmsd: a C version of the atormsd program for XYZ
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
 * Run it without arguments for the full list of options.
 *
 * Output (as in atormsd): one line per record with the RMSD and, with
 * -printmapping, the 1-based atom permutation as a comma-separated list. With
 * -coords xyz the coordinates of cluster 2 (aligned with -align, reordered
 * with -remap) are written to stdout instead, one XYZ frame per record with
 * the RMSD in its title line. With -coords xyz -savecoords PREFIX, each record
 * is written instead to the file PREFIX_<n>.xyz (n = 1, 2, ...) as two
 * frames, cluster 1 and then cluster 2, and the RMSD lines are printed as
 * usual.
 *
 * Dummy atoms (symbol "X") are skipped when reading, as atormsd does: the
 * library has no dummy element. As in atormsd, a file path given as just
 * the format name, "xyz", reads standard input instead of a file, e.g.
 *   ./atormsd_demo xyz file2.xyz -align -remap < file1.xyz
 *
 * Differences from atormsd: only XYZ files are read and written.
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
 * (0 if absent), e.g. "C", "fe", "H12". Groups only matter with -useatomtype.
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
 * base name, and the whole name is the format (e.g. "xyz"): atormsd then
 * reads standard input. Returns NULL, after printing a message, for a path
 * that ends in '/' or a file name with an empty base name or extension
 * (e.g. ".xyz" or "mol."), which atormsd rejects. */
static const char *path_format(const char *path, bool *has_basename)
{
    const char *filename = strrchr(path, '/');
    const char *dot;

    filename = filename ? filename + 1 : path;
    if (*filename == '\0') {
        fprintf(stderr, "Invalid path\n");
        return NULL;
    }
    dot = strrchr(filename, '.');
    *has_basename = (dot != NULL);
    if (dot && (dot == filename || dot[1] == '\0')) {
        fprintf(stderr, "Invalid file name %s\n", filename);
        return NULL;
    }
    return dot ? dot + 1 : filename;
}

/* Read an XYZ file into packed atom_data and coords arrays. Dummy atoms
 * ("X") are skipped, so *n_atoms_out counts the real atoms only.
 *
 * atom_data_out receives a heap-allocated array of length n_atoms*2:
 *   [elnum0, group0, elnum1, group1, ...]
 * coords_out receives a heap-allocated array of length n_atoms*3.
 * The caller is responsible for freeing both arrays. If title_out is not
 * NULL it receives the title line (at most 255 characters, without the
 * line break). As in atormsd, a path that is just a format name (no base
 * name, e.g. "xyz") reads standard input.
 *
 * Returns 0 on success, non-zero on error (message written to stderr). */
static int read_xyz(const char *path,
                    int *n_atoms_out,
                    int **atom_data_out,
                    double **coords_out,
                    char *title_out)
{
    FILE *fp;
    int n_atoms, n_real = 0, elnum, group;
    bool has_basename, is_dummy;
    char label[32], line[256];
    int *atom_data = NULL;
    double *coords = NULL;

    const char *in_format = path_format(path, &has_basename);
    if (!in_format) return 1;
    if (strcmp(in_format, "xyz") != 0) {
        fprintf(stderr, "File format \"%s\" is not supported (this demo only reads .xyz files)\n",
                in_format);
        return 1;
    }

    if (has_basename) {
        fp = fopen(path, "r");
        if (!fp) { fprintf(stderr, "Error opening %s for reading\n", path); return 1; }
    } else {
        fp = stdin;
    }

    if (fscanf(fp, " %d", &n_atoms) != 1) {
        fprintf(stderr, "Invalid XYZ format in %s\n", path);
        goto fail;
    }
    if (n_atoms <= 0) {
        fprintf(stderr, "File %s contains no atoms\n", path);
        goto fail;
    }

    /* Skip the rest of the count line and read the title line. */
    if (!fgets(line, sizeof(line), fp) || !fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "Invalid XYZ format in %s\n", path);
        goto fail;
    }
    if (title_out) {
        line[strcspn(line, "\r\n")] = '\0';
        strcpy(title_out, line);
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

    if (fp != stdin) fclose(fp);
    *n_atoms_out   = n_real;
    *atom_data_out = atom_data;
    *coords_out    = coords;
    return 0;

fail:
    free(atom_data);
    free(coords);
    if (fp != stdin) fclose(fp);
    return 1;
}

/* Write one XYZ frame of a cluster padded to n_padding atoms, transformed
 * by transform and with its atoms reordered by mapping, like atormsd
 * -coords. Line i holds atom mapping[i] (0-based) of the padded cluster:
 * indices >= n_atoms are padding atoms, written as "X" with the origin as
 * their input coordinates, as in the library. A NULL mapping keeps the
 * input order and a NULL transform keeps the input coordinates (as for
 * cluster 1 with -savecoords). */
static void write_xyz_frame(FILE *fp, const char *title, int n_padding,
                            int n_atoms, const int *atom_data, const double *coords,
                            const int *mapping, const double *transform)
{
    fprintf(fp, "%d\n%s\n", n_padding, title);

    for (int i = 0; i < n_padding; i++) {
        int j = mapping ? mapping[i] : i;
        const double origin[3] = {0.0, 0.0, 0.0};
        const double *p = j < n_atoms ? &coords[j*3] : origin;
        int elnum = j < n_atoms ? atom_data[j*2] : 0;
        double q[3];

        /* p_out = R * p_in + t, with the row-major 4x4 transform */
        for (int k = 0; k < 3; k++)
            q[k] = transform ? transform[k*4]*p[0] + transform[k*4+1]*p[1]
                             + transform[k*4+2]*p[2] + transform[k*4+3]
                             : p[k];

        fprintf(fp, "%-3s  %12.6f  %12.6f  %12.6f\n", atomic_symbols[elnum], q[0], q[1], q[2]);
    }
}

/* ---- Demo program ---- */

/* Option identifiers (also indices into opt_info) */
enum {
    OPT_ALIGN = 1, OPT_REMAP, OPT_PRUNE, OPT_USEATOMTYPE,
    OPT_HEAVY, OPT_MASSWEIGHT, OPT_MIRROR, OPT_MINFREQ,
    OPT_MAXTRIALS, OPT_MAXRECS, OPT_COORDS, OPT_SAVECOORDS, OPT_PRINTMAPPING,
    OPT_PRINTSTATS, OPT_RANDOM
};

/* Same options, in the same order, as the atormsd program */
static const long_opt_t long_options[] = {
    {"align",        0, OPT_ALIGN},
    {"remap",        0, OPT_REMAP},
    {"prune",        1, OPT_PRUNE},
    {"useatomtype",  0, OPT_USEATOMTYPE},
    {"heavy",        0, OPT_HEAVY},
    {"massweight",   0, OPT_MASSWEIGHT},
    {"mirror",       0, OPT_MIRROR},
    {"minfreq",      1, OPT_MINFREQ},
    {"maxtrials",    1, OPT_MAXTRIALS},
    {"maxrecs",      1, OPT_MAXRECS},
    {"coords",       1, OPT_COORDS},
    {"savecoords",   1, OPT_SAVECOORDS},
    {"printmapping", 0, OPT_PRINTMAPPING},
    {"printstats",   0, OPT_PRINTSTATS},
    {"random",       0, OPT_RANDOM},
    {NULL, 0, 0}
};

static const opt_info_t opt_info[] = {
    [OPT_ALIGN]        = {"Align atoms to minimise the RMSD",                                NULL  },
    [OPT_REMAP]        = {"Remap atoms to minimise the RMSD",                                NULL  },
    [OPT_PRUNE]        = {"Prune atom pairs with tolerance TOL (Angstrom)",                  "TOL" },
    [OPT_USEATOMTYPE]  = {"Only match atoms with the same label (e.g. C1, C2)",              NULL  },
    [OPT_HEAVY]        = {"Ignore hydrogen atoms",                                           NULL  },
    [OPT_MASSWEIGHT]   = {"Use mass-weighted coordinates",                                   NULL  },
    [OPT_MIRROR]       = {"Reflect cluster 2 (x -> -x) before comparison",                   NULL  },
    [OPT_MINFREQ]      = {"Stop once the best solution is found N times (default: 10)",      "N"   },
    [OPT_MAXTRIALS]    = {"Stop after at most N random orientations (default: 10000)",       "N"   },
    [OPT_MAXRECS]      = {"Record the N lowest RMSDs found (default: 1)",                    "N"   },
    [OPT_COORDS]       = {"Write the coordinates of cluster 2 to stdout in FORMAT (xyz)",   "FORMAT"},
    [OPT_SAVECOORDS]   = {"Write each cluster 1/cluster 2 pair to PREFIX_<n>.xyz (needs -coords)", "PREFIX"},
    [OPT_PRINTMAPPING] = {"Print the 1-based atom permutation after each RMSD",              NULL  },
    [OPT_PRINTSTATS]   = {"Print optimisation statistics (needs -align -remap)",             NULL  },
    [OPT_RANDOM]       = {"Seed the random-number generator from the system clock",          NULL  },
};

static void print_usage(const char *prog)
{
    fprintf(stderr, "Usage: %s file1.xyz file2.xyz [options]\n\n"
                    "Calculate the RMSD between two atom clusters (XYZ format).\n\n"
                    "A path that is just the format name, xyz, reads standard input.\n\n"
                    "Options:\n", prog);
    print_options(long_options, opt_info);
}

int main(int argc, char **argv)
{
    const char *arg = NULL;
    const char *option = NULL;
    bool align_flag = false, remap_flag = false, heavy_flag = false, massweight_flag = false;
    bool mirror_flag = false, useatomtype_flag = false;
    bool printstats_flag = false, random_flag = false;
    bool printmapping_flag = false, coords_flag = false, savecoords_flag = false;
    const char *pairs_prefix = NULL;
    bool prune_flag = false;
    double prune_tol = 0.0;
    int min_freq = 10, max_trials = 10000;
    int max_recs = 1;
    int argi = 1, opt;

    /* Print the usage when no arguments are given, as atormsd does */
    if (argc < 2) {
        print_usage(argv[0]);
        return 0;
    }

    /* Positional arguments (the two file paths) */
    const char *posargs[2] = {NULL, NULL};
    int npos = 0;

    while ((opt = parse_long_opt(argc, argv, &argi, &option, &arg, long_options,
                                 posargs, 2, &npos)) != -1) {
        switch (opt) {
        case OPT_ALIGN:      align_flag     = true; break;
        case OPT_REMAP:      remap_flag     = true; break;
        case OPT_PRUNE:
            prune_flag = true;
            if (read_real_optarg(option, arg, &prune_tol) != 0) return 1;
            break;
        case OPT_USEATOMTYPE:   useatomtype_flag = true; break;
        case OPT_HEAVY:      heavy_flag     = true; break;
        case OPT_MASSWEIGHT: massweight_flag      = true; break;
        case OPT_MIRROR:     mirror_flag    = true; break;
        case OPT_MINFREQ:
            if (read_int_optarg(option, arg, 1, &min_freq) != 0) return 1;
            break;
        case OPT_MAXTRIALS:
            if (read_int_optarg(option, arg, 1, &max_trials) != 0) return 1;
            break;
        case OPT_MAXRECS:
            if (read_int_optarg(option, arg, 1, &max_recs) != 0) return 1;
            break;
        case OPT_COORDS:
            /* The argument is a format, not a path; output goes to stdout */
            if (!equal_nocase(arg, "xyz")) {
                fprintf(stderr, "Option %s requires an output format "
                                "(this demo only writes xyz), not \"%s\"\n", option, arg);
                return 1;
            }
            coords_flag = true;
            break;
        case OPT_SAVECOORDS:
            savecoords_flag = true;
            pairs_prefix = arg;
            break;
        case OPT_PRINTMAPPING: printmapping_flag = true; break;
        case OPT_PRINTSTATS: printstats_flag      = true; break;
        case OPT_RANDOM:     random_flag     = true; break;
        default:             return 1;   /* message already printed */
        }
    }

    /* -savecoords takes the output format from -coords */
    if (savecoords_flag && !coords_flag) {
        fprintf(stderr, "Option -savecoords requires -coords\n");
        return 1;
    }

    switch (npos) {
    case 0:  fprintf(stderr, "File paths are missing\n"); return 1;
    case 1:  fprintf(stderr, "Too few file paths\n");     return 1;
    case 2:  break;
    default: fprintf(stderr, "Too many file paths\n");    return 1;
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
    char title1[256];

    if (read_xyz(posargs[0], &n_atoms1, &atom_data1, &coords1, title1) != 0) goto done;
    if (read_xyz(posargs[1], &n_atoms2, &atom_data2, &coords2, NULL) != 0) goto done;

    /* Output buffers for up to max_recs solutions. Each permutation has
     * n_padding = max(n_atoms1, n_atoms2) entries: the clusters may differ
     * in size with -heavy, and the smaller one is then padded with padding
     * atoms. */
    n_padding = n_atoms1 > n_atoms2 ? n_atoms1 : n_atoms2;
    rmsd_list      = malloc((size_t)max_recs * sizeof(double));
    mapping_list   = malloc((size_t)max_recs * (size_t)n_padding * sizeof(int));
    transform_list = malloc((size_t)max_recs * 16 * sizeof(double));
    if (!rmsd_list || !mapping_list || !transform_list) {
        fprintf(stderr, "Error: out of memory\n");
        goto done;
    }

    /* The library prints the -printstats output through Fortran I/O */
    fflush(stdout);

    ato_rmsd(
        n_atoms1, atom_data1, coords1,
        n_atoms2, atom_data2, coords2,
        align_flag, remap_flag, heavy_flag, massweight_flag,
        mirror_flag, useatomtype_flag,
        printstats_flag, random_flag,
        prune_flag, prune_tol, min_freq, max_trials,
        max_recs,
        rmsd_list, mapping_list,
        transform_list, &n_records, &error_code);

    if (error_code != MOLALIGN_SUCCESS) {
        switch (error_code) {
        case MOLALIGN_ERROR_INVALID_BOUND:
            fprintf(stderr, "Error: -minfreq, -maxtrials and -maxrecs must be at least 1\n"); break;
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
        const int *full_mapping1 = &mapping_list[(size_t)r * n_padding];

        char title[64];
        snprintf(title, sizeof(title), "rmsd=%.6f", rmsd);

        if (savecoords_flag) {
            /* As in atormsd -savecoords: cluster 1, then cluster 2, in the
             * file PREFIX_<n>.xyz, with n counted from 1 */
            size_t len = strlen(pairs_prefix) + 32;
            char *pair_path = malloc(len);
            FILE *pair_file;
            if (!pair_path) { fprintf(stderr, "Error: out of memory\n"); goto done; }
            snprintf(pair_path, len, "%s_%d.xyz", pairs_prefix, r + 1);
            pair_file = fopen(pair_path, "w");
            if (!pair_file) {
                fprintf(stderr, "Can't open %s for writing\n", pair_path);
                free(pair_path);
                goto done;
            }
            write_xyz_frame(pair_file, title1, n_padding,
                            n_atoms1, atom_data1, coords1, NULL, NULL);
            write_xyz_frame(pair_file, title, n_padding,
                            n_atoms2, atom_data2, coords2, full_mapping1, transform);
            fclose(pair_file);
            free(pair_path);
        }

        if (coords_flag && !savecoords_flag) {
            /* As in atormsd -coords: written to stdout with or without
             * -align, replacing the RMSD lines */
            write_xyz_frame(stdout, title, n_padding,
                            n_atoms2, atom_data2, coords2, full_mapping1, transform);
        } else {
            printf("%.6f", rmsd);
            if (printmapping_flag) {
                /* 1-based, comma-separated, as printed by atormsd */
                putchar(' ');
                for (int i = 0; i < n_padding; i++)
                    printf(i ? ",%d" : "%d", full_mapping1[i] + 1);
            }
            putchar('\n');
        }
    }
    status = 0;

done:
    free(atom_data1); free(coords1);
    free(atom_data2); free(coords2);
    free(rmsd_list); free(mapping_list); free(transform_list);
    return status;
}
