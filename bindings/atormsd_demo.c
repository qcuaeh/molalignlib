/* atormsd_demo.c
 *
 * Demo for atormsd_calculate (array-based interface).
 * Reads XYZ files and passes coordinate/element arrays to the Fortran library.
 *
 * Self-contained: nothing beyond molalign.h and libmolalign is needed to
 * build it.
 *
 * Compile:
 *   gcc atormsd_demo.c -o demo_atormsd -lm -lgfortran -lmolalign
 *
 * If molalign.h and libmolalign were installed to a user-local prefix (see
 * "User-local install" in the top-level README) rather than a system-wide
 * location, add explicit include and library search paths:
 *
 *   gcc atormsd_demo.c -o demo_atormsd \
 *       -I$HOME/.local/include/molalignlib -L$HOME/.local/lib \
 *       -lm -lgfortran -lmolalign
 *
 * Run:
 *   ./demo_atormsd atoms1.xyz atoms2.xyz [options]
 *
 * Pass -help for the full list of options.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include "molalign.h"
#include "error_codes.h"

/* ======================================================================
 * CLI - single-dash long-option parser
 * ====================================================================== */

typedef struct {
    const char *name;
    int has_arg;   /* 0 = no_argument, 1 = required_argument */
    int val;
} long_opt_t;

typedef struct {
    const char *desc;
    const char *arg;   /* argument metavar, or NULL */
} opt_info_t;

/* Parse the next option from argv, recognising options in any position
 * relative to non-option arguments.
 *
 * Non-option tokens are collected into posargs as they are encountered.
 * *npos counts every one of them (it must be 0 before the first call), but
 * only the first max_pos are stored, so the caller can detect too many
 * without posargs overflowing.
 *
 * Sets *arg to the option argument when has_arg=1.
 * Returns the matched val, 0 on error (unknown option or missing
 * argument), -1 when done. */
static int parse_long_opt(int argc, char **argv, int *argi, const char **arg,
                          const long_opt_t *opts,
                          const char **posargs, int max_pos, int *npos)
{
    /* Collect non-option tokens ("-" alone counts as one). */
    while (*argi < argc && (argv[*argi][0] != '-' || argv[*argi][1] == '\0')) {
        if (*npos < max_pos) posargs[*npos] = argv[*argi];
        (*npos)++;
        (*argi)++;
    }
    if (*argi >= argc) return -1;

    const char *p = argv[(*argi)++];
    for (int i = 0; opts[i].name; i++) {
        if (strcmp(p + 1, opts[i].name) != 0) continue;   /* exact match only */
        if (opts[i].has_arg) {
            if (*argi >= argc) {
                fprintf(stderr, "Error: -%s requires an argument\n", opts[i].name);
                return 0;
            }
            *arg = argv[(*argi)++];
        }
        return opts[i].val;
    }

    fprintf(stderr, "Error: unknown option: %s\n", p);
    return 0;
}

/* Print a formatted option list to stderr. */
static void print_options(const long_opt_t *opts, const opt_info_t *info)
{
    for (int i = 0; opts[i].name; i++) {
        const opt_info_t *oi = &info[opts[i].val];
        if (oi->arg)
            fprintf(stderr, "  -%-23s %s <%s>\n", opts[i].name, oi->desc, oi->arg);
        else
            fprintf(stderr, "  -%-23s %s\n", opts[i].name, oi->desc);
    }
}

/* ======================================================================
 * XYZ reader
 * ====================================================================== */

/* Lowercase atomic symbols indexed by atomic number. */
static const char *pt_symbols[] = {
    "",                                                                   /* 0      */
    "h",  "he",                                                           /* 1-2    */
    "li", "be", "b",  "c",  "n",  "o",  "f",  "ne",                     /* 3-10   */
    "na", "mg", "al", "si", "p",  "s",  "cl", "ar",                     /* 11-18  */
    "k",  "ca",                                                           /* 19-20  */
    "sc", "ti", "v",  "cr", "mn", "fe", "co", "ni", "cu", "zn",         /* 21-30  */
    "ga", "ge", "as", "se", "br", "kr",                                  /* 31-36  */
    "rb", "sr", "y",  "zr", "nb", "mo", "tc", "ru", "rh", "pd",         /* 37-46  */
    "ag", "cd", "in", "sn", "sb", "te", "i",  "xe",                     /* 47-54  */
    "cs", "ba",                                                           /* 55-56  */
    "la", "ce", "pr", "nd", "pm", "sm", "eu", "gd",                     /* 57-64  */
    "tb", "dy", "ho", "er", "tm", "yb", "lu",                           /* 65-71  */
    "hf", "ta", "w",  "re", "os", "ir", "pt", "au", "hg",              /* 72-80  */
    "tl", "pb", "bi", "po", "at", "rn",                                  /* 81-86  */
    "fr", "ra",                                                           /* 87-88  */
    "ac", "th", "pa", "u",  "np", "pu", "am", "cm",                     /* 89-96  */
    "bk", "cf", "es", "fm", "md", "no", "lr",                           /* 97-103 */
    "rf", "db", "sg", "bh", "hs", "mt", "ds", "rg", "cn",              /* 104-112*/
    "nh", "fl", "mc", "lv", "ts", "og"                                   /* 113-118*/
};
static const int pt_size = (int)(sizeof(pt_symbols) / sizeof(pt_symbols[0]));

/* Return atomic number for a lowercase element symbol (e.g. "c", "fe").
 * Returns 0 if not found. */
static int elnum_lookup(const char *elsym)
{
    for (int i = 1; i < pt_size; i++)
        if (strcmp(elsym, pt_symbols[i]) == 0) return i;
    return 0;
}

/* Split an atomic label into atomic number and group id.
 *
 * The label is matched case-insensitively. Its leading letters are the
 * element symbol and anything after them must be digits, giving the group
 * number (0 if absent), e.g. "C", "fe", "H12".
 *
 * Returns 0 on success, -1 on invalid label. */
static int parse_label(const char *sym, int *elnum, int *group)
{
    char label[32];
    int i, m = 0;

    for (i = 0; sym[i] && i < (int)sizeof(label) - 1; i++)
        label[i] = (char)tolower((unsigned char)sym[i]);
    label[i] = '\0';

    while (islower((unsigned char)label[m])) m++;
    for (i = m; label[i]; i++)
        if (!isdigit((unsigned char)label[i])) goto invalid;

    *group = atoi(&label[m]);   /* 0 when there is no digit suffix */
    label[m] = '\0';            /* keep just the element symbol */

    /* Also rejects an empty or over-long symbol: neither is in the table. */
    *elnum = elnum_lookup(label);
    if (*elnum == 0) goto invalid;
    return 0;

invalid:
    fprintf(stderr, "Invalid atomic label: %s\n", sym);
    return -1;
}

/* Check that path ends with ".xyz" (case-insensitive). */
static int has_xyz_extension(const char *path)
{
    size_t n = strlen(path);
    if (n < 4) return 0;
    for (int i = 0; i < 4; i++)
        if (tolower((unsigned char)path[n - 4 + i]) != ".xyz"[i]) return 0;
    return 1;
}

/* Read an XYZ file into packed atomdata and coords arrays.
 *
 * atomdata_out receives a heap-allocated array of length n_atoms*2:
 *   [elnum0, group0, elnum1, group1, ...]
 * coords_out receives a heap-allocated array of length n_atoms*3.
 * The caller is responsible for freeing both arrays.
 *
 * Returns 0 on success, non-zero on error (message written to stderr). */
static int read_xyz(const char *path,
                    int *n_atoms_out,
                    int **atomdata_out,
                    double **coords_out)
{
    FILE *fp;
    int n, elnum, group;
    char sym[32], line[256];
    int *atomdata = NULL;
    double *coords = NULL;

    if (!has_xyz_extension(path)) {
        fprintf(stderr, "Error: '%s' does not have a .xyz extension\n", path);
        return 1;
    }

    fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "Error: cannot open '%s'\n", path); return 1; }

    if (fscanf(fp, " %d", &n) != 1 || n <= 0) {
        fprintf(stderr, "Error: invalid atom count in '%s'\n", path);
        goto fail;
    }

    /* Skip the rest of the count line and the title line. */
    if (!fgets(line, sizeof(line), fp) || !fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "Error: missing title line in '%s'\n", path);
        goto fail;
    }

    atomdata = malloc((size_t)n * 2 * sizeof(int));
    coords   = malloc((size_t)n * 3 * sizeof(double));
    if (!atomdata || !coords) {
        fprintf(stderr, "Error: out of memory\n");
        goto fail;
    }

    for (int i = 0; i < n; i++) {
        double *xyz = &coords[i*3];
        if (fscanf(fp, " %31s %lf %lf %lf", sym, &xyz[0], &xyz[1], &xyz[2]) != 4) {
            fprintf(stderr, "Error: malformed atom line %d in '%s'\n", i+1, path);
            goto fail;
        }
        if (parse_label(sym, &elnum, &group) != 0) goto fail;
        atomdata[i*2]     = elnum;
        atomdata[i*2 + 1] = group;
    }

    fclose(fp);
    *n_atoms_out  = n;
    *atomdata_out = atomdata;
    *coords_out   = coords;
    return 0;

fail:
    free(atomdata);
    free(coords);
    fclose(fp);
    return 1;
}

/* ======================================================================
 * Demo program
 * ====================================================================== */

/* option value constants */
enum {
    OPT_ALIGN = 1, OPT_REMAP, OPT_HEAVY, OPT_MASS,
    OPT_MIRROR, OPT_LABEL, OPT_ASSIGNMENT, OPT_PRINTTRANS,
    OPT_STATS, OPT_RANDOM, OPT_PRUNE, OPT_FREQ,
    OPT_TRIALS, OPT_RECORDS, OPT_HELP
};

static const long_opt_t long_options[] = {
    {"align",           0, OPT_ALIGN},
    {"remap",           0, OPT_REMAP},
    {"heavy",           0, OPT_HEAVY},
    {"mass",            0, OPT_MASS},
    {"mirror",          0, OPT_MIRROR},
    {"label",           0, OPT_LABEL},
    {"stats",           0, OPT_STATS},
    {"random",          0, OPT_RANDOM},
    {"prune",           1, OPT_PRUNE},
    {"freq",            1, OPT_FREQ},
    {"trials",          1, OPT_TRIALS},
    {"records",         1, OPT_RECORDS},
    {"assignment",      0, OPT_ASSIGNMENT},
    {"printtrans",      0, OPT_PRINTTRANS},
    {"help",            0, OPT_HELP},
    {NULL, 0, 0}
};

static const opt_info_t opt_info[] = {
    [OPT_ALIGN]           = {"Enable structural alignment",                        NULL   },
    [OPT_REMAP]           = {"Enable atom remapping to minimise RMSD",             NULL   },
    [OPT_HEAVY]           = {"Use only heavy (non-hydrogen) atoms",                NULL   },
    [OPT_MASS]            = {"Weight atoms by their atomic masses",                NULL   },
    [OPT_MIRROR]          = {"Mirror the second molecule",                         NULL   },
    [OPT_LABEL]           = {"Use atom labels for matching",                       NULL   },
    [OPT_STATS]           = {"Print optimisation statistics",                      NULL   },
    [OPT_RANDOM]          = {"Use random algorithm",                               NULL   },
    [OPT_PRUNE]           = {"Set pruning tolerance",                              "TOL"  },
    [OPT_FREQ]            = {"Set convergence frequency (default: 10)",            "N"    },
    [OPT_TRIALS]          = {"Set maximum number of trials (default: 10000)",      "N"    },
    [OPT_RECORDS]         = {"Return up to N ranked solutions (default: 1)",       "N"    },
    [OPT_ASSIGNMENT]      = {"Print the atom permutation",                         NULL   },
    [OPT_PRINTTRANS]      = {"Print the 4x4 homogeneous transformation matrix",    NULL   },
    [OPT_HELP]            = {"Show this help message",                             NULL   },
};

static void print_usage(const char *prog)
{
    fprintf(stderr, "Usage: %s [OPTIONS] file1.xyz file2.xyz\n\n"
                    "Calculate RMSD between two atom clusters (XYZ format).\n\n"
                    "Options:\n", prog);
    print_options(long_options, opt_info);
    fputc('\n', stderr);
}

int main(int argc, char **argv)
{
    const char *optarg = NULL;
    bool align_flag = false, remap_flag = false, heavy_flag = false, mass_flag = false;
    bool mirror_flag = false, label_flag = false;
    bool print_assignment = false, print_transform = false;
    bool print_stats = false, random_flag = false;
    bool prune_flag = false;
    double prune_tol = 0.0;
    int conv_freq = 10, max_trials = 10000;
    int n_records = 1;
    int argi = 1, opt;

    /* positional arguments collected during the parse loop */
    const char *posargs[2] = {NULL, NULL};
    int npos = 0;

    while ((opt = parse_long_opt(argc, argv, &argi, &optarg, long_options,
                                 posargs, 2, &npos)) != -1) {
        switch (opt) {
        case OPT_ALIGN:           align_flag    = true;         break;
        case OPT_REMAP:           remap_flag    = true;         break;
        case OPT_HEAVY:           heavy_flag    = true;         break;
        case OPT_MASS:            mass_flag     = true;         break;
        case OPT_MIRROR:          mirror_flag   = true;         break;
        case OPT_LABEL:           label_flag    = true;         break;
        case OPT_STATS:           print_stats   = true;         break;
        case OPT_RANDOM:          random_flag   = true;         break;
        case OPT_PRUNE:           prune_flag = true; prune_tol = atof(optarg); break;
        case OPT_FREQ:            conv_freq     = atoi(optarg); break;
        case OPT_TRIALS:          max_trials    = atoi(optarg); break;
        case OPT_RECORDS:         n_records     = atoi(optarg); break;
        case OPT_ASSIGNMENT:      print_assignment = true;      break;
        case OPT_PRINTTRANS:      print_transform = true;       break;
        case OPT_HELP: print_usage(argv[0]); return 0;
        default:
            print_usage(argv[0]); return 1;
        }
    }

    if (npos != 2) {
        fprintf(stderr, "Error: expected exactly two XYZ file arguments\n");
        print_usage(argv[0]); return 1;
    }
    if (n_records < 1) {
        fprintf(stderr, "Error: -records must be at least 1\n");
        print_usage(argv[0]); return 1;
    }

    /* Everything below is released at `done`; free(NULL) is a no-op, so
     * every exit path can jump there regardless of how far it got. */
    int status = 1;
    int n1 = 0, n2 = 0;
    int *atom_data1 = NULL, *atom_data2 = NULL;
    double *coords1 = NULL, *coords2 = NULL;
    double *rmsd_list = NULL, *transform_list = NULL;
    int *atomperm_list = NULL;
    int occ_records = 0, error_code = MOLALIGN_SUCCESS;

    if (read_xyz(posargs[0], &n1, &atom_data1, &coords1) != 0) goto done;
    if (read_xyz(posargs[1], &n2, &atom_data2, &coords2) != 0) goto done;

    /* Output buffers hold up to n_records candidate solutions. Each
     * permutation has exactly n1 entries (one per atom of molecule 1). */
    rmsd_list      = malloc((size_t)n_records * sizeof(double));
    atomperm_list  = malloc((size_t)n_records * (size_t)n1 * sizeof(int));
    transform_list = malloc((size_t)n_records * 16 * sizeof(double));
    if (!rmsd_list || !atomperm_list || !transform_list) {
        fprintf(stderr, "Error: out of memory\n");
        goto done;
    }

    atormsd_calculate(
        n1, atom_data1, coords1,
        n2, atom_data2, coords2,
        align_flag, remap_flag, heavy_flag, mass_flag,
        mirror_flag, label_flag,
        print_stats, random_flag,
        prune_flag, prune_tol, conv_freq, max_trials,
        n_records,
        rmsd_list, atomperm_list,
        transform_list, &occ_records, &error_code);

    if (error_code != MOLALIGN_SUCCESS) {
        switch (error_code) {
        case MOLALIGN_ERROR_NOT_ISOMERS:        fprintf(stderr, "Error: molecules are not isomers\n"); break;
        case MOLALIGN_ERROR_ATOM_TYPE_MISMATCH: fprintf(stderr, "Error: atom types do not match\n");   break;
        case MOLALIGN_ERROR_PRUNED_ASSIGNMENT_FAILED:  fprintf(stderr, "Error: assignment failed (pruning tolerance might be too tight)\n"); break;
        default: fprintf(stderr, "Error: error code %d\n", error_code); break;
        }
        status = error_code;
        goto done;
    }

    for (int r = 0; r < occ_records; r++) {
        const double *rec_transform = &transform_list[r * 16];
        const int *rec_perm = &atomperm_list[r * n1];

        printf("RMSD: %.6f\n", rmsd_list[r]);

        if (print_assignment) {
            printf("Mapping:");
            for (int i = 0; i < n1; i++) printf(" %d", rec_perm[i] + 1); /* 1-based */
            putchar('\n');
        }

        if (print_transform) {
            printf("Transform:\n");
            for (int i = 0; i < 4; i++)
                printf("  %12.6f %12.6f %12.6f %12.6f\n",
                       rec_transform[i*4], rec_transform[i*4+1],
                       rec_transform[i*4+2], rec_transform[i*4+3]);
        }

        if (r < occ_records - 1) putchar('\n');
    }
    status = 0;

done:
    free(atom_data1); free(coords1);
    free(atom_data2); free(coords2);
    free(rmsd_list); free(atomperm_list); free(transform_list);
    return status;
}