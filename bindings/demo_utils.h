#ifndef DEMO_UTILS_H
#define DEMO_UTILS_H

/*
 * demo_utils.h - self-contained utilities for the demo programs.
 *
 * Provides:
 *
 *   CLI option parsing:
 *     long_opt_t, opt_info_t
 *     int  parse_long_opt(argc, argv, *argi, **arg, *opts, *posargs, *npos)
 *     void print_options(const long_opt_t *opts, const opt_info_t *info)
 *
 *   XYZ file reading:
 *     int read_xyz(path, *n_atoms_out, **atomdata_out, **coords_out)
 *
 *   atomdata is a packed flat array of length n_atoms*2:
 *     [elnum0, group0, elnum1, group1, ...]
 *     group = 0 means unlabelled.
 *
 * Include this header in exactly one translation unit per binary.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

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

/* Parse options from argv, recognising them in any position relative to
 * non-option arguments.
 *
 * Non-option tokens are collected into posargs[0..(*npos)-1] as they are
 * encountered; posargs must point to a caller-allocated array of sufficient
 * size and *npos must be initialised to 0 before the first call.
 *
 * Sets *arg to the option argument when has_arg=1 and advances *argi.
 * Returns the matched val, 0 for unknown option, -1 when done. */
static int parse_long_opt(int argc, char **argv, int *argi, const char **arg,
                          const long_opt_t *opts,
                          const char **posargs, int *npos)
{
    int i;
    const char *p, *name;
    size_t nlen;

    for (;;) {
        if (*argi >= argc) return -1;
        p = argv[*argi];

        /* Non-option token: collect as positional argument and keep scanning. */
        if (p[0] != '-' || p[1] == '\0') {
            posargs[(*npos)++] = p;
            (*argi)++;
            continue;
        }
        break;
    }

    name = p + 1;

    for (i = 0; opts[i].name; i++) {
        nlen = strlen(opts[i].name);
        if (strncmp(name, opts[i].name, nlen) != 0) continue;
        if (name[nlen] != '\0') continue;   /* exact match only */

        (*argi)++;
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
    (*argi)++;   /* advance past unknown option to avoid an infinite loop */
    return 0;
}

/* Print a formatted option list to stderr. */
static void print_options(const long_opt_t *opts, const opt_info_t *info)
{
    int i;
    for (i = 0; opts[i].name; i++) {
        int v = opts[i].val;
        if (info[v].arg)
            fprintf(stderr, "  -%-23s %s <%s>\n",
                    opts[i].name, info[v].desc, info[v].arg);
        else
            fprintf(stderr, "  -%-23s %s\n",
                    opts[i].name, info[v].desc);
    }
}

/* ======================================================================
 * XYZ reader
 * ====================================================================== */

/* -----------------------------------------------------------------------
 * Periodic table - lowercase atomic symbols indexed by atomic number
 * ----------------------------------------------------------------------- */
static const char *_pt_symbols[] = {
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
static const int _pt_size =
    (int)(sizeof(_pt_symbols) / sizeof(_pt_symbols[0]));

/* Return atomic number for a lowercase element symbol (e.g. "c", "fe").
 * Returns 0 if not found. */
static int _elnum_lookup(const char *elsym)
{
    int i;
    for (i = 1; i < _pt_size; i++)
        if (strcmp(elsym, _pt_symbols[i]) == 0) return i;
    return 0;
}

/* Split an atomic label into atomic number and group id.
 *
 * The label is first normalised to lowercase. The leading alphabetic run
 * becomes the element symbol and the trailing digit run becomes the group
 * number (0 if absent). Any other character layout is an error.
 *
 * Returns 0 on success, -1 on invalid token. */
static int _parse_label(const char *sym, int *elnum, int *group)
{
    char normalized_label[32];
    char elsym[32];
    const char *p;
    int i, m;

    /* Normalise to lowercase */
    for (i = 0; sym[i] && i < (int)(sizeof(normalized_label) - 1); i++)
        normalized_label[i] = tolower((unsigned char)sym[i]);
    normalized_label[i] = '\0';

    /* Find the length of the leading alphabetic run */
    p = normalized_label;
    m = 0;
    while (p[m] && islower((unsigned char)p[m])) m++;

    /* Validate suffix: must be all digits (or empty) */
    for (i = m; normalized_label[i]; i++) {
        if (!isdigit((unsigned char)normalized_label[i])) {
            fprintf(stderr, "Invalid atomic label: %s\n", normalized_label);
            return -1;
        }
    }

    /* Validate alpha prefix length: 1-3 characters */
    if (m == 0 || m > 3) {
        fprintf(stderr, "Invalid atomic label: %s\n", normalized_label);
        return -1;
    }

    /* Extract element symbol and group number */
    strncpy(elsym, normalized_label, m);
    elsym[m] = '\0';

    *group = (normalized_label[m] != '\0') ? atoi(&normalized_label[m]) : 0;

    /* Lookup atomic number */
    *elnum = _elnum_lookup(elsym);
    if (*elnum == 0) {
        fprintf(stderr, "Invalid atomic label: %s\n", normalized_label);
        return -1;
    }

    return 0;
}

/* Check that path ends with ".xyz" (case-insensitive). */
static int _check_xyz_extension(const char *path)
{
    size_t n = strlen(path);
    if (n < 4) return 0;
    const char *ext = path + n - 4;
    return (tolower((unsigned char)ext[0]) == '.' &&
            tolower((unsigned char)ext[1]) == 'x' &&
            tolower((unsigned char)ext[2]) == 'y' &&
            tolower((unsigned char)ext[3]) == 'z');
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
    int n, i, elnum, group;
    char sym[32], line[256];
    double x, y, z;
    int *atomdata = NULL;
    double *coords = NULL;

    if (!_check_xyz_extension(path)) {
        fprintf(stderr, "Error: '%s' does not have a .xyz extension\n", path);
        return 1;
    }

    fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "Error: cannot open '%s'\n", path); return 1; }

    if (fscanf(fp, " %d", &n) != 1 || n <= 0) {
        fprintf(stderr, "Error: invalid atom count in '%s'\n", path);
        fclose(fp); return 1;
    }

    fgets(line, sizeof(line), fp);   /* rest of count line */
    fgets(line, sizeof(line), fp);   /* title line         */

    atomdata = malloc(n * 2 * sizeof(int));
    coords   = malloc(n * 3 * sizeof(double));
    if (!atomdata || !coords) {
        fprintf(stderr, "Error: out of memory\n");
        free(atomdata); free(coords);
        fclose(fp); return 1;
    }

    for (i = 0; i < n; i++) {
        if (fscanf(fp, " %31s %lf %lf %lf", sym, &x, &y, &z) != 4) {
            fprintf(stderr, "Error: malformed atom line %d in '%s'\n", i+1, path);
            free(atomdata); free(coords);
            fclose(fp); return 1;
        }
        if (_parse_label(sym, &elnum, &group) != 0) {
            free(atomdata); free(coords);
            fclose(fp); return 1;
        }
        atomdata[i*2]     = elnum;
        atomdata[i*2 + 1] = group;
        coords[i*3]       = x;
        coords[i*3 + 1]   = y;
        coords[i*3 + 2]   = z;
    }

    fclose(fp);
    *n_atoms_out  = n;
    *atomdata_out = atomdata;
    *coords_out   = coords;
    return 0;
}

#endif /* DEMO_UTILS_H */
