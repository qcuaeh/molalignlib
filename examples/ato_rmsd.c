/* ato_rmsd.c
 *
 * Minimal demo of ato_rmsd: reads two XYZ files, passes their atoms and
 * coordinates to the library, and prints the results. It shows how to pack
 * the input arrays, size the output buffers, interpret the returned
 * permutations, and apply the returned transforms.
 *
 * Compile:
 *   gcc ato_rmsd.c -o ato_rmsd \
 *       -I path/to/molalign.h/directory \
 *       -L path/to/libmolalign.a/directory \
 *       -lmolalign -lgfortran -lm
 *
 * Run:
 *   ./ato_rmsd file1.xyz file2.xyz [options]
 *
 * Run it without arguments for the list of options. Each option sets the
 * ato_rmsd argument of the same name (see molalign.h).
 *
 * Output: one line per record with the RMSD and the 1-based atom
 * permutation that maps cluster 2 onto cluster 1. With -coords, cluster 2,
 * reordered and transformed onto cluster 1, is written to stdout as one
 * XYZ frame per record instead.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include "molalign.h"

/* ---- Elements ---- */

/* Symbols indexed by atomic number, as in the library: elements 1-103,
 * then the Lennard-Jones pseudo-element "LJ" (104). Entry 0 is the symbol
 * written for padding atoms. */
static const char *symbols[] = {
    "X",
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "LJ"
};
#define N_ELEMENTS 104

/* Atomic number of an element symbol (case-insensitive), or -1 */
static int atomic_number(const char *symbol)
{
    for (int z = 1; z <= N_ELEMENTS; z++) {
        const char *s = symbols[z], *t = symbol;
        while (*s && *t && tolower((unsigned char)*s) == tolower((unsigned char)*t)) s++, t++;
        if (!*s && !*t) return z;
    }
    return -1;
}

/* ---- XYZ files ---- */

/* Read an XYZ file into the flat arrays ato_rmsd takes:
 *   atom_data: [elnum0, label0, elnum1, label1, ...]  (n_atoms*2)
 *   coords:    [x0, y0, z0, x1, ...]                  (n_atoms*3)
 * An atom name is an element symbol optionally followed by digits, the
 * label (e.g. "C2"; 0 without digits), used with -useatomtype.
 * Returns 0 on success; the caller frees both arrays. */
static int read_xyz(const char *path, int *n_atoms, int **atom_data, double **coords)
{
    FILE *fp = fopen(path, "r");
    char line[256], name[32];
    int n;

    if (!fp) { fprintf(stderr, "Cannot open %s\n", path); return 1; }
    /* Atom count, then the title line, which is not used */
    if (!fgets(line, sizeof line, fp) || sscanf(line, "%d", &n) != 1 || n < 1
            || !fgets(line, sizeof line, fp)) {
        fprintf(stderr, "Invalid XYZ file %s\n", path);
        fclose(fp);
        return 1;
    }

    *atom_data = malloc((size_t)n * 2 * sizeof(int));
    *coords    = malloc((size_t)n * 3 * sizeof(double));
    if (!*atom_data || !*coords) { fprintf(stderr, "Out of memory\n"); goto fail; }

    for (int i = 0; i < n; i++) {
        double *xyz = &(*coords)[i*3];
        size_t len;

        if (fscanf(fp, "%31s %lf %lf %lf", name, &xyz[0], &xyz[1], &xyz[2]) != 4) {
            fprintf(stderr, "Invalid atom line %d in %s\n", i + 1, path);
            goto fail;
        }
        len = strcspn(name, "0123456789");   /* symbol, then the label digits */
        (*atom_data)[i*2 + 1] = atoi(name + len);
        name[len] = '\0';
        (*atom_data)[i*2] = atomic_number(name);
        if ((*atom_data)[i*2] < 0) {
            fprintf(stderr, "Unknown element %s in %s\n", name, path);
            goto fail;
        }
    }

    fclose(fp);
    *n_atoms = n;
    return 0;

fail:
    free(*atom_data);
    free(*coords);
    fclose(fp);
    return 1;
}

/* Write one XYZ frame of cluster 2 reordered and transformed onto cluster 1.
 * Line i holds atom mapping[i] (0-based) of cluster 2; indices >= n_atoms
 * are padding atoms (only with -heavy), written as "X" at the origin. The
 * transform is a row-major 4x4 matrix: p_out = R*p_in + t. */
static void write_xyz(int n_padding, int n_atoms, const int *atom_data, const double *coords,
                      const int *mapping, const double *transform, double rmsd)
{
    printf("%d\nrmsd=%.6f\n", n_padding, rmsd);
    for (int i = 0; i < n_padding; i++) {
        int j = mapping[i];
        double p[3] = {0.0, 0.0, 0.0}, q[3];

        if (j < n_atoms) memcpy(p, &coords[j*3], sizeof(p));
        for (int k = 0; k < 3; k++)
            q[k] = transform[k*4]*p[0] + transform[k*4+1]*p[1] + transform[k*4+2]*p[2]
                 + transform[k*4+3];
        printf("%-3s  %12.6f  %12.6f  %12.6f\n",
               symbols[j < n_atoms ? atom_data[j*2] : 0], q[0], q[1], q[2]);
    }
}

/* ---- Program ---- */

static void print_usage(const char *prog)
{
    fprintf(stderr,
        "Usage: %s file1.xyz file2.xyz [options]\n\n"
        "Calculate the RMSD between two atom clusters with ato_rmsd.\n\n"
        "Options:\n"
        "  -align         Align cluster 2 onto cluster 1\n"
        "  -remap         Find the atom permutation that minimises the RMSD\n"
        "  -heavy         Ignore hydrogen atoms\n"
        "  -massweight    Use mass-weighted coordinates\n"
        "  -mirror        Reflect cluster 2 (x -> -x) before comparison\n"
        "  -useatomtype   Only match atoms with the same label (e.g. C1, C2)\n"
        "  -prune TOL     Prune atom pairings with tolerance TOL (Angstrom)\n"
        "  -minfreq N     Stop once the best solution is found N times (default: 10)\n"
        "  -maxtrials N   Stop after at most N random orientations (default: 10000)\n"
        "  -maxrecs N     Return up to N ranked solutions (default: 1)\n"
        "  -printstats    Print search statistics (with -align -remap)\n"
        "  -random        Seed the random-number generator from the clock\n"
        "  -coords        Write cluster 2, reordered and transformed, as XYZ\n",
        prog);
}

static const char *error_message(int error_code)
{
    switch (error_code) {
    case MOLALIGN_ERROR_INVALID_BOUND:          return "-minfreq, -maxtrials and -maxrecs must be at least 1";
    case MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER:  return "Atomic number out of range";
    case MOLALIGN_ERROR_NOT_ISOMERS:            return "The clusters have different atoms";
    case MOLALIGN_ERROR_ATOM_TYPE_MISMATCH:     return "The atoms are not in the same order (use -remap)";
    case MOLALIGN_ERROR_ASSIGNMENT_FAILED:      return "Assignment failed (pruning tolerance might be too tight)";
    default:                                    return "Unexpected error";
    }
}

int main(int argc, char **argv)
{
    /* ato_rmsd arguments, with the defaults of the atormsd program */
    bool align = false, remap = false, heavy = false, massweight = false;
    bool mirror = false, useatomtype = false, printstats = false, random = false;
    bool prune = false;
    double prune_tol = 0.0;
    int min_freq = 10, max_trials = 10000, max_recs = 1;
    bool coords_out = false;
    const char *paths[2];
    int n_paths = 0;

    if (argc < 2) { print_usage(argv[0]); return 0; }

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        bool has_value = i + 1 < argc;

        if      (!strcmp(a, "-align"))       align = true;
        else if (!strcmp(a, "-remap"))       remap = true;
        else if (!strcmp(a, "-heavy"))       heavy = true;
        else if (!strcmp(a, "-massweight"))  massweight = true;
        else if (!strcmp(a, "-mirror"))      mirror = true;
        else if (!strcmp(a, "-useatomtype")) useatomtype = true;
        else if (!strcmp(a, "-printstats"))  printstats = true;
        else if (!strcmp(a, "-random"))      random = true;
        else if (!strcmp(a, "-coords"))      coords_out = true;
        else if (!strcmp(a, "-prune") && has_value)     { prune = true; prune_tol = atof(argv[++i]); }
        else if (!strcmp(a, "-minfreq") && has_value)   min_freq = atoi(argv[++i]);
        else if (!strcmp(a, "-maxtrials") && has_value) max_trials = atoi(argv[++i]);
        else if (!strcmp(a, "-maxrecs") && has_value)   max_recs = atoi(argv[++i]);
        else if (a[0] != '-' && n_paths < 2)            paths[n_paths++] = a;
        else { fprintf(stderr, "Invalid argument: %s\n", a); return 1; }
    }
    if (n_paths != 2) { fprintf(stderr, "Two XYZ files are required\n"); return 1; }
    if (max_recs < 1) { fprintf(stderr, "%s\n", error_message(MOLALIGN_ERROR_INVALID_BOUND)); return 1; }

    int n_atoms1, n_atoms2, *atom_data1, *atom_data2;
    double *coords1, *coords2;
    if (read_xyz(paths[0], &n_atoms1, &atom_data1, &coords1)) return 1;
    if (read_xyz(paths[1], &n_atoms2, &atom_data2, &coords2)) return 1;

    /* Output buffers for up to max_recs records. Each permutation has
     * n_padding = max(n_atoms1, n_atoms2) entries: with -heavy the clusters
     * may have different numbers of hydrogens, and the smaller one is then
     * padded with padding atoms. */
    int n_padding = n_atoms1 > n_atoms2 ? n_atoms1 : n_atoms2;
    double *rmsd_list      = malloc((size_t)max_recs * sizeof(double));
    int    *mapping_list   = malloc((size_t)max_recs * n_padding * sizeof(int));
    double *transform_list = malloc((size_t)max_recs * 16 * sizeof(double));
    int n_records, error_code;
    if (!rmsd_list || !mapping_list || !transform_list) { fprintf(stderr, "Out of memory\n"); return 1; }

    fflush(stdout);   /* -printstats output is written by the library */

    ato_rmsd(
        n_atoms1, atom_data1, coords1,
        n_atoms2, atom_data2, coords2,
        align, remap, heavy, massweight,
        mirror, useatomtype,
        printstats, random,
        prune, prune_tol, min_freq, max_trials,
        max_recs,
        rmsd_list, mapping_list,
        transform_list, &n_records, &error_code);

    if (error_code != MOLALIGN_SUCCESS) {
        fprintf(stderr, "Error: %s\n", error_message(error_code));
        return 1;
    }

    /* Record r occupies mapping_list[r*n_padding ...] (0-based: entry i is
     * the atom of cluster 2 placed on line i of cluster 1) and
     * transform_list[r*16 ...] */
    for (int r = 0; r < n_records; r++) {
        const int *mapping = &mapping_list[r * n_padding];
        const double *transform = &transform_list[r * 16];

        if (coords_out) {
            write_xyz(n_padding, n_atoms2, atom_data2, coords2, mapping, transform, rmsd_list[r]);
        } else {
            printf("%.6f", rmsd_list[r]);
            for (int i = 0; i < n_padding; i++)
                printf(i ? ",%d" : " %d", mapping[i] + 1);
            putchar('\n');
        }
    }

    free(atom_data1); free(coords1);
    free(atom_data2); free(coords2);
    free(rmsd_list); free(mapping_list); free(transform_list);
    return 0;
}
