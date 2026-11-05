/* confo_rmsd.c
 *
 * Minimal demo of confo_rmsd: reads two MOL/SDF files (V2000, first record),
 * passes their atoms, coordinates and bonds to the library, and prints the
 * results. It shows how to pack the input arrays (bonds included, with bond
 * type codes), size the output buffers, interpret the returned permutations,
 * and apply the returned transforms.
 *
 * Compile:
 *   gcc confo_rmsd.c -o confo_rmsd \
 *       -I path/to/molalign.h/directory \
 *       -L path/to/libmolalign.a/directory \
 *       -lmolalign -lgfortran -lm
 *
 * Run:
 *   ./confo_rmsd file1.sdf file2.sdf [options]
 *
 * Run it without arguments for the list of options. Each option sets the
 * confo_rmsd argument of the same name (see molalign.h).
 *
 * Output: one line per record with the RMSD and the 1-based atom
 * permutation that maps molecule 2 onto molecule 1. With -coords, molecule
 * 2, reordered and transformed onto molecule 1, is written to stdout as one
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

/* ---- Bond types ---- */

/* Bond type codes (see "Bond types" in molalign.h): a single digit d has
 * code d, and a two-letter label "ab" has code 10 + 36*a + b, with
 * a = 0..25 for a..z and b = 0..35 for 0..9, a..z. */
#define BOND_SINGLE    1
#define BOND_DOUBLE    2
#define BOND_TRIPLE    3
#define BOND_AROMATIC  (10 + 36*('a' - 'a') + 10 + ('r' - 'a'))   /* "ar" = 37 */
#define BOND_UNDEFINED 1446

/* Code of a V2000 bond type: 1-3 single to triple, 4 aromatic; the query
 * types (5-8) have no fixed meaning, so they are of undefined type, which
 * is only an error with -usebondtype. */
static int bond_type_code(int v2000_type)
{
    switch (v2000_type) {
    case 1:  return BOND_SINGLE;
    case 2:  return BOND_DOUBLE;
    case 3:  return BOND_TRIPLE;
    case 4:  return BOND_AROMATIC;
    default: return BOND_UNDEFINED;
    }
}

/* ---- MOL/SDF files ---- */

/* Read the first record of a V2000 MOL/SDF file into the flat arrays
 * confo_rmsd takes:
 *   atom_data: [elnum0, label0, elnum1, label1, ...]  (n_atoms*2, labels 0)
 *   coords:    [x0, y0, z0, x1, ...]                  (n_atoms*3)
 *   bond_data: [atom1_0, atom2_0, type_0, ...]        (n_bonds*3, 1-based)
 * Returns 0 on success; the caller frees the three arrays. */
static int read_sdf(const char *path, int *n_atoms, int **atom_data, double **coords,
                    int *n_bonds, int **bond_data)
{
    FILE *fp = fopen(path, "r");
    char line[256], symbol[4];
    int na, nb;

    *atom_data = NULL; *coords = NULL; *bond_data = NULL;
    if (!fp) { fprintf(stderr, "Cannot open %s\n", path); return 1; }

    /* Three header lines, then the counts line: atoms and bonds, 3 columns each */
    for (int i = 0; i < 4; i++)
        if (!fgets(line, sizeof line, fp)) goto invalid;
    if (strstr(line, "V3000") || sscanf(line, "%3d%3d", &na, &nb) != 2 || na < 1 || nb < 0)
        goto invalid;

    *atom_data = malloc((size_t)na * 2 * sizeof(int));
    *coords    = malloc((size_t)na * 3 * sizeof(double));
    *bond_data = malloc((size_t)(nb > 0 ? nb : 1) * 3 * sizeof(int));
    if (!*atom_data || !*coords || !*bond_data) { fprintf(stderr, "Out of memory\n"); goto fail; }

    /* Atom block: x, y, z, element symbol, ... */
    for (int i = 0; i < na; i++) {
        double *xyz = &(*coords)[i*3];
        if (!fgets(line, sizeof line, fp)
                || sscanf(line, "%lf %lf %lf %3s", &xyz[0], &xyz[1], &xyz[2], symbol) != 4)
            goto invalid;
        (*atom_data)[i*2]     = atomic_number(symbol);
        (*atom_data)[i*2 + 1] = 0;   /* unlabelled */
        if ((*atom_data)[i*2] < 0) {
            fprintf(stderr, "Unknown element %s in %s\n", symbol, path);
            goto fail;
        }
    }

    /* Bond block: first atom, second atom, type, 3 columns each */
    for (int i = 0; i < nb; i++) {
        int *bond = &(*bond_data)[i*3], type;
        if (!fgets(line, sizeof line, fp)
                || sscanf(line, "%3d%3d%3d", &bond[0], &bond[1], &type) != 3)
            goto invalid;
        bond[2] = bond_type_code(type);
    }

    fclose(fp);
    *n_atoms = na;
    *n_bonds = nb;
    return 0;

invalid:
    fprintf(stderr, "Invalid or unsupported (not V2000) MOL/SDF file %s\n", path);
fail:
    free(*atom_data); free(*coords); free(*bond_data);
    fclose(fp);
    return 1;
}

/* Write one XYZ frame of molecule 2 reordered and transformed onto
 * molecule 1. Line i holds atom mapping[i] (0-based) of molecule 2;
 * indices >= n_atoms are padding atoms (only with -heavy), written as "X"
 * at the origin. The transform is a row-major 4x4 matrix:
 * p_out = R*p_in + t. */
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
        "Usage: %s file1.sdf file2.sdf [options]\n\n"
        "Calculate the symmetry-corrected RMSD between two conformers with confo_rmsd.\n"
        "Files are MOL/SDF V2000; only the first record is read.\n\n"
        "Options:\n"
        "  -align            Align molecule 2 onto molecule 1\n"
        "  -remap            Find the atom permutation that minimises the RMSD\n"
        "  -heavy            Ignore hydrogen atoms\n"
        "  -massweight       Use mass-weighted coordinates\n"
        "  -mirror           Reflect molecule 2 (x -> -x) before comparison\n"
        "  -rebond           Derive the bonds from the geometry instead of the files\n"
        "  -bondtol TOL      Bond-detection tolerance of -rebond (default: %.1f Angstrom)\n"
        "  -usebondtype      Also match bond types (single, double, triple, aromatic)\n"
        "  -minfreq N        Convergence threshold of the search (default: 100)\n"
        "  -maxtrials N      Stop after at most N random orientations (default: 10000)\n"
        "  -maxfrags N       Maximum number of molecular fragments (default: 1)\n"
        "  -maxrecs N        Return up to N ranked solutions (default: 1)\n"
        "  -printstats       Print search statistics (with -align -remap)\n"
        "  -printassigntree  Print the assignment tree (with -remap)\n"
        "  -random           Seed the random-number generator from the clock\n"
        "  -coords           Write molecule 2, reordered and transformed, as XYZ\n",
        prog, MOLALIGN_BOND_TOL_DEFAULT);
}

static const char *error_message(int error_code)
{
    switch (error_code) {
    case MOLALIGN_ERROR_INVALID_BOUND:          return "-minfreq, -maxtrials, -maxfrags and -maxrecs must be at least 1";
    case MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER:  return "Atomic number out of range";
    case MOLALIGN_ERROR_NOT_ISOMERS:            return "The molecules have different atoms";
    case MOLALIGN_ERROR_ATOM_TYPE_MISMATCH:     return "The atoms are not in the same order (use -remap)";
    case MOLALIGN_ERROR_BOND_MISMATCH:          return "The bonds are not in the same order (use -remap)";
    case MOLALIGN_ERROR_NOT_CONFORMERS:         return "The molecules have different bond graphs";
    case MOLALIGN_ERROR_TOO_MANY_FRAGMENTS:     return "Too many molecular fragments (see -maxfrags)";
    case MOLALIGN_ERROR_INVALID_BOND_TYPE:      return "Invalid bond type";
    case MOLALIGN_ERROR_UNDEFINED_BOND_TYPE:    return "Bonds of undefined type cannot be compared (see -usebondtype)";
    default:                                    return "Unexpected error";
    }
}

int main(int argc, char **argv)
{
    /* confo_rmsd arguments, with the defaults of the conformsd program */
    bool align = false, remap = false, heavy = false, massweight = false;
    bool mirror = false, rebond = false, usebondtype = false;
    bool printstats = false, printassigntree = false, random = false;
    double bond_tol = MOLALIGN_BOND_TOL_DEFAULT;
    int min_freq = 100, max_trials = 10000, max_frags = 1, max_recs = 1;
    bool coords_out = false;
    const char *paths[2];
    int n_paths = 0;

    if (argc < 2) { print_usage(argv[0]); return 0; }

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        bool has_value = i + 1 < argc;

        if      (!strcmp(a, "-align"))           align = true;
        else if (!strcmp(a, "-remap"))           remap = true;
        else if (!strcmp(a, "-heavy"))           heavy = true;
        else if (!strcmp(a, "-massweight"))      massweight = true;
        else if (!strcmp(a, "-mirror"))          mirror = true;
        else if (!strcmp(a, "-rebond"))          rebond = true;
        else if (!strcmp(a, "-usebondtype"))     usebondtype = true;
        else if (!strcmp(a, "-printstats"))      printstats = true;
        else if (!strcmp(a, "-printassigntree")) printassigntree = true;
        else if (!strcmp(a, "-random"))          random = true;
        else if (!strcmp(a, "-coords"))          coords_out = true;
        else if (!strcmp(a, "-bondtol") && has_value)   bond_tol = atof(argv[++i]);
        else if (!strcmp(a, "-minfreq") && has_value)   min_freq = atoi(argv[++i]);
        else if (!strcmp(a, "-maxtrials") && has_value) max_trials = atoi(argv[++i]);
        else if (!strcmp(a, "-maxfrags") && has_value)  max_frags = atoi(argv[++i]);
        else if (!strcmp(a, "-maxrecs") && has_value)   max_recs = atoi(argv[++i]);
        else if (a[0] != '-' && n_paths < 2)            paths[n_paths++] = a;
        else { fprintf(stderr, "Invalid argument: %s\n", a); return 1; }
    }
    if (n_paths != 2) { fprintf(stderr, "Two MOL/SDF files are required\n"); return 1; }
    if (max_recs < 1) { fprintf(stderr, "%s\n", error_message(MOLALIGN_ERROR_INVALID_BOUND)); return 1; }

    int n_atoms1, n_atoms2, n_bonds1, n_bonds2;
    int *atom_data1, *atom_data2, *bond_data1, *bond_data2;
    double *coords1, *coords2;
    if (read_sdf(paths[0], &n_atoms1, &atom_data1, &coords1, &n_bonds1, &bond_data1)) return 1;
    if (read_sdf(paths[1], &n_atoms2, &atom_data2, &coords2, &n_bonds2, &bond_data2)) return 1;

    /* Output buffers for up to max_recs records. Each permutation has
     * n_padding = max(n_atoms1, n_atoms2) entries: with -heavy the molecules
     * may have different numbers of hydrogens, and the smaller one is then
     * padded with padding atoms. */
    int n_padding = n_atoms1 > n_atoms2 ? n_atoms1 : n_atoms2;
    double *rmsd_list      = malloc((size_t)max_recs * sizeof(double));
    int    *mapping_list   = malloc((size_t)max_recs * n_padding * sizeof(int));
    double *transform_list = malloc((size_t)max_recs * 16 * sizeof(double));
    int n_records, error_code;
    if (!rmsd_list || !mapping_list || !transform_list) { fprintf(stderr, "Out of memory\n"); return 1; }

    fflush(stdout);   /* -printstats and -printassigntree output is written by the library */

    confo_rmsd(
        n_atoms1, atom_data1, coords1, n_bonds1, bond_data1,
        n_atoms2, atom_data2, coords2, n_bonds2, bond_data2,
        align, remap, heavy, massweight,
        mirror, /*useatomtype_flag=*/false, rebond, bond_tol,
        usebondtype,
        printstats, printassigntree, random,
        min_freq, max_trials, max_frags,
        max_recs,
        rmsd_list, mapping_list,
        transform_list, &n_records, &error_code);

    if (error_code != MOLALIGN_SUCCESS) {
        fprintf(stderr, "Error: %s\n", error_message(error_code));
        return 1;
    }

    /* Record r occupies mapping_list[r*n_padding ...] (0-based: entry i is
     * the atom of molecule 2 placed on line i of molecule 1) and
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

    free(atom_data1); free(coords1); free(bond_data1);
    free(atom_data2); free(coords2); free(bond_data2);
    free(rmsd_list); free(mapping_list); free(transform_list);
    return 0;
}
