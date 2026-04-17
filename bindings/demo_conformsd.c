/* demo_conformsd.c
 *
 * Demo for conformsd_calculate (array-based interface).
 * Reads XYZ files and passes coordinate/element arrays to the Fortran library.
 * Because XYZ carries no bond table, bond_flag is set to true so that
 * connectivity is derived from atomic geometry inside the library.
 *
 * Compile:
 *   gcc demo_conformsd.c -o demo_conformsd -lm -lgfortran -lmolalignlib
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "demo_utils.h"
#include "conformsd.h"

/* option value constants */
enum {
    OPT_ALIGN = 1, OPT_REMAP, OPT_HEAVY, OPT_MASS,
    OPT_MIRROR, OPT_LABEL,
    OPT_PRINTMAP, OPT_PRINTTRANS,
    OPT_STATS, OPT_RANDOM,
    OPT_FREQ, OPT_TRIALS,
    OPT_HELP
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
    {"freq",            1, OPT_FREQ},
    {"trials",          1, OPT_TRIALS},
    {"printmap",        0, OPT_PRINTMAP},
    {"printtrans",      0, OPT_PRINTTRANS},
    {"help",            0, OPT_HELP},
    {NULL, 0, 0}
};

static const opt_info_t opt_info[] = {
    [OPT_ALIGN]           = {"Enable structural alignment",                        NULL },
    [OPT_REMAP]           = {"Enable atom remapping to minimise RMSD",             NULL },
    [OPT_HEAVY]           = {"Use only heavy (non-hydrogen) atoms",                NULL },
    [OPT_MASS]            = {"Weight atoms by their atomic masses",                NULL },
    [OPT_MIRROR]          = {"Mirror the second molecule",                         NULL },
    [OPT_LABEL]           = {"Use atom labels for matching",                       NULL },
    [OPT_STATS]           = {"Print optimisation statistics",                      NULL },
    [OPT_RANDOM]          = {"Use random algorithm",                               NULL },
    [OPT_FREQ]            = {"Set convergence frequency (default: 100)",           "N"  },
    [OPT_TRIALS]          = {"Set maximum number of trials (default: 10000)",      "N"  },
    [OPT_PRINTMAP]        = {"Print the atom permutation",                         NULL },
    [OPT_PRINTTRANS]      = {"Print the 4x4 homogeneous transformation matrix",    NULL },
    [OPT_HELP]            = {"Show this help message",                             NULL },
};

static void print_usage(const char *prog)
{
    fprintf(stderr, "Usage: %s [OPTIONS] file1.xyz file2.xyz\n\n"
                    "Calculate RMSD between two molecular conformers (XYZ format).\n"
                    "Connectivity is derived from geometry (equivalent to -bond).\n\n"
                    "Options:\n", prog);
    print_options(long_options, opt_info);
    fputc('\n', stderr);
}

int main(int argc, char **argv)
{
    const char *optarg = NULL;
    bool align_flag = false, remap_flag = false, heavy_flag = false, mass_flag = false;
    bool mirror_flag = false, label_flag = false;
    bool print_mapping = false, print_transform = false;
    bool stats_flag = false, random_flag = false;
    int conv_freq = 100, max_trials = 10000;
    int argi = 1, opt;

    /* positional arguments collected during the parse loop */
    const char *posargs[2] = {NULL, NULL};
    int npos = 0;

    /* molecule data */
    int n1 = 0, n2 = 0;
    int *atom_data1 = NULL, *atom_data2 = NULL;
    double *coords1 = NULL, *coords2 = NULL;

    /* outputs */
    double rmsd, transform[16];
    int natoms = 0, error_code = 0;
    int *atomperm = NULL;

    while ((opt = parse_long_opt(argc, argv, &argi, &optarg,
                                 long_options, posargs, &npos)) != -1) {
        switch (opt) {
        case OPT_ALIGN:           align_flag      = true;         break;
        case OPT_REMAP:           remap_flag      = true;         break;
        case OPT_HEAVY:           heavy_flag      = true;         break;
        case OPT_MASS:            mass_flag       = true;         break;
        case OPT_MIRROR:          mirror_flag     = true;         break;
        case OPT_LABEL:           label_flag      = true;         break;
        case OPT_STATS:           stats_flag      = true;         break;
        case OPT_RANDOM:          random_flag     = true;         break;
        case OPT_FREQ:           conv_freq     = atoi(optarg); break;
        case OPT_TRIALS:          max_trials      = atoi(optarg); break;
        case OPT_PRINTMAP:   print_mapping   = true;         break;
        case OPT_PRINTTRANS: print_transform = true;         break;
        case OPT_HELP: print_usage(argv[0]); return 0;
        default:
            print_usage(argv[0]); return 1;
        }
    }

    if (npos != 2) {
        fprintf(stderr, "Error: expected exactly two XYZ file arguments\n");
        print_usage(argv[0]); return 1;
    }
    const char *file1 = posargs[0];
    const char *file2 = posargs[1];

    if (read_xyz(file1, &n1, &atom_data1, &coords1) != 0) return 1;
    if (read_xyz(file2, &n2, &atom_data2, &coords2) != 0) {
        free(atom_data1); free(coords1); return 1;
    }

    /* allocate permutation array - worst case is the larger molecule */
    atomperm = malloc((n1 > n2 ? n1 : n2) * sizeof(int));
    if (!atomperm) { fprintf(stderr, "Error: out of memory\n"); return 1; }

    /* XYZ has no bond table; pass n_bonds=0 and bond_flag=true so the library
     * derives connectivity from atomic geometry. */
    conformsd_calculate(
        n1, atom_data1, coords1, /*n_bonds1=*/0, /*bond_data1=*/NULL,
        n2, atom_data2, coords2, /*n_bonds2=*/0, /*bond_data2=*/NULL,
        align_flag, remap_flag, heavy_flag, mass_flag,
        mirror_flag, label_flag, /*bond_flag=*/true,
        stats_flag, random_flag,
        conv_freq, max_trials,
        &rmsd, &natoms, atomperm,
        transform, &error_code);

    if (error_code != 0) {
        switch (error_code) {
        case 1: fprintf(stderr, "Error: molecules are not isomers\n");  break;
        case 2: fprintf(stderr, "Error: molecules have no bonds\n");    break;
        case 3: fprintf(stderr, "Error: atom types do not match\n");    break;
        default: fprintf(stderr, "Error: error code %d\n", error_code); break;
        }
        free(atom_data1); free(coords1);
        free(atom_data2); free(coords2);
        free(atomperm);
        return error_code;
    }

    printf("RMSD: %.6f\n", rmsd);

    if (print_mapping) {
        int i;
        printf("Mapping:");
        for (i = 0; i < natoms; i++) printf(" %d", atomperm[i] + 1); /* 1-based */
        putchar('\n');
    }

    if (print_transform) {
        int i;
        printf("Transform:\n");
        for (i = 0; i < 4; i++)
            printf("  %12.6f %12.6f %12.6f %12.6f\n",
                   transform[i*4], transform[i*4+1],
                   transform[i*4+2], transform[i*4+3]);
    }

    free(atom_data1); free(coords1);
    free(atom_data2); free(coords2);
    free(atomperm);
    return 0;
}
