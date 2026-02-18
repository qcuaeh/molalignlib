#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atormsd_c_binding.h"

int main(int argc, char *argv[]) {
    const char *file1 = NULL;
    const char *file2 = NULL;
    int align = 0;
    int remap = 0;
    int heavy = 0;
    int mass = 0;
    int mirror = 0;
    int label = 0;
    int mapping = 0;
    int stats = 0;
    int random = 0;
    double prune_tol = -1.0;  // Default: use -near (no pruning)
    int ato_thres = 10;
    int max_trials = 10000;
    int num_records = 1;

    double rmsd;
    int natoms;
    int *atomperm = NULL;
    int error_code;

    /* Parse command line arguments */
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <file1> <file2> [options]\n", argv[0]);
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -align      Enable alignment\n");
        fprintf(stderr, "  -remap      Enable atom remapping\n");
        fprintf(stderr, "  -heavy      Use only heavy atoms\n");
        fprintf(stderr, "  -mass       Use atomic masses as weights\n");
        fprintf(stderr, "  -mirror     Mirror second molecule\n");
        fprintf(stderr, "  -label      Use atom labels\n");
        fprintf(stderr, "  -mapping    Print atom mapping\n");
        fprintf(stderr, "  -stats      Print statistics\n");
        fprintf(stderr, "  -random     Use random algorithm\n");
        fprintf(stderr, "  -prune <tol> Set pruning tolerance\n");
        fprintf(stderr, "  -near       Disable pruning (default)\n");
        fprintf(stderr, "  -thres <n>  Set atom threshold (default: 10)\n");
        fprintf(stderr, "  -trials <n> Set max trials (default: 10000)\n");
        fprintf(stderr, "  -records <n> Set number of records (default: 1)\n");
        return 1;
    }

    file1 = argv[1];
    file2 = argv[2];

    /* Parse options */
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "-align") == 0) {
            align = 1;
        } else if (strcmp(argv[i], "-remap") == 0) {
            remap = 1;
        } else if (strcmp(argv[i], "-heavy") == 0) {
            heavy = 1;
        } else if (strcmp(argv[i], "-mass") == 0) {
            mass = 1;
        } else if (strcmp(argv[i], "-mirror") == 0) {
            mirror = 1;
        } else if (strcmp(argv[i], "-label") == 0) {
            label = 1;
        } else if (strcmp(argv[i], "-mapping") == 0) {
            mapping = 1;
        } else if (strcmp(argv[i], "-stats") == 0) {
            stats = 1;
        } else if (strcmp(argv[i], "-random") == 0) {
            random = 1;
        } else if (strcmp(argv[i], "-near") == 0) {
            prune_tol = -1.0;
        } else if (strcmp(argv[i], "-prune") == 0 && i + 1 < argc) {
            prune_tol = atof(argv[++i]);
        } else if (strcmp(argv[i], "-thres") == 0 && i + 1 < argc) {
            ato_thres = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-trials") == 0 && i + 1 < argc) {
            max_trials = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-records") == 0 && i + 1 < argc) {
            num_records = atoi(argv[++i]);
        }
    }

    /* Allocate space for atom permutation (assume max 1000 atoms) */
    atomperm = (int *)malloc(1000 * sizeof(int));
    if (atomperm == NULL) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return 1;
    }

    /* Call the wrapper function */
    atormsd_calculate(file1, file2,
                      align, remap, heavy, mass,
                      mirror, label, mapping,
                      stats, random,
                      prune_tol, ato_thres, max_trials, num_records,
                      &rmsd, &natoms, atomperm, &error_code);

    /* Check for errors */
    if (error_code != 0) {
        switch (error_code) {
            case 1:
                fprintf(stderr, "Error: Molecules are not isomers\n");
                break;
            case 2:
                fprintf(stderr, "Error: Atom types do not match\n");
                break;
            default:
                fprintf(stderr, "Error: Unknown error code %d\n", error_code);
                break;
        }
        free(atomperm);
        return error_code;
    }

    /* Print results */
    printf("RMSD: %.6f\n", rmsd);

    if (mapping) {
        printf("Atom mapping (%d atoms): ", natoms);
        for (int i = 0; i < natoms; i++) {
            printf("%d ", atomperm[i] + 1);  // Convert back to 1-based for display
        }
        printf("\n");
    }

    /* Clean up */
    free(atomperm);

    return 0;
}
