#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "conformsd_c_binding.h"

#define MAX_ATOMS 10000

void print_usage(const char *prog_name) {
    printf("Usage: %s [OPTIONS] file1 file2\n", prog_name);
    printf("\nCalculate RMSD between two molecular conformers\n\n");
    printf("Options:\n");
    printf("  -align              Enable structural alignment\n");
    printf("  -remap              Enable atom remapping to minimize RMSD\n");
    printf("  -mapping          Print atom permutation order\n");
    printf("  -exhaustive         Use exhaustive search algorithm\n");
    printf("  -stochastic         Use stochastic algorithm (no adaptiveness)\n");
    printf("  -label              Use atom labels for matching\n");
    printf("  -heavy              Consider only heavy (non-hydrogen) atoms\n");
    printf("  -mass               Weight atoms by their masses\n");
    printf("  -mirror             Mirror the second molecule\n");
    printf("  -thres <n>          Set conformer threshold (default: 100)\n");
    printf("  -trials <n>         Set maximum number of trials (default: 10000)\n");
    printf("  -records <n>        Set number of records to keep (default: 1)\n");
    printf("  -tree               Use tree-based algorithm\n");
    printf("  -stats              Print optimization statistics\n");
    printf("  -random             Use random initialization\n");
    printf("  -bond               Use bond information for matching\n");
    printf("  -h, -help           Show this help message\n");
}

int parse_int_arg(const char *option, const char *arg, int *argc_ptr, char **argv[], int *i_ptr) {
    int value;
    
    if (*i_ptr + 1 >= *argc_ptr) {
        fprintf(stderr, "Error: Option %s requires an argument\n", option);
        return -1;
    }
    
    (*i_ptr)++;
    if (sscanf((*argv)[*i_ptr], "%d", &value) != 1) {
        fprintf(stderr, "Error: Option %s requires an integer argument\n", option);
        return -1;
    }
    
    return value;
}

void print_permutation(int *atomperm, int natoms) {
    int i;
    
    if (natoms > 0) {
        printf("%d", atomperm[0] + 1);  // First element, no leading comma
        for (i = 1; i < natoms; i++) {
            printf(",%d", atomperm[i] + 1);  // Subsequent elements with comma separator
        }
    }
}

int main(int argc, char **argv) {
    // Initialize parameters with defaults
    int align = 0, remap = 0, heavy = 0, mass = 0;
    int mirror = 0, label = 0, bond = 0, mapping = 0;
    int stochastic = 0, exhaustive = 0, stats = 0;
    int tree = 0, random_flag = 0;
    int confo_thres = 100;
    int max_trials = 10000;
    int num_records = 1;
    char *file1 = NULL, *file2 = NULL;
    double rmsd;
    int natoms;
    int *atomperm;
    int error_code;
    int pos_arg_count = 0;
    int i, temp_int;
    
    // Allocate atom permutation array
    atomperm = (int *)malloc(MAX_ATOMS * sizeof(int));
    if (atomperm == NULL) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return 1;
    }
    
    // Parse command-line arguments
    for (i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            // Process option
            if (strcmp(argv[i], "-align") == 0) {
                align = 1;
            } else if (strcmp(argv[i], "-remap") == 0) {
                remap = 1;
            } else if (strcmp(argv[i], "-mapping") == 0) {
                mapping = 1;
            } else if (strcmp(argv[i], "-exhaustive") == 0) {
                exhaustive = 1;
                stochastic = 0;
            } else if (strcmp(argv[i], "-stochastic") == 0) {
                stochastic = 1;
                exhaustive = 0;
            } else if (strcmp(argv[i], "-label") == 0) {
                label = 1;
            } else if (strcmp(argv[i], "-heavy") == 0) {
                heavy = 1;
            } else if (strcmp(argv[i], "-mass") == 0) {
                mass = 1;
            } else if (strcmp(argv[i], "-mirror") == 0) {
                mirror = 1;
            } else if (strcmp(argv[i], "-tree") == 0) {
                tree = 1;
            } else if (strcmp(argv[i], "-stats") == 0) {
                stats = 1;
            } else if (strcmp(argv[i], "-random") == 0) {
                random_flag = 1;
            } else if (strcmp(argv[i], "-bond") == 0) {
                bond = 1;
            } else if (strcmp(argv[i], "-thres") == 0) {
                temp_int = parse_int_arg(argv[i], argv[i+1], &argc, &argv, &i);
                if (temp_int < 0) {
                    free(atomperm);
                    return 1;
                }
                confo_thres = temp_int;
            } else if (strcmp(argv[i], "-trials") == 0) {
                temp_int = parse_int_arg(argv[i], argv[i+1], &argc, &argv, &i);
                if (temp_int < 0) {
                    free(atomperm);
                    return 1;
                }
                max_trials = temp_int;
            } else if (strcmp(argv[i], "-records") == 0) {
                temp_int = parse_int_arg(argv[i], argv[i+1], &argc, &argv, &i);
                if (temp_int < 0) {
                    free(atomperm);
                    return 1;
                }
                num_records = temp_int;
            } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-help") == 0) {
                print_usage(argv[0]);
                free(atomperm);
                return 0;
            } else {
                fprintf(stderr, "Error: Unknown option: %s\n", argv[i]);
                print_usage(argv[0]);
                free(atomperm);
                return 1;
            }
        } else {
            // Positional argument
            if (pos_arg_count == 0) {
                file1 = argv[i];
                pos_arg_count++;
            } else if (pos_arg_count == 1) {
                file2 = argv[i];
                pos_arg_count++;
            } else {
                fprintf(stderr, "Error: Too many positional arguments\n");
                print_usage(argv[0]);
                free(atomperm);
                return 1;
            }
        }
    }
    
    // Validate required arguments
    if (pos_arg_count == 0) {
        fprintf(stderr, "Error: File paths are missing\n");
        print_usage(argv[0]);
        free(atomperm);
        return 1;
    } else if (pos_arg_count == 1) {
        fprintf(stderr, "Error: Too few file paths\n");
        print_usage(argv[0]);
        free(atomperm);
        return 1;
    }
    
    // Call Fortran wrapper
    conformsd_calculate(file1, file2,
                       align, remap, heavy, mass, mirror, label, bond,
                       stochastic, exhaustive, stats,
                       confo_thres, max_trials, num_records,
                       &rmsd, &natoms, atomperm, &error_code);
    
    // Handle errors
    if (error_code != 0) {
        switch (error_code) {
            case 1:
                fprintf(stderr, "Error: These molecules are not isomers\n");
                break;
            case 2:
                fprintf(stderr, "Error: Molecules have no bonds!\n");
                break;
            case 3:
                fprintf(stderr, "Error: Atom types do not match\n");
                break;
            default:
                fprintf(stderr, "Error: Unknown error code %d\n", error_code);
        }
        free(atomperm);
        return error_code;
    }
    
    // Print result (format matches Fortran output)
    printf("%.6f", rmsd);
    
    if (mapping && remap) {
        printf(" ");
        print_permutation(atomperm, natoms);
    }
    
    printf("\n");
    
    // Clean up
    free(atomperm);
    
    return 0;
}
