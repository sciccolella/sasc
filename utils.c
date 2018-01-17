#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include "utils.h"
#include "mt19937ar.h"

void
shuffle(int *array, int n) {
    for (int i = 0; i < n - 1; i++) {
        size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
        int t = array[j];
        array[j] = array[i];
        array[i] = t;
    }
}

char
*remove_extension (char* path) {
    char *retstr, *lastdot, *lastsep;
    if (path == NULL)
        return NULL;
    if ((retstr = malloc (strlen (path) + 1)) == NULL)
        return NULL;

    strcpy (retstr, path);
    lastdot = strrchr (retstr, '.');
    lastsep = ('/' == 0) ? NULL : strrchr (retstr, '/');

    if (lastdot != NULL) {
        if (lastsep != NULL) {
            if (lastsep < lastdot) {
                *lastdot = '\0';
            }
        } else {
            *lastdot = '\0';
        }
    }

    return retstr;
}

void
import_input(int **input_matrix, int rows, int columns, char *path){
    FILE *fp;
    char buff[rows*2 + 10];

    char seps[] = " ,\t";
    char *token;

    fp = fopen(path, "r");
    int i = 0;
    if (fp == NULL) {
        fprintf(stderr, "ERROR: Cannot open INFILE.\n");
        exit(EXIT_FAILURE);
    }
        
    while (fgets(buff, sizeof buff, fp) != NULL) {
        int j = 0;

        buff[strcspn(buff, "\n")] = '\0';
        
        token = strtok(buff, seps);
        while (token != NULL) {
            input_matrix[i][j] = atoi(token);

            token = strtok(NULL, seps);
            j++;
        }
        
        

        i++;
    }
    fclose(fp);
}

void
print_help() {
    printf("SASC -- Simulated Anneling for cancer progression inference via Single Cell Sequencing.\n");
    
    printf("\vRequired arguments:\n");
    printf("\t-n CELLS\t\tNumber of cells in the input file.\n");
    printf("\t-m MUTATIONS\t\tNumber of mutations in the input file.\n");
    printf("\t-k DOLLOK\t\tK value of Dollo(k) model.\n");
    printf("\t-a ALPHA\t\tFalse Negative rate in the input file.\n");
    printf("\t-b BETA\t\t\tFalse Positive rate in the input file.\n");
    printf("\t-i INFILE\t\tPath of the input file.\n");

    printf("\vOptional arguments:\n");
    printf("\t-e MUTFILE\t\tPath of the file containing mutations' names.\n");
    exit(EXIT_SUCCESS);
}

args_t*
get_arguments(int cargc, char **cargsv) {
    char *argp_program_version = "SACS -- version: 0.1";

    args_t *arguments = malloc(sizeof(args_t));

    arguments->max_del = INT_MAX;

    strcpy(arguments->mut_file, "NULL");
    int inserted_args = 0;

    char *cvalue = NULL;
    int c;

    opterr = 0;

    while ((c = getopt(cargc, cargsv, "hVm:n:a:b:k:i:e:d:")) != - 1) {
        switch(c) {
            case 'm':
                arguments->m = atoi(optarg);
                inserted_args++;
                break;
            case 'n':
                arguments->n = atoi(optarg);
                inserted_args++;
                break;
            case 'a':
                sscanf(optarg, "%lf", &arguments->alpha);
                inserted_args++;
                break;
            case 'b':
                sscanf(optarg, "%lf", &arguments->beta);
                inserted_args++;
                break;
            case 'k':
                arguments->k = atoi(optarg);
                inserted_args++;
                break;
            case 'i':
                strcpy(arguments->infile, optarg);
                inserted_args++;
                break;
            case 'e':
                strcpy(arguments->mut_file, optarg);
                break;
            case 'h':
                print_help();
                break;
            case 'd':
                arguments->max_del = atoi(optarg);
                break;
            case 'V':
                printf("%s\n", argp_program_version);
                exit(0);
                break;
            case '?':
                if (isprint (optopt))
                    fprintf (stderr, "ERROR: Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr,
                            "ERROR: Unknown option character `\\x%x'.\n",
                            optopt);
            default:
                exit(EXIT_FAILURE);
        }
    }

    if (inserted_args != 6) {
        fprintf (stderr, "ERROR: Not all required arguments were passed.\n");
        exit(EXIT_FAILURE);
    }

    return arguments;
}
