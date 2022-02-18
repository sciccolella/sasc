/*
    MIT License

    Copyright (c) 2017-2019 Simone Ciccolella

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

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
    char buff[columns*2 + 10];

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

    //TODO: update this
    
    printf("\vRequired arguments:\n");
    printf("\t-n CELLS\t\tNumber of cells in the input file.\n");
    printf("\t-m MUTATIONS\t\tNumber of mutations in the input file.\n");
    printf("\t-k DOLLOK\t\tK value of Dollo(k) model.\n");
    printf("\t-a ALPHA\t\tFalse Negative rate in the input file or path to FN rates file.\n");
    printf("\t-b BETA\t\t\tFalse Positive rate in the input file.\n");
    printf("\t-i INFILE\t\tPath of the input file.\n");

    printf("\vOptional arguments:\n");
    printf("\vModel parameters (optional):\n");
    printf("\t-d DELETIONS\t\tMaximum number of total deletions allowed in the solution (default: INT_MAX).\n");
    printf("\t-e MUTFILE\t\tPath of the file containing mutations' names.\n");
    printf("\t-E CELLFILE\t\tPath of the file containing cells' names.\n");
    printf("\t-g GAMMA\t\tLoss rate in the input file or path of the file containing different GAMMA rates for each mutations.\n");
    printf("\t-r REPETITIONS\t\tSet the total number of Simulated Annealing repetitions (default: 5).\n");
    printf("\t-M \t\t\tForce SASC to infer a monoclonal tree (default: false).\n");
    
    printf("\vOutput parameters (optional):\n");
    printf("\t-l \t\t\tOutput a mutational tree with cells attached to it. Otherwise cells will not be present.\n");
    printf("\t-x \t\t\tIf this option is use, SASC will also output the expected matrix E.\n");
    printf("\t-P \t\t\tPrefix for the output files.\n");

    printf("\vSimulated Annealing parameters (optional):\n");
    printf("\t-S STARTTEMP\t\tStarting temperature of the Simulated Annealing algorithm.\n");
    printf("\t-C COOLINGRATE\t\tCooling rate of the Simulated Annealing algorithm.\n");

    printf("\vError parameters (optional):\n");
    printf("\t-A ALPHASTDV\t\tStandard deviation for new FN discovery.\n");
    printf("\t-B BETASTDV\t\tStandard deviation for new FP discovery.\n");
    printf("\t-G GAMMASTDV\t\tStandard deviation for new GAMMA discovery.\n");

    exit(EXIT_SUCCESS);
}

args_t*
get_arguments(int cargc, char **cargsv) {
    char *argp_program_version = "SACS -- version: 3.0";

    args_t *arguments = malloc(sizeof(args_t));

    arguments->max_del = INT_MAX;
    arguments->print_leaves = 0;
    arguments->print_expected = 0;
    arguments->monoclonal = 0;
    arguments->repetitions = 5;
    arguments->start_temp = 10000.0;
    arguments->cooling_rate = 0.01;

    arguments->gamma = 1;

    arguments->el_a_variance = 0;
    arguments->el_b_variance = 0;
    arguments->el_g_variance = 0;

    arguments->cores = 1;

    strcpy(arguments->mut_file, "NULL");
    strcpy(arguments->cell_file, "NULL");
    int inserted_args = 0;

    char *cvalue = NULL;
    int c;
    int prefix = 0;

    double single_alpha;
    double single_gamma;
    int res;

    opterr = 0;

    while ((c = getopt(cargc, cargsv, "hVm:n:a:b:g:k:i:e:E:d:lxMr:S:C:A:B:G:p:P:")) != - 1) {
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
                res = sscanf(optarg, "%lf", &single_alpha);
                if (res == 1) {
                    arguments->alpha = single_alpha;
                } else {
                    arguments->alpha = -1;
                    strcpy(arguments->alpha_file, optarg);
                }
                inserted_args++;
                break;
            case 'b':
                sscanf(optarg, "%lf", &arguments->beta);
                inserted_args++;
                break;
            case 'g':
                res = sscanf(optarg, "%lf", &single_gamma);
                if (res == 1) {
                    arguments->gamma = single_gamma;
                } else {
                    arguments->gamma = -1;
                    strcpy(arguments->gamma_file, optarg);
                }
                break;
            case 'k':
                arguments->k = atoi(optarg);
                inserted_args++;
                break;
            case 'i':
                strcpy(arguments->infile, optarg);
                inserted_args++;
                break;
            case 'P':
                strcpy(arguments->prefix, optarg);
                prefix = 1;
                break;
            case 'e':
                strcpy(arguments->mut_file, optarg);
                break;
            case 'E':
                strcpy(arguments->cell_file, optarg);
                break;
            case 'h':
                print_help();
                break;
            case 'd':
                arguments->max_del = atoi(optarg);
                break;
            case 'l':
                arguments->print_leaves = 1;
                break;
            case 'x':
                arguments->print_expected = 1;
                break;
            case 'M':
                arguments->monoclonal = 1;
                break;
            case 'r':
                arguments->repetitions = atoi(optarg);
                break;
            case 'v':
                printf("%s\n", argp_program_version);
                exit(0);
            case 'S':
                sscanf(optarg, "%lf", &arguments->start_temp);
                break;
            case 'C':
                sscanf(optarg, "%lf", &arguments->cooling_rate);
                break;
            case 'A':
                sscanf(optarg, "%lf", &arguments->el_a_variance);
                break;
            case 'B':
                sscanf(optarg, "%lf", &arguments->el_b_variance);
                break;
            case 'G':
                sscanf(optarg, "%lf", &arguments->el_g_variance);
                break;
            case 'p':
                arguments->cores = atoi(optarg);
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

    if (arguments->k == 0)
        arguments->max_del = 0;
    if (arguments->max_del == 0)
        arguments->k = 0;

    if (prefix == 0) {
        strcpy(arguments->prefix, remove_extension(arguments->infile));
    }

    return arguments;
}
