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

#ifndef UTILS_H
#define UTILS_H

typedef struct Arguments {
    char infile[255];
    int n;
    int m;
    int k;
    int max_del;
    double alpha;
    char alpha_file[255];
    double beta;
    double gamma;
    char gamma_file[255];
    char mut_file[255];
    char cell_file[255];
    int print_leaves;
    int print_expected;
    int monoclonal;
    int repetitions;
    double start_temp;
    double cooling_rate;
    double el_a_variance;
    double el_b_variance;
    double el_g_variance;
} args_t;

void shuffle(int *array, int n);
// unsigned int get_randint();
char *remove_extension (char* path);
void import_input(int **input_matrix, int rows, int columns, char *path);
args_t* get_arguments(int cargc, char **cargsv);

#endif