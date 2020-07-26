#pragma clang diagnostic ignored "-Wpragma-once-outside-header"
#pragma once
#include <algorithm>
#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>

using namespace std;

class aco
{
public:
    // algo name
    string algoname;

    // basic parameter
    int runs;
    int generations;
    string input_filename;
    string output_filename;

    // aco parameter
    FILE *fp;
    int nodes;
    int ants;
    double **distance;
    double *walk_length;
    double min_length;
    double pheromone;
    double ***pheromone_table;
    int **place_table;
    int **solutions;
    int *min_solution;
    int *have_not_been;
    double *wheel;
    double alpha;
    double beta;
    double rho;
    double *convergence;
    double global_min_length;
    int *global_min_solution;
};