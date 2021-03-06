#pragma clang diagnostic ignored "-Wunused-private-field"

#include "acotsp.h"
#include <ctime>
#include <iostream>
#include <string>

using namespace std;
/////
int main(int argc, char *argv[])
{
    // basic parameter
    int runs = 30;
    int iterations = 500;
    string input_filename = "eil51.txt";

    // aco parameter
    int nodes = 51;
    int ants = 51;
    double pheromone = 1;
    double alpha = 1;
    double beta = 3;
    double rho = 0.85;

    if (argc > 5)
    {
        runs = atoi(argv[1]);
        iterations = atoi(argv[2]);
        nodes = atoi(argv[3]);
        ants = atoi(argv[4]);
        pheromone = atof(argv[5]);
        alpha = atof(argv[6]);
        beta = atof(argv[7]);
        rho = atof(argv[8]);
    }
    time_t start, end;

    start = time(NULL);
    AcoTSP aco = AcoTSP(runs, iterations, nodes, ants, pheromone, alpha, beta, rho, input_filename);
    aco.run();
    end = time(NULL);
    cout << "cost : " << end - start << "s " << endl;
    //
    return 0;
}