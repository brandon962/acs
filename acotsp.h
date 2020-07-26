#pragma clang diagnostic ignored "-Wpragma-once-outside-header"
#pragma once
#include "aco.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>

using namespace std;

class AcoTSP : public aco
{
public:
    double wheel_max = 0;
    double choice_rand;
    int choice;
    double dtemp;
    double spray;
    double tau, eta;
    int itemp, itemp1, itemp2;
    int count;
    int walk_rand;
    double avg_length;
    double ever_min_length;
    char *ctemp;
    int keep_fix;
    int min_solutione_run;
    int min_ant;
    double min_ant_length;
    int min_generation;

public:
    AcoTSP(int _r, int _g, int _n, int _a, double _p, double _alpha, double _beta, double _rho, string _ifilename)
    {
        runs = _r;
        iterations = _g;
        nodes = _n;
        ants = _a;
        pheromone = _p;
        alpha = _alpha;
        beta = _beta;
        rho = _rho;
        input_filename = _ifilename;
        min_length = __DBL_MAX__;

        srand(time(NULL));

        have_not_been = (int *)calloc(nodes + 1, sizeof(int));
        wheel = (double *)calloc(nodes + 1, sizeof(double));
        distance = (double **)calloc(nodes + 1, sizeof(double *));
        for (int i = 0; i < nodes + 1; i++)
            distance[i] = (double *)calloc(nodes + 1, sizeof(double));

        walk_length = (double *)calloc(ants, sizeof(double));
        for (int i = 0; i < ants; i++)
            walk_length[i] = 0.0;

        pheromone_table = (double ***)calloc(nodes + 1, sizeof(double **));
        for (int i = 0; i < nodes + 1; i++)
        {
            pheromone_table[i] = (double **)calloc(nodes + 1, sizeof(double *));
            for (int j = 0; j < nodes + 1; j++)
                pheromone_table[i][j] = (double *)calloc(nodes + 1, sizeof(double));
        }
        for (int i = 0; i < nodes + 1; i++)
            for (int j = 0; j < nodes + 1; j++)
                pheromone_table[i][j][0] = 1.0;

        place_table = (int **)calloc(nodes + 1, sizeof(int *));
        for (int i = 0; i < nodes + 1; i++)
            place_table[i] = (int *)calloc(2, sizeof(int));

        solutions = (int **)calloc(ants, sizeof(int *));
        for (int i = 0; i < ants; i++)
            solutions[i] = (int *)calloc(nodes + 2, sizeof(int *));

        min_solution = (int *)calloc(nodes + 2, sizeof(int));
        for (int i = 0; i < nodes + 2; i++)
            min_solution[i] = 1;

        global_min_solution = (int *)calloc(nodes + 2, sizeof(int));
        global_min_length = __DBL_MAX__;

        convergence = (double *)calloc(iterations, sizeof(double));

        // read in the place table
        fp = fopen(input_filename.c_str(), "r");
        for (int i = 1; i < nodes + 1; i++)
            fscanf(fp, "%d %d %d", &place_table[i][0], &place_table[i][0], &place_table[i][1]);

        fclose(fp);

        // calculate the distance of the nodes
        for (int i = 1; i < nodes + 1; i++)
            for (int j = i; j < nodes + 1; j++)
                distance[i][j] = pow((pow((place_table[i][0] - place_table[j][0]), 2) + pow((place_table[i][1] - place_table[j][1]), 2)), 0.5);

        for (int i = 1; i < nodes + 1; i++)
            for (int j = 0; j < i; j++)
                distance[i][j] = distance[j][i];

        ctemp = (char *)calloc(10, sizeof(char));
    }

    int run()
    {
        mssg();
        avg_length = 0.0;
        ever_min_length = __DBL_MAX__;
        for (int i = 0; i < runs; i++)
        {
            count = 0;
            //cout << "start " << i << " runs" << endl;
            min_length = __DBL_MAX__;

            initialize();
            probabilityUpdate();
            for (int j = 0; j < iterations; j++)
            {
                count++;
                solutionConstruction();
                checkDistance();
                changeNeighbor();
                checkDistance();
                min_ant = 1;
                min_ant_length = __DBL_MAX__;
                for (int k = 0; k < ants; k++)
                {
                    if (walk_length[k] < min_ant_length)
                    {
                        min_ant_length = walk_length[k];
                    }
                }
                checkDistance();
                if (walk_length[min_ant] < 500)
                {
                    keep_fix = fix2OPTAll(min_ant);
                    while (keep_fix)
                    {
                        keep_fix = fix2OPTAll(min_ant);
                    }
                }
                checkDistance();

                min_ant_length = walk_length[min_ant];

                checkDistance();
                for (int k = 0; k < ants; k++)
                {

                    if (walk_length[k] < 480)
                    {
                        keep_fix = fix2OPTAll(k);
                        while (keep_fix)
                        {
                            keep_fix = fix2OPTAll(k);
                        }
                    }

                    checkDistance();
                    if (walk_length[k] < min_length)
                    {
                        min_length = walk_length[k];
                        min_generation = j;
                        for (int l = 1; l < nodes + 2; l++)
                            min_solution[l] = solutions[k][l];

                        /*cout << "now generation : " << count << endl;
                        cout << "now min length : " << min_length << endl;
                        cout << "solution : " << endl;

                        for (int l = 1; l < nodes + 2; l++)
                        {
                            itemp = min_solution[l];
                            cout << min_solution[l] << " ";
                        }
                        cout << endl
                             << endl;*/
                    }
                }

                checkDistance();
                pheromoneUpdateACS();
                probabilityUpdate();
                checkDistance();
                convergence[j] += min_length;
            }
            cout << "the " << i << " runs best length is : " << min_length << ", at " << min_generation << " runs" << endl;
            avg_length += min_length;
            if (min_length < global_min_length)
            {
                global_min_length = min_length;
                min_solutione_run = i;
                for (int l = 1; l < nodes + 2; l++)
                    global_min_solution[l] = min_solution[l];
            }

            output_filename = "data/output";
            itoa(i, ctemp, 10);
            output_filename += ctemp;
            output_filename += ".txt";

            fp = fopen(output_filename.c_str(), "w");
            for (int j = 1; j < nodes + 2; j++)
            {
                fprintf(fp, "%d %d %d%s", min_solution[j], place_table[min_solution[j]][1], place_table[min_solution[j]][0], "\n");
            }

            fclose(fp);
        }
        cout << endl;
        cout << "avg length is : " << avg_length / runs << endl;
        cout << "best length is : " << global_min_length;
        cout << ", at : " << min_solutione_run << " run" << endl;

        string convergence_output_file = "data/convergence.txt";
        fp = fopen(convergence_output_file.c_str(), "w");
        for (int i = 0; i < iterations; i++)
        {
            convergence[i] /= runs;
            fprintf(fp, "%d %f%s", i * ants, convergence[i], "\n");
        }
        fclose(fp);
        string global_min_solution_output_file = "data/global.txt";
        fp = fopen(global_min_solution_output_file.c_str(), "w");

        for (int i = 1; i < nodes + 2; i++)
        {
            fprintf(fp, "%d %d %d%s", global_min_solution[i], place_table[global_min_solution[i]][1], place_table[global_min_solution[i]][0], "\n");
        }
        fclose(fp);

        string global_min_avglength_output_file = "data/minavglength.txt";
        fp = fopen(global_min_avglength_output_file.c_str(), "w");

        for (int i = 0; i < iterations; i++)
        {
            fprintf(fp, "%d %f%s", i * ants, avg_length / runs, "\n");
        }
        fclose(fp);
        return 0;
    }

private:
    void mssg()
    {
        cout << "ant colony system." << endl;
        cout << "runs : " << runs << endl;
        cout << "iterations : " << iterations << endl;
        cout << "nodes : " << nodes << endl;
        cout << "ants number : " << ants << endl;
        cout << "pheromone : " << pheromone << endl;
        cout << "alpha(pheromone ratio) : " << alpha << endl;
        cout << "beta(distance ratio) : " << beta << endl;
        cout << "rho(evaporation) : " << rho << endl;
        cout << "input filename : " << input_filename << endl;
        cout << endl;
    }

    void initialize()
    {
        min_length = __DBL_MAX__;
        int itemp = 0;
        int irandom;
        int spray;

        for (int i = 0; i < ants; i++)
        {
            for (int j = 1; j < nodes + 1; j++)
                solutions[i][j] = j;

            solutions[i][nodes + 1] = solutions[i][1];
            // swapping for random path
            for (int j = 2; j < nodes + 1; j++)
            {
                irandom = rand() % (nodes - 1) + 2;
                itemp = solutions[i][j];
                solutions[i][j] = solutions[i][irandom];
                solutions[i][irandom] = itemp;
            }
        }

        for (int i = 0; i < nodes + 1; i++)
            for (int j = 0; j < nodes + 1; j++)
                pheromone_table[i][j][0] = 0.00001;

        for (int i = 0; i < ants; i++)
        {
            for (int j = 1; j < nodes + 1; j++)
                walk_length[i] += distance[solutions[i][j]][solutions[i][j + 1]];

            spray = pheromone / walk_length[i];
            for (int j = 1; j < nodes + 1; j++)
            {
                pheromone_table[solutions[i][j]][solutions[i][j + 1]][0] += spray;
                pheromone_table[solutions[i][j + 1]][solutions[i][j]][0] += spray;
            }
        }
        probabilityUpdate();
    }

    void showAnts()
    {
        for (int i = 0; i < ants; i++)
        {
            cout << "ant" << i << " : " << endl;
            for (int j = 1; j < nodes + 2; j++)
                cout << solutions[i][j] << " ";

            cout << endl;
            cout << "path length : " << walk_length[i] << endl
                 << endl;
        }
    }

    void probabilityUpdate()
    {

        for (int i = 1; i < nodes + 1; i++)
            for (int j = 1; j < nodes + 1; j++)
            {
                tau = pheromone_table[i][j][0];
                eta = 1 / distance[i][j];
                pheromone_table[i][j][1] = pow(tau, alpha) * pow(eta, beta);
            }
    }

    void solutionConstruction()
    {

        for (int i = 0; i < ants; i++)
        {
            walk_length[i] = 0;
            wheel_max = 0;
            for (int j = 2; j < nodes + 1; j++)
            {
                solutions[i][j] = 0;
                have_not_been[j] = 1;
            }

            for (int j = 2; j < nodes; j++)
            {
                walkWheel(i, j);
            }

            for (int j = 2; j < nodes + 1; j++)
                if (have_not_been[j] == 1)
                {
                    solutions[i][nodes] = j;
                    walk_length[i] += distance[solutions[i][nodes - 1]][solutions[i][nodes]];
                    walk_length[i] += distance[solutions[i][nodes]][solutions[i][nodes + 1]];
                    break;
                }
        }
    }

    void walkWheel(int i, int j)
    {
        wheel_max = 0;
        for (int k = 2; k < nodes + 1; k++)
            wheel[k] = 0;

        for (int k = 2; k < nodes + 1; k++)
        {
            itemp = have_not_been[k];
            if (have_not_been[k] == 1)
            {
                wheel_max += pheromone_table[solutions[i][j - 1]][k][1];
                wheel[k] = wheel_max;
            }
        }
        for (int k = 2; k < nodes + 1; k++)
            wheel[k] /= wheel_max;

        choice_rand = rand() / (double)(RAND_MAX + 1);
        for (choice = 2; choice < nodes + 1; choice++)
            if (choice_rand < wheel[choice])
                break;

        if (choice > nodes)
            for (choice = 2; choice < nodes + 1; choice++)
                if (have_not_been[choice] == 1)
                    break;

        solutions[i][j] = choice;
        have_not_been[choice] = 0;
        walk_length[i] += distance[solutions[i][j - 1]][solutions[i][j]];
    }

    void pheromoneUpdate()
    {
        for (int i = 1; i < nodes + 1; i++)
            for (int j = 1; j < nodes + 1; j++)
            {
                pheromone_table[i][j][0] *= rho;
                dtemp = pheromone_table[i][j][0] *= rho;
            }
        for (int i = 0; i < ants; i++)
        {
            walk_length[i] = 0;
            for (int j = 1; j < nodes + 1; j++)
                walk_length[i] += distance[solutions[i][j]][solutions[i][j + 1]];

            spray = pheromone / (double)walk_length[i];
            for (int j = 1; j < nodes + 1; j++)
            {
                pheromone_table[solutions[i][j]][solutions[i][j + 1]][0] += (1 - rho) * spray;
                pheromone_table[solutions[i][j + 1]][solutions[i][j]][0] += (1 - rho) * spray;
                //pheromone_table[solutions[i][j]][solutions[i][j + 1]][0] += spray;
                //pheromone_table[solutions[i][j + 1]][solutions[i][j]][0] += spray;
            }
        }
    }

    void pheromoneUpdateACS()
    {

        /*int best_ant = 0;
        double best_ant_length = __DBL_MAX__;

        for (int i = 0; i < ants; i++)
        {
            if (walk_length[i] < best_ant_length)
            {
                best_ant_length = walk_length[i];
                best_ant = i;
            }
        }
*/
        for (int i = 1; i < nodes + 1; i++)
            for (int j = 1; j < nodes + 1; j++)
            {
                pheromone_table[i][j][0] *= rho;
                if (min_length < 450)
                    dtemp = pheromone_table[i][j][0];
                if (pheromone_table[i][j][0] < 0.000001)
                    pheromone_table[i][j][0] = 0.000001;
            }

        spray = pheromone / min_ant_length;
        for (int j = 1; j < nodes + 1; j++)
        {
            pheromone_table[solutions[min_ant][j]][solutions[min_ant][j + 1]][0] += (1 - rho) * spray;
            pheromone_table[solutions[min_ant][j + 1]][solutions[min_ant][j]][0] += (1 - rho) * spray;
            if (min_length < 450)
                dtemp = pheromone_table[solutions[min_ant][j + 1]][solutions[min_ant][j]][0];
        }
    }

    int fix2OPT()
    {
        double temp_length = min_length;
        int *temp_solution = (int *)calloc(nodes + 2, sizeof(int));
        for (int i = 2; i < nodes + 1; i++)
        {
            for (int j = i + 3; j < nodes + 1; j++)
            {
                temp_length = min_length;
                if (distance[min_solution[i]][min_solution[j - 1]] + distance[min_solution[i + 1]][min_solution[j]] > distance[min_solution[i]][min_solution[i + 1]] + distance[min_solution[j - 1]][min_solution[j]])
                {
                    for (int k = 1; k < nodes + 2; k++)
                        temp_solution[k] = min_solution[k];

                    temp_length = 0;

                    for (int k = 2; k < i; k++)
                        temp_solution[k] = min_solution[k];

                    for (int k = i; k < j + 1; k++)
                        temp_solution[k] = min_solution[i + j - k];

                    for (int k = j + 1; k < nodes + 1; k++)
                        temp_solution[k] = min_solution[k];

                    for (int k = 1; k < nodes + 1; k++)
                        temp_length += distance[temp_solution[k]][temp_solution[k + 1]];
                }
                if (temp_length < min_length)
                {
                    min_length = temp_length;

                    for (int k = 1; k < nodes + 2; k++)
                        min_solution[k] = temp_solution[k];
                    //cout << "fixed 2-opt" << endl;

                    /*cout << "now generation : " << count << endl;
                    cout << "now min length : " << min_length << endl;
                    cout << "solution : " << endl;
                    for (int l = 1; l < nodes + 2; l++)
                    {
                        itemp = min_solution[l];
                        cout << min_solution[l] << " ";
                    }
                    cout << endl
                         << endl;*/

                    return 1;
                }
            }
        }
        return 0;
    }

    int fix2OPTAll(int a)
    {
        double temp_length = walk_length[a];
        int *temp_solution = (int *)calloc(nodes + 2, sizeof(int));
        for (int i = 2; i < nodes + 1; i++)
        {
            for (int j = i + 3; j < nodes + 1; j++)
            {
                temp_length = walk_length[a];
                if (distance[solutions[a][i]][solutions[a][j - 1]] + distance[solutions[a][i + 1]][solutions[a][j]] > distance[solutions[a][i]][solutions[a][i + 1]] + distance[solutions[a][j - 1]][solutions[a][j]])
                {
                    for (int k = 1; k < nodes + 2; k++)
                        temp_solution[k] = solutions[a][k];

                    temp_length = 0;

                    for (int k = 2; k < i; k++)
                        temp_solution[k] = solutions[a][k];

                    for (int k = i; k < j + 1; k++)
                        temp_solution[k] = solutions[a][i + j - k];

                    for (int k = j + 1; k < nodes + 1; k++)
                        temp_solution[k] = solutions[a][k];

                    for (int k = 1; k < nodes + 1; k++)
                        temp_length += distance[temp_solution[k]][temp_solution[k + 1]];
                }
                if (temp_length < walk_length[a])
                {
                    walk_length[a] = temp_length;

                    for (int k = 1; k < nodes + 2; k++)
                        solutions[a][k] = temp_solution[k];
                    //cout << "fixed 2-opt" << endl;
                    //spray = pheromone / (double)walk_length[a];

                    /*for (int k = 1; k < nodes + 1; k++)
                        for (int l = 1; l < nodes + 1; l++)
                            pheromone_table[k][l][0] *= rho;

                    for (int k = 1; k < nodes + 1; k++)
                    {
                        pheromone_table[solutions[a][k]][solutions[a][k + 1]][0] += (1 - rho) * spray;
                        pheromone_table[solutions[a][k + 1]][solutions[a][k]][0] += (1 - rho) * spray;
                    }*/

                    //probabilityUpdate();

                    /*cout << "now generation : " << count << endl;
                    cout << "now min length : " << walk_length[a] << endl;
                    cout << "solution : " << endl;
                    for (int l = 1; l < nodes + 2; l++)
                    {
                        itemp = solutions[a][l];
                        cout << solutions[a][l] << " ";
                    }
                    cout << endl
                         << endl;
*/
                    return 1;
                }
            }
        }
        return 0;
    }

    void changeNeighbor()
    {
        int a, b, c, d;
        int keep = 1;
        for (int i = 0; i < ants; i++)
        {
            keep = 1;
            while (keep)
            {
                keep = 0;
                for (int j = 2; j < nodes; j++)
                {
                    a = distance[solutions[i][j]][solutions[i][j + 1]];
                    b = distance[solutions[i][j + 1]][solutions[i][j + 2]];
                    c = distance[solutions[i][j - 1]][solutions[i][j + 1]];
                    d = distance[solutions[i][j]][solutions[i][j + 2]];
                    if (c + d < a + b)
                    {
                        //cout << "change" << endl;
                        walk_length[i] -= ((a + b) - (c + d));
                        itemp = solutions[i][j];
                        solutions[i][j] = solutions[i][j + 1];
                        solutions[i][j + 1] = itemp;
                        //keep = 1;
                    }
                }
            }
        }
    }

    void checkDistance()
    {
        double length = 0;
        for (int i = 0; i < ants; i++)
        {
            length = 0;
            for (int j = 1; j < nodes + 1; j++)
            {
                length += distance[solutions[i][j]][solutions[i][j + 1]];
            }
            walk_length[i] = length;
        }
    }
};