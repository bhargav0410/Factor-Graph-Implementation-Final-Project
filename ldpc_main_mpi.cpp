#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <random>
#include "ldpc_bp_mpi.h"

using namespace std::chrono;

int main (int argc, char* argv[]) {

    //Initialzing MPI
    MPI_Init(NULL, NULL);
    //getting size and rank
    //int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //To measure execution time
    double serial_time = 1e30;
    duration<double> timediff;
    high_resolution_clock::time_point start, finish;
    float snr = 10; //SNR in dB

    int n = 12, m = 1, k = 4, iter = 10, num_syms = std::max(1024,n);
    float rate = 0;
    if (argc >= 4) {
        n = atoi(argv[1]);
        m = atoi(argv[2]);
        k = atoi(argv[3]);
        if (argc >= 5) {
            snr = strtof(argv[4], NULL);
        }
        if (argc >= 6) {
            iter = atoi(argv[5]);
        }
        if (argc >= 7) {
            num_syms = atoi(argv[6]);
        }
    }
    //printf("n: %d, m: %d, k: %d at processor %d\n", n, m, k, rank);

    MPI_Finalize();
}