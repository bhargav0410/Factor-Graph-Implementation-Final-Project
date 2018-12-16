#include <mpi.h>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <random>
#include "ldpc_bp_mpi.h"
#include "ldpc_bp.h"

using namespace std::chrono;

int main (int argc, char* argv[]) {

    //Initialzing MPI
    MPI_Init(NULL, NULL);
    //getting size and rank
    int gsize, grank;
    MPI_Comm_size(MPI_COMM_WORLD, &gsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &grank);

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
    printf("n: %d, m: %d, k: %d at processor %d\n", n, m, k, grank);

    ldpc_bp_mpi ldpc(grank, gsize);
    if (grank == 0) {
        std::cout << "Creating parity check matrix...\n";
        start = high_resolution_clock::now();
        ldpc.create_H_mat(n, m, k);
        finish = high_resolution_clock::now();
        timediff = duration_cast<duration<double>>(finish - start);
        serial_time = std::min(serial_time, timediff.count());
        std::cout << "Time: " << serial_time << std::endl;

        std::cout << "Anticipated Rate >= " << ldpc.getRate() << "\n";
        std::cout << "Creating generator matrix...\n";
        serial_time = 1e30;
        start = high_resolution_clock::now();
        ldpc.gen_mat_from_H_mat();
        finish = high_resolution_clock::now();
        timediff = duration_cast<duration<double>>(finish - start);
        serial_time = std::min(serial_time, timediff.count());
        std::cout << "Time: " << serial_time << std::endl;

        std::cout << "Converting to standard form...\n";
        serial_time = 1e30;
        start = high_resolution_clock::now();
        ldpc.standard_form();
        finish = high_resolution_clock::now();
        timediff = duration_cast<duration<double>>(finish - start);
        serial_time = std::min(serial_time, timediff.count());
        std::cout << "Time: " << serial_time << std::endl;

        ldpc.check_matrices();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    ldpc.bcast_H_mat(0);
    ldpc.bcast_G_mat(0);
    /*
    for (int proc = 0; proc < gsize; proc++) {
        if (proc == grank) {
            ldpc.print_matrices();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    */
    ldpc.H_mat_comp_form();
   // printf("Created H matrix in comp form ar proc %d\n", grank);
    ldpc.create_list_from_mat();
   // printf("Created H matrix in comp form and list from matrix at proc %d\n", grank);
    //Creating input random vector and encoding it
    std::vector<int> in(num_syms*ldpc.get_num_input_syms()), out;    
    if (grank == 0) {
        srand(time(NULL));
        for (int i = 0; i < in.size(); i++) {
            in[i] = rand()%2;
        }
        print_vector(in);
        ldpc.print_matrices();
        
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //printf("Created input vector...\n");
    MPI_Bcast((void *)&in[0], (int)in.size(), MPI_INT, 0, MPI_COMM_WORLD);
    //print_vector(in);
    MPI_Barrier(MPI_COMM_WORLD);
    ldpc.encode_using_G_mat_mpi(in, out);
    printf("Processor %d encoding done.\n", grank);
    if (ldpc.check_vector_mpi(out) != 0) {
        printf("Processor %d encoding incorrect.\n", grank);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (grank == 0) {
        std::cout << "Final Rate = " << ldpc.getGenMatRate() << "\n";
        print_vector(out);
    }
    
    

    MPI_Finalize();
}