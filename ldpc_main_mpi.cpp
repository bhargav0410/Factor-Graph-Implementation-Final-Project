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
    double serial_encode = 0, mpi_encode = 0, mpi_decode = 0, serial_decode = 0;
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
        ldpc.create_H_mat(n, m, k);
        std::cout << "Anticipated Rate >= " << ldpc.getRate() << "\n";
        std::cout << "Creating generator matrix...\n";
    }
    ldpc.bcast_H_mat(0);
    MPI_Barrier(MPI_COMM_WORLD);
    ldpc.gen_mat_from_H_mat_mpi();
    std::cout << "Converting to standard form...\n";
    ldpc.standard_form();
    ldpc.check_matrices();
    MPI_Barrier(MPI_COMM_WORLD);
    //ldpc.bcast_G_mat(0);
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
        //print_vector(in);
        //ldpc.print_matrices();
        
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //printf("Created input vector...\n");
    MPI_Bcast((void *)&in[0], (int)in.size(), MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    start = high_resolution_clock::now();
    ldpc.encode_using_G_mat_mpi(in, out);
    finish = high_resolution_clock::now();
    if (grank == 0) {
        printf("MPI encoding time with %d procs for %d vectors: %f secs\n", gsize, num_syms, duration_cast<duration<double>>(finish - start).count());
    }
    printf("Processor %d encoding done.\n", grank);
    if (ldpc.check_vector_mpi(out) != 0) {
        printf("Processor %d encoding incorrect.\n", grank);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (grank == 0) {
        std::cout << "Final Rate = " << ldpc.getGenMatRate() << "\n";
        //print_vector(out);
    }
    
    std::vector<float> awgn(out.size()), chan_in(out.size()), chan_out(out.size());
    if (grank == 0) {
        float std_dev = pow((float)10.0, -((float)snr/(float)10.0));
        std::cout << "Noise power: " << std_dev << "\n";
        std::default_random_engine generator;
        std::normal_distribution<float> distribution(0.0, std_dev);

        //Passing through AWGN channel
        for (int i = 0; i < awgn.size(); i++) {
            //Creating vector ready to be transmitted through channel
            chan_in[i] = 2*(float)out[i] - 1;
            //Creating noise for channel emulation
            awgn[i] = distribution(generator);
            //Passing through AWGN channel
            chan_out[i] = chan_in[i] + awgn[i];
        }
        //print_vector(chan_in);
        //print_vector(chan_out);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast((void *)&chan_out[0], (int)chan_out.size(), MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<int> final_out;

    //Decode noise signal
    start = high_resolution_clock::now();
    ldpc.sum_product_decode_mpi_block(chan_out, final_out, iter, snr);
    finish = high_resolution_clock::now();
    if (grank == 0) {
        printf("MPI decoding time with %d procs for %d vectors: %f secs\n", gsize, num_syms, duration_cast<duration<double>>(finish - start).count());
    }

    if (grank == 0) {
        //print_vector(in);
        //print_vector(final_out);

        float ber = 0;
        for (int i = 0; i < final_out.size(); i++) {
            ber += abs(in[i] - final_out[i]);
        }
        std::cout << "BER: " << ber/(float)final_out.size() << "\n";
    }

    MPI_Finalize();
}