#include <mpi.h>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <random>
#include <fstream>
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
    double serial_construct = 0, mpi_construct = 0, serial_encode = 0, mpi_encode = 0, mpi_decode = 0, serial_decode = 0;
    duration<double> timediff;
    high_resolution_clock::time_point start, finish;
    float snr = 10; //SNR in dB

    int n = 12, m = 1, k = 4, iter = 10, num_syms = 1, num_times = 1, qam_size = 2;
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
        if (argc >= 8) {
            num_times = atoi(argv[7]);
        }
        if (argc >= 9) {
            qam_size = atoi(argv[8]);
        }
    }
    printf("n: %d, m: %d, k: %d at processor %d\n", n, m, k, grank);

    for (int times = 0; times < num_times; times++) {
        ldpc_bp_mpi ldpc(grank, gsize);
        ldpc.set_contellation(qam_size);
        if (grank == 0) {
            std::cout << "Creating parity check matrix...\n";
            ldpc.create_H_mat(n, m, k);
        //    ldpc.H_mat_to_salt_form();
            std::cout << "Anticipated Rate >= " << ldpc.getRate() << "\n";
            std::cout << "Creating generator matrix...\n";
            
            start = high_resolution_clock::now();
            ldpc.gen_mat_from_H_mat();
            finish = high_resolution_clock::now();
            serial_construct += duration_cast<duration<double>>(finish - start).count();
            
        }
        MPI_Barrier(MPI_COMM_WORLD);
        ldpc.bcast_H_mat(0);
        ldpc.bcast_G_mat(0);
        start = high_resolution_clock::now();
        //ldpc.gen_mat_from_H_mat_mpi();
        finish = high_resolution_clock::now();
        mpi_construct += duration_cast<duration<double>>(finish - start).count();
        std::cout << "Converting to standard form...\n";
        ldpc.standard_form();
        ldpc.check_matrices();
        MPI_Barrier(MPI_COMM_WORLD);

        ldpc.H_mat_comp_form();
    // printf("Created H matrix in comp form ar proc %d\n", grank);
        ldpc.create_list_from_mat();
    // printf("Created H matrix in comp form and list from matrix at proc %d\n", grank);
        //Creating input random vector and encoding it
        std::vector<int> in(num_syms*ldpc.get_num_input_syms()), out;
        std::vector<std::complex<float>> qam_out;
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
        
        if (grank == 0) {
         //   printf("Serial encode...\n");
            printf("In size: %d\n", (int)in.size());
            start = high_resolution_clock::now();
            ldpc.encode_using_G_mat(in, out);
            finish = high_resolution_clock::now();
            serial_encode += duration_cast<duration<double>>(finish - start).count();
        }
        MPI_Barrier(MPI_COMM_WORLD);
      //  printf("MPI encode...\n");
        start = high_resolution_clock::now();
        ldpc.encode_using_G_mat_mpi(in, out);
        printf("Processor %d encoding done.\n", grank);
        ldpc.gray_to_qam_mpi(out, qam_out);
        finish = high_resolution_clock::now();
        mpi_encode += duration_cast<duration<double>>(finish - start).count();
        printf("Processor %d QAM encoding done.\n", grank);
        MPI_Barrier(MPI_COMM_WORLD);
        if (ldpc.check_vector_mpi(out) != 0) {
            printf("Processor %d encoding incorrect.\n", grank);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (grank == 0) {
            std::cout << "Final Rate = " << ldpc.getGenMatRate() << "\n";
            //print_vector(out);
        }
        
        std::vector<std::complex<float>> awgn(qam_out.size()), chan_in(qam_out.size()), chan_out(qam_out.size());
        float std_dev = (pow((float)10.0, -((float)snr/(float)10.0)));
        if (grank == 0) {
            std::cout << "Noise power: " << std_dev*std_dev << "\n";
            std::default_random_engine generator;
            std::normal_distribution<float> distribution(0.0, std_dev);

            //Passing through AWGN channel
            for (int i = 0; i < awgn.size(); i++) {
                //Creating vector ready to be transmitted through channel
               // chan_in[i] = 2*(float)out[i] - 1;
                //Creating noise for channel emulation
                awgn[i] = std::complex<float>(distribution(generator), distribution(generator));
                //Passing through AWGN channel
                chan_out[i] = qam_out[i] + awgn[i];
      //          awgn[i] = std::abs(awgn[i]);
            }
      //      print_vector(awgn);
      //      printf("Chan out size: %d\n", (int)chan_out.size());
      //      print_vector(chan_out);
            printf("Signal passed through channel...\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast((void *)&chan_out[0], (int)chan_out.size(), MPI_COMPLEX, 0, MPI_COMM_WORLD);
        std::vector<float> llr_out;
        std::vector<int> final_out;

        //Decode noise signal
        /*
        if (grank == 0) {
            start = high_resolution_clock::now();
            ldpc.sum_product_decode(chan_out, final_out, iter, snr);
            finish = high_resolution_clock::now();
            serial_decode += duration_cast<duration<double>>(finish - start).count();
        }
        */
        MPI_Barrier(MPI_COMM_WORLD);
        start = high_resolution_clock::now();
        ldpc.get_llr_mpi(chan_out, llr_out, (int)out.size(), std_dev);
     //   printf("Got LLR values from QAM...\n");
    //    printf("LLR out size: %d\n", (int)llr_out.size());
    //    printf("In size: %d\n", (int)in.size());
    //    if (grank == 0) {
     //       print_vector(llr_out);
     //   }
     //   MPI_Barrier(MPI_COMM_WORLD);
        ldpc.sum_product_decode_mpi_block(llr_out, final_out, iter, snr);
    //    printf("LDPC decoding done...\n");
        finish = high_resolution_clock::now();
        mpi_decode += duration_cast<duration<double>>(finish - start).count();

        if (grank == 0) {
      //      print_vector(in);
       //     print_vector(out);
      //      printf("Out size: %d\n", (int)out.size());
      //      print_vector(llr_out);
      //      printf("LLR out size: %d\n", (int)llr_out.size());
       //     print_vector(final_out);
            float ber = 0;
            for (int i = 0; i < final_out.size(); i++) {
                ber += abs(in[i] - final_out[i]);
            }
            std::cout << "BER: " << ber/(float)final_out.size() << "\n";

        }
    }
    if (grank == 0) {
        printf("Serial construction time: %f secs\n", serial_construct);
        printf("Serial encoding time for %d vectors: %f secs\n", num_syms, serial_encode);
        printf("Serial decoding time for %d vectors: %f secs\n", num_syms, serial_decode);
        printf("MPI construction time with %d procs: %f secs\n", gsize, mpi_construct);
        printf("MPI encoding time with %d procs for %d vectors: %f secs\n", gsize, num_syms, mpi_encode);
        printf("MPI decoding time with %d procs for %d vectors: %f secs\n", gsize, num_syms, mpi_decode);
        std::ofstream ofs;
        std::string file = "mpi_time_" + std::to_string(gsize) + "_procs_" + std::to_string(n) + "_" + std::to_string(m) + "_" + std::to_string(k) + "_" + std::to_string(num_syms) + "_syms_" + std::to_string(num_times) + "_times.csv";
        ofs.open (file.c_str(), std::ofstream::out | std::ofstream::trunc);
        ofs << serial_construct << ",";
        ofs << serial_encode << ",";
        ofs << serial_decode << "\n";
        ofs << mpi_construct << ",";
        ofs << mpi_encode << ",";
        ofs << mpi_decode << "\n";
    }

    MPI_Finalize();
}