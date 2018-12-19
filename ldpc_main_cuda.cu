#include <iostream>
#include "factor_graph.h"
#include "ldpc_bp.h"
#include <cstdlib>
#include <chrono>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <random>
#include <iostream>
#include <string>
#include "ldpc_bp_cuda.cuh"

using namespace std::chrono;

int main (int argc, char* argv[]) {
    //To measure execution time
    double serial_encode = 0, cuda_encode = 0, serial_decode = 0, cuda_decode = 0;
    duration<double> timediff;
    high_resolution_clock::time_point start, finish;
    float snr = 10; //SNR in dB

    int n = 12, m = 1, k = 4, iter = 10, num_syms = 1, num_times = 1;
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
    }

    double serial_construct = 0, cuda_construct = 0;

    for (int times = 0; times < num_times; times++) {
        std::cout << "Creating parity check matrix...\n";
        ldpc_bp_cuda ldpc;
        ldpc.create_H_mat(n, m, k);

        std::cout << "Anticipated Rate >= " << ldpc.getRate() << "\n";
        std::cout << "Creating generator matrix...\n";
        start = high_resolution_clock::now();
        ldpc.gen_mat_from_H_mat();
        finish = high_resolution_clock::now();
        serial_construct += duration_cast<duration<double>>(finish - start).count();
        cuda_construct += ldpc.gen_mat_from_H_mat_cu();
        std::cout << "Converting to standard form...\n";
        ldpc.standard_form();
        ldpc.H_mat_comp_form();
        ldpc.create_list_from_mat();
        ldpc.check_matrices();

        //Creating input random vector and encoding it
        std::vector<int> in(num_syms*ldpc.get_num_input_syms()), out;
        srand(time(NULL));
        for (int i = 0; i < num_syms*ldpc.get_num_input_syms(); i++) {
            in[i] = (rand()%2);
        }
        
        //print_vector(in);
        //ldpc.print_matrices();
        start = high_resolution_clock::now();
        ldpc.encode_using_G_mat(in, out);
        finish = high_resolution_clock::now();
        serial_encode += duration_cast<duration<double>>(finish - start).count();
        cuda_encode += ldpc.encode_using_G_mat_cuda(in, out);
        std::cout << "Encoding done...\n";
        if (ldpc.check_vector(out) != 0) {
            std::cout << "Encoding incorrect...\n";
        }
        //print_vector(out);
        std::cout << "Final Rate = " << ldpc.getGenMatRate() << "\n";

        //Noise generation (equivalent to passing through a channel)

        //Noise variance based on input SNR
        float std_dev = pow((float)10.0, -((float)snr/(float)10.0));
        std::cout << "Noise power: " << std_dev << "\n";
        std::default_random_engine generator;
        std::normal_distribution<float> distribution(0.0, std_dev);
        std::vector<float> awgn(out.size()), chan_in(out.size()), chan_out(out.size());

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

        std::vector<int> final_out;

        //Decode noisy signal
        start = high_resolution_clock::now();
        ldpc.sum_product_decode(chan_out, final_out, snr, iter);
        finish = high_resolution_clock::now();
        final_out.clear();
        serial_decode += duration_cast<duration<double>>(finish - start).count();
        cuda_decode += ldpc.sum_product_decoding_cuda(chan_out, final_out, snr, iter);
        //print_vector(in);
        //print_vector(final_out);

        float ber = 0;
        for (int i = 0; i < final_out.size(); i++) {
            ber += abs(in[i] - final_out[i]);
        }
        std::cout << "BER: " << ber/(float)final_out.size() << "\n";
    }
    printf("Serial construction time: %f secs\n", serial_construct);
    printf("Serial encoding time: %f secs\n", serial_encode);
    printf("Serial decoding time: %f secs\n", serial_decode);
    printf("CUDA construction time: %f secs\n", cuda_construct);
    printf("CUDA encoding time: %f secs\n", cuda_encode);
    printf("CUDA decoding time: %f secs\n", cuda_decode);

    std::ofstream ofs;
    std::string file = "cuda_time_" + std::to_string(n) + "_" + std::to_string(m) + "_" + std::to_string(k) + "_" + std::to_string(num_syms) + "_syms_" + std::to_string(num_times) + "_times.csv";
    ofs.open (file.c_str(), std::ofstream::out | std::ofstream::trunc);
    ofs << serial_construct << ",";
    ofs << serial_encode << ",";
    ofs << serial_decode << "\n";
    ofs << cuda_construct << ",";
    ofs << cuda_encode << ",";
    ofs << cuda_decode << "\n";
    ofs.close();
    return 0;

}