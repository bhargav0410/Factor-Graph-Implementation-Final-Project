#include "qam_llr_mpi.h"
#include "ldpc_bp_mpi.h"
#include "ofdm_mpi.h"
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <ctime>
#include <cmath>
#include <random>
#include <cmath>
#include <algorithm>
#include <mpi.h>

using namespace std::chrono;

int main(int argc, char* argv[]) {

    //Initialzing MPI
    int provided;
    MPI_Init_thread(NULL, NULL, MPI_THREAD_SERIALIZED, &provided);
	//getting size and rank
	int gsize, grank;
	MPI_Comm_size(MPI_COMM_WORLD, &gsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &grank);

	double serial_construct = 0, mpi_construct = 0, serial_encode = 0, mpi_encode = 0, mpi_decode = 0, serial_decode = 0;
    duration<double> timediff;
    high_resolution_clock::time_point start, finish;

	int m = 1, n = 12, k = 4, snr = 10, iter = 10, qam_size = 4, num_syms = 10, num_times = 1, fft_size = 64, prefix_size = 16, num_ants = 1;
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
            qam_size = atoi(argv[7]);
        }
		if (argc >= 9) {
			fft_size = atoi(argv[8]);
		}
		if (argc >= 10) {
			prefix_size = atoi(argv[9]);
		}
		if (argc >= 11) {
			num_ants = atoi(argv[10]);
		}
		if (argc >= 12) {
			num_times = atoi(argv[11]);
		}
    }
	MPI_Barrier(MPI_COMM_WORLD);

	for (int times = 0; times < num_times; times++) {
        ldpc_bp_mpi ldpc(grank, gsize);
		ofdm_mpi ofdm(grank, gsize, fft_size, prefix_size, num_ants);
        ldpc.set_contellation(qam_size);
        std::vector<std::vector<std::complex<float>>> constel = ldpc.get_constellation_vals();
        std::cout << "Contellation dim 1 size: " << constel.size() << "\n";
        std::cout << "Contellation dim 2 size: " << constel[0].size() << "\n";
        for (int i = 0; i < constel.size(); i++) {
            for (int j = 0; j < constel[i].size(); j++) {
                std::cout << constel[i][j] << " ";
            }
            std::cout << "\n";
        }

        if (grank == 0) {
        //    std::cout << "Creating parity check matrix...\n";
            ldpc.create_H_mat(n, m, k);
        //    ldpc.H_mat_to_salt_form();
            std::cout << "Anticipated Rate >= " << ldpc.getRate() << "\n";
        //    std::cout << "Creating generator matrix...\n";
            start = high_resolution_clock::now();
            ldpc.gen_mat_from_H_mat();
            finish = high_resolution_clock::now();
            serial_construct += duration_cast<duration<double>>(finish - start).count();
            
        }
        MPI_Barrier(MPI_COMM_WORLD);
        ldpc.bcast_H_mat(0);
        ldpc.bcast_G_mat(0);
        std::cout << "Converting to standard form...\n";
        ldpc.standard_form();
    //    ldpc.check_matrices();
        MPI_Barrier(MPI_COMM_WORLD);

        ldpc.H_mat_comp_form();
        printf("Created H matrix in comp form ar proc %d\n", grank);
        ldpc.create_list_from_mat();
         printf("Created H matrix in comp form and list from matrix at proc %d\n", grank);
        //Creating input random vector and encoding it
        std::vector<int> in(num_syms*ldpc.get_num_input_syms()), out;
        std::vector<std::complex<float>> qam_out, ofdm_in;
		std::vector<std::vector<std::complex<float>>> ofdm_out;
        if (grank == 0) {
            srand(time(NULL));
            for (int i = 0; i < in.size(); i++) {
                in[i] = rand()%2;
            }
            print_vector(in);
            //ldpc.print_matrices();
            
        }
        MPI_Barrier(MPI_COMM_WORLD);
        printf("Created input vector...\n");
        MPI_Bcast((void *)&in[0], (int)in.size(), MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        printf("MPI encode...\n");
        start = high_resolution_clock::now();
        ldpc.encode_using_G_mat_mpi(in, out);
        printf("Processor %d encoding done.\n", grank);
        ldpc.gray_to_qam_mpi(out, qam_out);
        ofdm_in = qam_out;

		//OFDM based modulation of QAM symbols
		std::vector<std::complex<float>> pilot_in(fft_size - 1, std::complex<float>(1,0));
        print_vector(pilot_in);
        printf("Pilot vector size: %d\n", (int)pilot_in.size());
		std::vector<std::vector<std::complex<float>>> chan_est_in(num_ants, std::vector<std::complex<float>> (fft_size - 1, std::complex<float>(1,0)));
		print_vector(chan_est_in[0]);
        ofdm.set_chan_est(chan_est_in);
		ofdm.mod_one_frame_mpi(ofdm_in, ofdm_out, pilot_in);
        print_vector(ofdm_out[0]);
        finish = high_resolution_clock::now();
        mpi_encode += duration_cast<duration<double>>(finish - start).count();
        printf("Processor %d QAM encoding and OFDM modulation done.\n", grank);
        MPI_Barrier(MPI_COMM_WORLD);
        if (ldpc.check_vector_mpi(out) != 0) {
            printf("Processor %d encoding incorrect...\n", grank);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
        if (grank == 0) {
        //    std::cout << "Final Rate = " << ldpc.getGenMatRate() << "\n";
            //print_vector(out);
        }
        
		std::vector<std::complex<float>> awgn(ofdm_out[0].size()), chan_in(ofdm_out[0].size());
		std::vector<std::vector<std::complex<float>>> chan_out(ofdm_out.size(), std::vector<std::complex<float>> (ofdm_out[0].size()));
        float std_dev = (pow((float)10.0, -((float)snr/(float)10.0)));
        if (grank == 0) {
            std::cout << "Noise power: " << std_dev*std_dev << "\n";
            std::default_random_engine generator;
            std::normal_distribution<float> distribution(0.0, std_dev);

            //Passing through AWGN channel
			for (int n = 0; n < num_ants; n++) {
				for (int i = 0; i < awgn.size(); i++) {
					//Creating vector ready to be transmitted through channel
				   // chan_in[i] = 2*(float)out[i] - 1;
					//Creating noise for channel emulation
					awgn[i] = std::complex<float>(distribution(generator), distribution(generator));
					//Passing through AWGN channel
					chan_out[n][i] = ofdm_out[n][i] + awgn[i];
					//          awgn[i] = std::abs(awgn[i]);
				}
			}
      //      print_vector(awgn);
      //      printf("Chan out size: %d\n", (int)chan_out.size());
      //      print_vector(chan_out);
            printf("Signal passed through channel...\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
		for (int n = 0; n < num_ants; n++) {
			MPI_Bcast((void *)&chan_out[n][0], (int)chan_out[n].size(), MPI_COMPLEX, 0, MPI_COMM_WORLD);
		}
		std::vector<std::complex<float>> ofdm_demod_out;
        std::vector<float> llr_out;
        std::vector<int> final_out;

        print_vector(chan_out[0]);

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
		//Demodulating OFDM symbols after passing through AWGN channel
		ofdm.demod_one_frame_mpi(chan_out, ofdm_demod_out, pilot_in);
        //print_vector(ofdm_demod_out);
		//Resizing according to QAM input
		ofdm_demod_out.resize(qam_out.size());
        printf("QAM decoding input size: %d\n", (int)ofdm_demod_out.size());
        print_vector(ofdm_demod_out);
        ldpc.get_llr_mpi(ofdm_demod_out, llr_out, (int)out.size(), std_dev);
     //   printf("Got LLR values from QAM...\n");
    //    printf("LLR out size: %d\n", (int)llr_out.size());
    //    printf("In size: %d\n", (int)in.size());
    //    if (grank == 0) {
     //       print_vector(llr_out);
     //   }
     //   MPI_Barrier(MPI_COMM_WORLD);
        print_vector(llr_out);
        ldpc.sum_product_decode_mpi_block(llr_out, final_out, iter, snr);
        print_vector(final_out);
    //    printf("LDPC decoding done...\n");
        finish = high_resolution_clock::now();
        mpi_decode += duration_cast<duration<double>>(finish - start).count();
        MPI_Bcast((void *)&final_out[0], (int)final_out.size(), MPI_INT, 0, MPI_COMM_WORLD);
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
        MPI_Barrier(MPI_COMM_WORLD);
    }



	MPI_Finalize();
}