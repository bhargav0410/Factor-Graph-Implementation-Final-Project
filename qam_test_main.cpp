#include "qam_llr_mpi.h"
#include "ldpc_bp.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <cmath>
#include <algorithm>
#include <mpi.h>

using namespace std::chrono;

int main (int argc, char* argv[]) {

    //Initialzing MPI
    MPI_Init(NULL, NULL);
    //getting size and rank
    int gsize, grank;
    MPI_Comm_size(MPI_COMM_WORLD, &gsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &grank);

    float snr = 10; //SNR in dB

    int qam_size = 16, num_syms = 10;
    qam_llr_mpi qam(gsize, grank);
    qam.set_contellation(qam_size);

    std::vector<int> in(qam_size*num_syms), out;    
    if (grank == 0) {
        srand(time(NULL));
        for (int i = 0; i < in.size(); i++) {
            in[i] = rand()%2;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast((void *)&in[0], (int)in.size(), MPI_INT, 0, MPI_COMM_WORLD);

    //QAM encoding
    qam.gray_to_qam_mpi(in, out);
    if (grank == 0) {
        print_vector(in);
        print_vector(out);
    }

    //Passing through channel
    std::vector<std::complex<float>> awgn(out.size()), chan_in(out.size()), chan_out(out.size());
    float std_dev = pow((float)10.0, -((float)snr/(float)10.0));
    std::cout << "Noise power: " << std_dev << "\n";
    if (grank == 0) {
        std::default_random_engine generator;
        std::normal_distribution<float> distribution(0.0, std_dev);

        //Passing through AWGN channel
        for (int i = 0; i < awgn.size(); i++) {
            //Creating noise for channel emulation
            awgn[i] = (distribution(generator), distribution(generator));
            //Passing through AWGN channel
            chan_out[i] = out[i] + awgn[i];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast((void *)&chan_out[0], (int)chan_out.size(), MPI_COMPLEX, 0, MPI_COMM_WORLD);
    std::vector<float> final_out;

    //QAM to LLR conversion
    qam.get_llr_mpi(chan_out, final_out, (int)in.size(), std_dev);
    if (grank == 0) {
        print_vector(chan_out);
        print_vector(final_out);
    }

    MPI_Finalize();
}

