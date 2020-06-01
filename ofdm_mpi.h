#ifndef OFDM_MPI_H
#define OFDM_MPI_H

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <mpi.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <omp.h>

#ifndef NUM_THREADS
    #define NUM_THREADS 32
#endif

class ofdm_mpi {

public:
    ofdm_mpi(int, int);
    ofdm_mpi(int, int, int, int);
    ~ofdm_mpi();
    void load_balancing_mpi(int *, int *, int, int);
    void create_pilot_vector_mpi(int);
    void fft_mpi(std::vector<std::complex<float>> *);
    void ifft_mpi(std::vector<std::complex<float>> *);
    void fft_omp_mpi(std::complex<float> *, std::complex<float> *, int);
    void ifft_omp_mpi(std::complex<float> *, std::complex<float> *, int);
    void chan_est_update_mpi(std::vector<std::vector<std::complex<float>>> &);
    void maximal_ratio_combining_mpi(std::vector<std::vector<std::complex<float>>> &, std::vector<std::complex<float>> &);
    //void maximal_ratio_transmission_mpi(std::vector<std::complex<float>> &, std::vector<std::vector<std::complex<float>>> &);
    void mod_one_frame_mpi(std::vector<std::vector<std::complex<float>>> &, std::vector<std::vector<std::complex<float>>> &);
    void demod_one_frame_mpi(std::vector<std::vector<std::complex<float>>> &, std::vector<std::vector<std::complex<float>>> &);
    void swap_halves(std::complex<float> *, int);
    void divide_by_vec_omp(std::complex<float> *, std::complex<float> *, std::complex<float> *, int, int);
    int get_fft_size();
    void set_fft_size(int);
    int get_prefix_size();
    void set_prefix_size(int);
    void set_chan_est(std::vector<std::vector<std::complex<float>>>);
    std::vector<std::vector<std::complex<float>>> get_chan_est();
    void set_pilot(std::vector<std::vector<std::complex<float>>>);
    std::vector<std::vector<std::complex<float>>> get_pilot();

protected:
    int orank, osize, fft_size, prefix_size, num_ants;
    std::vector<std::vector<std::complex<float>>> chan_est, pilot;

};

#endif