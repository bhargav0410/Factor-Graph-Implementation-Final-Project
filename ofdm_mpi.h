#ifndef OFDM_MPI_H
#define OFDM_MPI_H

#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include <fftw3.h>
#include <fftw3-mpi.h>

class ofdm_mpi : public qam_llr_mpi {

public:
    ofdm_mpi(int, int);
    ofdm_mpi(int, int, int, int);
    ~ofdm_mpi();
    void load_balancing_mpi(int *, int *, int, int);
    void create_pilot_vector_mpi(std::vector<std::vector<std::complex<float>>> &);
    void fft_mpi(std::vector<std::complex<float>> *);
    void ifft_mpi(std::vector<std::complex<float>> *);
    void chan_est_update_mpi(std::vector<std::vector<std::complex<float>>> &);
    void maximal_ratio_combining_mpi(std::vector<std::vector<std::complex<float>>> &, std::vector<std::complex<float>> &);
    void maximal_ratio_transmission_mpi(std::vector<std::complex<float>> &, std::vector<std::vector<std::complex<float>>> &);
    void mod_one_frame_mpi(std::vector<std::vector<std::complex<float>>> &, std::vector<std::vector<std::complex<float>>> &);
    void demod_one_frame_mpi(std::vector<std::vector<std::complex<float>>> &, std::vector<std::vector<std::complex<float>>> &);
    int get_fft_size();
    void set_fft_size(int);
    int get_prefix_size();
    void set_prefix_size(int);

protected:
    int orank, osize, fft_size, prefix_size;
    std::vector<std::vector<std::complex<float>>> chan_est, pilot_vector;

};

#endif