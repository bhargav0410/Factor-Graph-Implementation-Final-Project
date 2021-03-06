#ifndef QAM_LLR_MPI_H
#define QAM_LLR_MPI_H

#include <vector>
#include <cstdlib>
#include <string>
#include <iostream>
#include <complex>
#include <cmath>
#include <algorithm>
#include <mpi.h>

#ifndef NUM_THREADS
    #define NUM_THREADS 32
#endif

struct constel_val {
    std::complex<float> const_place;
    std::vector<int> gray_str;
};

class qam_llr_mpi {

public:
    qam_llr_mpi();
    qam_llr_mpi(int, int);
    ~qam_llr_mpi();
    void qam_to_gray_mpi(std::vector<std::complex<float>> &, std::vector<int> &, int);
    void gray_to_qam_mpi(std::vector<int> &, std::vector<std::complex<float>> &);
    void get_llr_mpi(std::vector<std::complex<float>> &, std::vector<float> &, int, float);
    void set_contellation(int);
    std::vector<std::vector<std::complex<float>>> get_constellation_vals();
    void load_balancing_mpi(int*, int*, int, int);

protected:
    std::vector<std::vector<constel_val>> constellation;
    int bits_per_sym, qam_size;
    int qrank, qsize;

};

#endif