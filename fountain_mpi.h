#ifndef FOUNTAIN_MPI_H
#define FOUNTAIN_MPI_H

#include <mpi.h>
#include <algorithm>
#include "ldpc_bp.h"
#include "qam_llr_mpi.h"

#define min(a,b)  ((a < b) ? a : b)

class fountain_mpi : public ldpc_bp,qam_llr_mpi {

public:
    fountain_mpi() {}
    fountain_mpi(int, int);
    fountain_mpi(int, int, int, int);
    ~fountain_mpi();
    void create_encoding_mat(int, int);
    void encode_using_G_mat_mpi_f(std::vector<int> &, std::vector<int> &);
    void encode_using_H_mat(std::vector<int> &, std::vector<int> &);
    void sum_product_decode_mpi_f(std::vector<float> &, std::vector<int> &, int, float);
    void sum_product_decode_mpi_block_f(std::vector<float> &, std::vector<int> &, int, float);
    void belief_propagation_mpi_f(int, float);

protected:
    int num_msg, num_coded, gsize, grank;

};

#endif