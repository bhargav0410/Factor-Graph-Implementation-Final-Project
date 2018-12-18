#ifndef FOUNTAIN_MPI_H
#define FOUNTAIN_MPI_H

//#include <mpi.h>
#include "ldpc_bp_mpi.h"

class fountain_mpi : public ldpc_bp_mpi {

public:
    fountain_mpi() {}
    fountain_mpi(int, int);
    fountain_mpi(int, int, int, int);
    ~fountain_mpi();
    void create_encoding_mat(int, int);
    void encode_using_H_mat(std::vector<int> &, std::vector<int> &);
    void sum_product_decode_mpi_f(std::vector<float> &, std::vector<int> &, int, float);
    void sum_product_decode_mpi_block_f(std::vector<float> &, std::vector<int> &, int, float);
    void belief_propagation_mpi_f(int, float);

protected:
    int num_msg, num_coded;

};

#endif