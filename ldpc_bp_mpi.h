#ifndef LDPC_BP_MPI_H
#define LDPC_BP_MPI_H

#include "ldpc_bp.h"
#include <mpi.h>
#include <thread>

class ldpc_bp_mpi : public ldpc_bp {

public:
    ldpc_bp_mpi() {}
    ldpc_bp_mpi(int,int);
    ~ldpc_bp_mpi();
    void bcast_H_mat(int);
    void bcast_G_mat(int);
    void H_mat_to_rref_form_mpi();
    void gen_mat_from_H_mat_mpi();
    int check_standard_form_mpi();
    int check_vector_mpi(std::vector<int> &);
    void encode_using_G_mat_mpi(std::vector<int> &, std::vector<int> &);
    void sum_product_decode_mpi(std::vector<float> &, std::vector<int> &, int, float);
    void sum_product_decode_mpi_min_sum(std::vector<float> &, std::vector<int> &, int, float);
    void sum_product_decode_mpi_block(std::vector<float> &, std::vector<int> &, int, float);
    void add_input_to_list_mpi(std::vector<float> &);
    std::vector<int> get_output_from_list_mpi();
    void belief_propagation_mpi(int, float);
    void belief_propagation_nonblock_mpi(int, float);
    void belief_propagation_mpi_min_sum(int, float);

protected:
    int grank, gsize;

};

#endif