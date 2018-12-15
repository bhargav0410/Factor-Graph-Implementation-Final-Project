#ifndef LDPC_BP_MPI_H
#define LDPC_BP_MPI_H

#include "factor_graph.h"
#include "ldpc_bp.h"
#include <mpi.h>

int rank, size;

class ldpc_bp_mpi : private ldpc_bp, public factor_graph {

public:
    ldpc_bp_mpi();
    ~ldpc_bp_mpi();
    int check_standard_form_mpi();
    int check_vector_mpi(std::vector<int> &);
    void encode_using_G_mat_mpi(std::vector<int> &, std::vector<int> &);
    void sum_product_decode_mpi(std::vector<float> &, std::vector<int> &, int, float);
    void add_input_to_list_mpi(std::vector<float> &);
    std::vector<int> get_output_from_list_mpi();
    void belief_propagation_mpi(int, float);

private:

};

#endif