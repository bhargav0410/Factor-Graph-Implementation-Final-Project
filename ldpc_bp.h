#ifndef LDPC_BP_H
#define LDPC_BP_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <fstream>
#include <thread>
#include "factor_graph.h"

#define INF_VAL 10000000.0

//Struct for compressed form of H matrix
struct comp_form {
    int col;
    int val;
};

struct llr_mats {
    std::vector<std::vector<float>> extrin_llr;
    std::vector<std::vector<float>> intrin_llr;
    std::vector<float> llr;
};

template <typename type>
void print_vector(std::vector<type> &vec) {
    std::cout << std::endl;
    std::cout << "|| "; 
    for (int i = 0; i < vec.size(); i++) {
        std::cout << vec[i] << " ";
    }
    std::cout << "||" << std::endl;
}

class ldpc_bp : public Factor_graph {

public:
    ldpc_bp();
    ldpc_bp(int, int, int);
    ~ldpc_bp();
    void create_adj_mat();
    void create_list_from_mat();
    void setNMK(int _n, int _m, int _k);
    void create_H_mat(int, int, int);
    void create_H_mat_diff_form(int, int, int);
    void create_H_mat_based_on_rate(float, int);
    void set_H_mat_from_file(std::string, int, int);
    void set_G_mat_from_file(std::string, int, int);
    void store_H_mat_in_file(std::string, int, int);
    void sort_H_mat_based_on_G_mat();
    void H_mat_to_syst_form();
    void H_mat_to_rref_form();
    void H_mat_to_alt_form();
    void H_mat_comp_form();
    void gen_mat_from_H_mat();
    void gen_mat_from_H_mat_inv();
    void store_H_mat_in_file(std::string);
    void store_G_mat_in_file(std::string);
    void standard_form();
    int check_standard_form();
    int check_matrices();
    template<typename type> void print_matrix(std::vector<std::vector<type> > &);
    void print_matrices();
    void set_H_mat(std::vector<std::vector<int> > &);
    void shuffle_vec(std::vector<int> &);
    void setRateAndPuncGenMat(int);
    float getGenMatRate();
    float getRate();
    int get_num_input_syms();
    void encode_using_G_mat(std::vector<int> &, std::vector<int> &);
    void sum_product_decode(std::vector<float> &, std::vector<int> &, int, float);
    void sum_product_encode(std::vector<int> &, std::vector<int> &, int);
    void add_input_to_list(std::vector<float> &);
    std::vector<int> get_output_from_list();
    void belief_propagation(int, float);
    void belief_propagation_omp(int, float);
    int check_vector(std::vector<int> &);

protected:
    std::vector<Conn> var;  //Variable node list
    std::vector<Conn> check;    //Check node list
    llr_mats llr;
    std::vector<std::vector<int> > & H_mat = getAdjMat();   //Parity check matrix
    std::vector<std::vector<int> > H_syst, H_rref;
    std::vector<std::vector<comp_form> > H_comp;
    std::vector<std::vector<int> > G_mat;   //Generator matrix
    int n = 0, m = 0, k = 0;    //n specifies length of codeword, m specifies the number of check nodes per variable node, and k specifies number of variable nodes per check node
    float rate = 0;
    int standard_form_var = 0;
};

#endif