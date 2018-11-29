#ifndef LDPC_BP_H
#define LDPC_BP_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "factor_graph.h"

struct comp_form {
    int col, val;
};

class ldpc_bp : public Factor_graph {

public:

    ldpc_bp();
    ldpc_bp(int, int, int);
    ~ldpc_bp();
    void create_adj_mat();
    void create_list_from_mat();
    void setNMK(int, int, int);
    void create_H_mat(int, int, int);
    void create_H_mat_diff_form(int, int, int);
    void create_H_mat_based_on_rate(float, int);
    void set_H_mat_from_file(std::string, int, int);
    void store_H_mat_in_file(std::string, int, int);
    void sort_H_mat_based_on_G_mat();
    void H_mat_to_syst_form();
    void H_mat_to_rref_form();
    void H_mat_to_alt_form();
    void H_mat_comp_form();
    void gen_mat_from_H_mat();
    void gen_mat_from_H_mat_inv();
    void standard_form();
    int check_matrices();
    void print_matrix(std::vector<std::vector<int> >);
    void print_matrices();
    void set_H_mat(std::vector<std::vector<int> > &);
    void shuffle_vec(std::vector<int> &);
    void setRateAndPuncGenMat(int);
    float getGenMatRate();
    float getRate();

private:
    std::vector<Conn> var;  //Variable node list
    std::vector<Conn> check;    //Check node list
    std::vector<std::vector<int> > & H_mat = getAdjMat();   //Parity check matrix
    std::vector<std::vector<int> > H_syst, A, B, C, D, E, T, H_rref;
    std::vector<std::vector<comp_form> > H_comp;
    std::vector<std::vector<int> > G_mat;   //Generator matrix
    int n = 0, m = 0, k = 0;    //n specifies length of codeword, m specifies the number of check nodes per variable node, and k specifies number of variable nodes per check node
    float rate = 0;
};

#endif