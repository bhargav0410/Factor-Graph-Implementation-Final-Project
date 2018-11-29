#include <iostream>
#include "factor_graph.h"
#include "ldpc_bp.h"
#include <cstdlib>
#include <chrono>
#include <algorithm>

using namespace std::chrono;

int main (int argc, char* argv[]) {
    //To measure execution time
    double serial_time = 1e30;
    duration<double> timediff;
    high_resolution_clock::time_point start, finish;

    int n = 12, m = 1, k = 4;
    float rate = 0;
    if (argc >= 4) {
        n = atoi(argv[1]);
        m = atoi(argv[2]);
        k = atoi(argv[3]);
        if (argc == 5) {
            rate = strtof(argv[4], NULL);
        }
    }

    std::cout << "Creating parity check matrix...\n";
    ldpc_bp ldpc;
    //ldpc.create_H_mat_diff_form(n, m, k);
    //ldpc.print_matrices();
    start = high_resolution_clock::now();
    if (argc == 5) {
        ldpc.create_H_mat_based_on_rate(rate, n);
    } else {
        ldpc.create_H_mat(n, m, k);
    } 
    finish = high_resolution_clock::now();
    timediff = duration_cast<duration<double>>(finish - start);
    serial_time = std::min(serial_time, timediff.count());
    std::cout << "Time: " << serial_time << std::endl;

    std::cout << "Anticipated Rate >= " << ldpc.getRate() << "\n";
    //ldpc.print_matrices();
    //std::cout << "Getting parity check matrix to systematic form matrix...\n";
    std::cout << "Creating generator matrix...\n";
    serial_time = 1e30;
    start = high_resolution_clock::now();
  //  ldpc.H_mat_to_rref_form();
    ldpc.gen_mat_from_H_mat();
    finish = high_resolution_clock::now();
    timediff = duration_cast<duration<double>>(finish - start);
    serial_time = std::min(serial_time, timediff.count());
    std::cout << "Time: " << serial_time << std::endl;

    std::cout << "Converting to standard form...\n";
    serial_time = 1e30;
    start = high_resolution_clock::now();
    ldpc.standard_form();
    finish = high_resolution_clock::now();
    timediff = duration_cast<duration<double>>(finish - start);
    serial_time = std::min(serial_time, timediff.count());
    std::cout << "Time: " << serial_time << std::endl;
    //std::cout << "Printing matrices...\n";

    ldpc.H_mat_comp_form();
    ldpc.create_list_from_mat();

    ldpc.check_matrices();
    ldpc.print_matrices();
    std::cout << "Final Rate = " << ldpc.getGenMatRate() << "\n";

    //ldpc.gen_mat_from_H_mat();

    //std::cin.get();
    return 0;

}