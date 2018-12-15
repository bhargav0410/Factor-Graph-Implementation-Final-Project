#include "ldpc_bp_mpi.h"

ldpc_bp_mpi::ldpc_bp_mpi() {}

ldpc_bp_mpi::~ldpc_bp_mpi() {}

/*
void ldpc_bp_mpi::H_mat_to_rref_form_mpi() {
    //std::vector<std::vector<int> > H_mat = this->H_mat;

    //NUmber of elements per processor is different
    int *size_of_proc_data, *displ;
    size_of_proc_data = (int *)malloc(size*sizeof(*size_of_proc_data));
    displ = (int *)malloc(size*sizeof(*displ));

    int numCols = getNumCols(), numRows = getNumRows(), c;
    int elems_per_proc = (int)ceil((float)numRows/(float)size);

    for (int i = 0; i < size; i++) {
        size_of_proc_data[i] = (std::min(numRows, (rank+1)*elems_per_proc) - rank*elems_per_proc)*numCols;
        displ[i] = i*(numCols*elems_per_proc);
    }

    int *H_mat = 0;
    H_mat = (int *)malloc(numRows*numCols*sizeof(*H_mat));
    for (int i = 0; i < numRows; i++) {
        memcpy(&H_mat[i*numCols], &this->H_mat[i][0], numCols*sizeof(*H_mat));
    }


    for (int i = 0; i < numRows; i++) {
       // std::cout << "Checking diagonal value...\n";
        //Checking if the diagonal value of the parity check matrix is 0, and swapping with a row which has value 1
        c = i;
        if (rank == 0) {
            if (H_mat[i*numCols + i] == 0) {
                for (int ii = i+1; ii < numRows; ii++) {
                    if (H_mat[ii*numCols + i] > 0) {
                        //Swapping rows
                        std::swap(H_mat[i*numCols], H_mat[ii*numCols]);
                        c = i;
                        break;
                    }
                }
            }
            if (H_mat[i*numCols + i] == 0) {
                int flag = 0;
                for (int ii = i+1; ii < numCols; ii++) {
                    if (H_mat[i*numCols + ii] > 0) {
                        c = ii;
                        break;
                    } else {
                        for (int jj = i+1; jj < numRows; jj++) {
                            if (H_mat[jj*numCols + ii] > 0) {
                                std::swap(H_mat[jj*numCols], H_mat[i*numCols]);
                                c = ii;
                                flag = 1;
                                break;
                            }
                        }
                        if (flag == 1) {
                            break;
                        }
                    }
                }   
            }
        }
        
        MPI_Bcast(H_mat, numCols*numRows, MPI_INT, 0, MPI_COMM_WORLD);
        
//    for (int i = 0; i < getNumRows(); i++) {
      //  std::cout << "Elimination of lower triangular values...\n";
        //Forward elimination step and back elimination step. Eliminates any non-zero elements in the lower triangular part of the parity check matrix.
        for (int j = rank*elems_per_proc; j < std::min(numRows, (rank+1)*elems_per_proc); j++) {
            if (j == i) {
                continue;
            }
            if (H_mat[j*numCols + c] > 0) {
                for (int jj = 0; jj < numCols; jj++) {
                    H_mat[j*numCols + jj] = (H_mat[j*numCols + jj] + H_mat[i*numCols + jj]) % 2;
                }
            }
        }
        MPI_Gatherv(&H_mat[rank*elems_per_proc*numCols], size_of_proc_data[rank], MPI_INT, H_mat, size_of_proc_data, displ, MPI_INT, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(H_mat, numCols*numRows, MPI_INT, 0, MPI_COMM_WORLD);
    //print_matrix(H_mat);


    this->H_rref = H_mat;
    for (int i = 0; i < numRows; i++) {
        memcpy(&this->H_rref[i][0], &H_mat[i*numCols], numCols*sizeof(*H_mat));
    }
}

void gen_mat_from_H_mat_mpi() {
    //Converting the H matrix to reduced row echelon form
    H_mat_to_rref_form();
    std::vector<std::vector<int> > H_mat = this->H_rref;
    std::vector<std::vector<int> > H_temp;
    H_temp.resize(n);
    //Creating null space of H matrix from its rref form
    for (int i = 0; i < H_temp.size(); i++) {
        H_temp[i].resize(n + H_mat.size());
        for (int j = 0; j < H_temp[i].size(); j++) {
            if (j < H_mat.size()) {
                H_temp[i][j] = H_mat[j][i];
            } else {
                if ((j-H_mat.size()) == i) {
                    H_temp[i][j] = 1;
                }
            }
        }
    }


    int c;
    for (int i = 0; i < getNumRows(); i++) {
     //   std::cout << "Checking diagonal value...\n";
        //Checking if the diagonal value of the I part of the parity check matrix is 0, and swapping with a row which has value 1
        for (int j = 0; j < getNumCols(); j++) {
            if (H_temp[j][i] > 0) {
                c = j;
                break;
            }
        }
        

     //   std::cout << "Elimination of lower triangular values...\n";
        //Forward elimination step. Eliminates any non-zero elements in the lower triangular part to find the null of H matrix.
        for (int j = c + 1; j < getNumCols(); j++) {
            if (H_temp[j][i] > 0) {
                for (int jj = 0; jj < H_temp[i].size(); jj++) {
                    H_temp[j][jj] = (H_temp[j][jj] + H_temp[c][jj]) % 2;
                }
            }
        }

        if (H_temp[i][i] == 0) {
            while (c < getNumCols()) {
                if (H_temp[c][i] > 0) {
                    //Swapping rows
                    std::swap(H_temp[i], H_temp[c]);
                    break;
                }
                c++;
            }
        }
    }

    //Finding where the null space starts
    int mm = getNumRows();
    for (int i = 0; i < getNumRows(); i++) {
        if (H_temp[i][i] == 0) {
            mm = i;
            break;
        }
    }

    //Copying the null of the H matrix into the generator matrix
    int in = n - mm;
    G_mat.resize(in);
   // std::cout << in << "\n";
    for (int i = 0; i < in; i++) {
        G_mat[i].resize(n);
        for (int j = 0; j < n; j++) {
            G_mat[i][j] = H_temp[i + (H_temp.size()) - in][j + (H_temp[0].size()) - n];
        }
    }
}
*/

//Checks if vector is error free by cultiplying it with H matrix
int ldpc_bp_mpi::check_vector_mpi(std::vector<int> &vec) {

}

//Checks if generator matrix is in standard form
int ldpc_bp_mpi::check_standard_form_mpi() {
    int elems_per_proc = ceil((float)G_mat.size()/(float)size);
    int flag = 0;
    for (int i = rank*elems_per_proc; i < std::min((int)G_mat.size(), (rank+1)*elems_per_proc); i++) {
        for (int j = 0; j < G_mat.size(); j++) {
            if (i == j) {
                if (G_mat[i][j] == 0) {
                    flag = 1;
                }
            } else {
                if (G_mat[i][j] != 0) {
                    flag = 1;
                }
            }
        }
    }
    return 0;
}

//Encodes the input symbols using the generator matrix
void ldpc_bp_mpi::encode_using_G_mat_mpi(std::vector<int> &in_vec, std::vector<int> &out_vec) {
    if (G_mat.size() == 0) {
        if (H_mat.size() == 0) {
            std::cout << "Create matrices first...\n";
            return;
        } else {
            gen_mat_from_H_mat();
        }
    }
    int len = in.size();
    int num_msg_bits = G_mat.size();

    std::cout << "Num msg bits: " << num_msg_bits << "\n";
    std::cout << "N: " << n << "\n";

    //Padding zeros to input vector (changes the length)
    while (fmod((float)len/(float)num_msg_bits, 1.0) != 0.0) {
      //  std::cout << "Len: " << len << std::endl;
        //for (int i = 0; i < fmod((float)len/(float)num_msg_bits, 1.0); i++) {
        in.push_back(0);
        len = in.size();
    }
    //std::cout << "Len: " << len << "\n";
    //Resizing the output vector
    if (out.size() != ceil(n*((float)len/(float)num_msg_bits))) {
        out.resize(ceil(n*((float)len/(float)num_msg_bits)));
    }
    //std::cout << "Out size: " << out.size() << "\n";
    //Encoding the input vector
    for (int i = 0; i < ceil((float)len/(float)num_msg_bits); i++) {
        std::copy(in.begin() + i*num_msg_bits, in.begin() + (i+1)*num_msg_bits, out.begin() + i*n);
        for (int j = num_msg_bits; j < n; j++) {
            out[j + i*n] = 0;
            for (int jj = 0; jj < num_msg_bits; jj++) {
                out[j + i*n] = (out[j + i*n] + (in[jj + i*num_msg_bits] * G_mat[jj][j])) % 2; 
            }
        }
    }
}

void sum_product_decode_mpi(std::vector<float> &, std::vector<int> &, int, float) {}
void add_input_to_list_mpi(std::vector<float> &) {}
std::vector<int> get_output_from_list_mpi() {}
void belief_propagation_mpi(int, float) {}