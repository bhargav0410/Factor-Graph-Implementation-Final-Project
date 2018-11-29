#include "ldpc_bp.h"

ldpc_bp::ldpc_bp() {}

ldpc_bp::ldpc_bp(int _n, int _m, int _k) : n(_n), m(_m), k(_k) {
    setMatSize(ceil((_n * _m)/_k), _n);
    std::cout << "H mat dims: " << getNumRows() << "x" << getNumCols() << "\n";
    create_adj_mat();
    gen_mat_from_H_mat();
}

ldpc_bp::~ldpc_bp() {}

//Function to create H matrix (calls the function to create adjacency matrix of the graph).
void ldpc_bp::create_H_mat(int _n, int _m, int _k) {
    setNMK(_n, _m, _k);
    setMatSize(ceil((_n * _m)/_k), _n);
    create_adj_mat();
}

//If a parity check matrix is created using some other function, then the parity check matrix can be set using this function.
//Then all other functions related to the parity check matrix can be used for the matrix created.
void ldpc_bp::set_H_mat(std::vector<std::vector<int> > &vec) {
    setAdjMat(vec);
}

//Used to create random permutations of the first sub-matrix of H matrix.
void ldpc_bp::shuffle_vec(std::vector<int> &vec) {
    int idx_1, idx_2, temp;
    srand(time(NULL));
    for (int i = 0; i < vec.size(); i++) {
        idx_1 = rand() % vec.size();
        idx_2 = rand() % vec.size();
        temp = vec[idx_1];
        vec[idx_1] = vec[idx_2];
        vec[idx_2] = temp;
    }
}

//Creates the parity check matrix in a unique manner
void ldpc_bp::create_H_mat_diff_form(int _n, int _m, int _k) {
    setNMK(_n, _m, _k);
    setMatSize(ceil((n * m)/k), n);
   // std::cout << "Starting "<< H_mat.size() << "x" << H_mat.at(0).size() << " matrix creation...\n";
    int numCols = getNumCols();
    int numRows = getNumRows();
    int stride  = numCols - numRows;
    //std::cout << "Here...\n";
    for (int i = stride; i < numCols; i++) {
        for (int j = 0; j < numRows; j++) {
            if ((i-stride) == j) {
                H_mat[j][i] = 1;
            } else {
                H_mat[j][i] = 0;
            }
        }
    }

  //  for (int i = 0; i < )

}

//Gets the parity check matrix from file
void ldpc_bp::set_H_mat_from_file(std::string H_file, int rows, int cols) {
    mat_from_file(H_file, rows, cols);
    rate = (float)(cols - rows)/(float)cols;
}

//Creates the parity check martix using Gallager's construction method in this case, but is used to create the adjacency matrix of a graph.
void ldpc_bp::create_adj_mat() {
    //setMatSize(ceil(n/k) * m, n);
    
    //Creating the first sub-matrix of the H matrix
   // std::cout << "Creating first sub-matrix...\n";
    for (int i = 0; i < ceil(n/k); i++) {
        for (int j = 0; j < getNumCols(); j++) {
            if (j >= i*k and j < (i+1)*k) {
                H_mat[i][j] = 1;
            } else {
                H_mat[i][j] = 0;
            }
        }
    }

    //std::cout << "Creating rest of sub-matrices...\n";
    std::vector<int> shuffled;
    for (int i = 0; i < getNumCols(); i++) {
        shuffled.push_back(i);
    }

    //Creating the rest of H matrix
    for (int sub = 1; sub < m; sub++) {
        //std::cout << "Creating sub-matrix "<< sub+1 <<"...\n";
        //Creating random permutations of the first sub-matrix        
        shuffle_vec(shuffled);
        
        int kk = 0;
        for (int i = sub * ceil(n/k); i < (sub+1) * ceil(n/k); i++) {
            //std::cout << i << std::endl;
            for (int j = 0; j < getNumCols(); j++) {
                H_mat[i][j] = H_mat[kk][shuffled[j]];
            }
            kk++;
        }
    }

    /*
    shuffle_vec(shuffled);
    int temp;
    for (int i = 0; i < getNumCols(); i++) {
        for (int j = 0; j < getNumRows(); j++) {
            temp = H_mat[j][i];
            H_mat[j][i] = H_mat[j][shuffled[i]];
            H_mat[j][shuffled[i]] = temp;
        }
    }
    */
    
    
}

//Converts the H matrix to reduced row echelon form
void ldpc_bp::H_mat_to_rref_form() {
    std::vector<std::vector<int> > H_mat = getAdjMat();

    int numCols = getNumCols(), numRows = getNumRows(), c;

    for (int i = 0; i < numRows; i++) {
       // std::cout << "Checking diagonal value...\n";
        //Checking if the diagonal value of the I part of the parity check matrix is 0, and swapping with a row which has value 1
        c = i;
        
        if (H_mat[i][i] == 0) {
            for (int ii = i+1; ii < numRows; ii++) {
                if (H_mat[ii][i] > 0) {
                  //  int temp;
                    //Swapping rows
                    std::swap(H_mat[i], H_mat[ii]);
                    /*
                    for (int j = 0; j < numCols; j++) {
                        temp = H_mat[i][j];
                        H_mat[i][j] = H_mat[c][j];
                        H_mat[c][j] = temp;
                    }
                    */
                    c = i;
                    break;
                }
            }
        }
        if (H_mat[i][i] == 0) {
            int flag = 0;
            for (int ii = i+1; ii < numCols; ii++) {
                if (H_mat[i][ii] > 0) {
                    c = ii;
                    break;
                } else {
                    for (int jj = i+1; jj < numRows; jj++) {
                        if (H_mat[jj][ii] > 0) {
                            std::swap(H_mat[jj], H_mat[i]);
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
      //  print_matrix(H_mat);
        
//    for (int i = 0; i < getNumRows(); i++) {
      //  std::cout << "Elimination of lower triangular values...\n";
        //Forward elimination step and back elimination step. Eliminates any non-zero elements in the lower trianular part of the I part of the parity check matrix.
        for (int j = 0; j < numRows; j++) {
            if (j == i) {
                continue;
            }
            if (H_mat[j][c] > 0) {
                for (int jj = 0; jj < numCols; jj++) {
                    H_mat[j][jj] = (H_mat[j][jj] + H_mat[i][jj]) % 2;
                }
            }
        }
    }
    //print_matrix(H_mat);


    this->H_rref = H_mat;
    
}

//Converting the parity check matrix to systematic form i.e. [A; I] form using Gauss-Jordan elimination
//Not completed as of yet (not required for now)
void ldpc_bp::H_mat_to_syst_form() {
    std::vector<std::vector<int> > H_mat = getAdjMat();
    int numRows = getNumRows(), numCols = getNumCols();

    //Gets the position of first 1 in each column
//    std::pair<int, int> first_one;
/*
    std::vector<std::pair<int, int> > list_of_first_ones;

    for(int i = 0; i < numCols; i++) {
        for (int j = 0; j < numRows; j++) {
            if (H_mat[j][i] > 0) {
                list_of_first_ones.push_back(std::make_pair(i, j));
                break;
            }
        }
        //Keeps the vector sorted in descending order
        for (int j = list_of_first_ones.size()-2; j >= 0 ; j--) {
            if (list_of_first_ones[list_of_first_ones.size()-1].second > list_of_first_ones[j].second) {
                std::swap(list_of_first_ones[list_of_first_ones.size()-1], list_of_first_ones[j]);
            } else {
                break;
            }
        }
    }
*/
    //Converts the H matrix to approximate lower triangular form
    int temp;
    int one_present = 0, kk = 0;
    for (int i = 0; i < ceil(n/k); i++) {
        for (int j = numCols-1; j >= 0; j--) {
            if (H_mat[i][j] > 0) {
                for (int ii = i-1; ii >= 0; ii++) {

                }
                for (int ii = 0; ii < numRows; ii++) {
                    temp = H_mat[ii][j];
                    H_mat[ii][j] = H_mat[ii][numCols-(ceil(n/k))+i];
                    H_mat[ii][numCols-(ceil(n/k))+i] = temp;
                }
            }
        }
    }

    for (int i = numCols-ceil(n/k); i < numCols; i++) {
        for (int j = ceil(n/k); j < numRows; j++) {
            if (H_mat[j][i] > 0) {
                for (int jj = 0; jj < numCols; jj++) {
                    //temp = H_mat[j][jj];
                    H_mat[j][jj] = (H_mat[j][jj] + H_mat[i-(numCols-ceil(n/k))][jj]) % 2;
                  //  H_mat[i-(numCols-(ceil(n/k)))][jj] = temp;
                }
            }
        }
    }

    this->H_syst = H_mat;

    print_matrix(H_syst);

}

//Converts H matrix to approximate lower triangle form using row and column permutations
//Very useful for linar time encoding
void ldpc_bp::H_mat_to_alt_form() {}

//Sorts the H matrix based on G matrix such that the encoding time is linear
void ldpc_bp::sort_H_mat_based_on_G_mat() {
    if (G_mat.size() == 0) {
        std::cout << "Sorting not possible. Create generator matrix...\n";
        return;
    }
    int count, num_msg_bits = G_mat.size()-1;
    std::vector<int> count_vec;
    if (H_comp.size() == 0) {
        H_mat_comp_form();
    }

    for (int i = 0; i < H_comp.size(); i++) {
        int count  = 0;
        for (int j = 0; j < H_comp[i].size(); j++) {
            if (H_comp[i][j].col > num_msg_bits) {
                count++;
            }
        }
        count_vec.push_back(count);
        int count_vec_size = count_vec.size(), temp;
        for (int j = i; j > 0; j--) {
            if (count_vec[j] < count_vec[j-1]) {
                std::swap(count_vec[j], count_vec[j-1]);
                std::swap(H_comp[j], H_comp[j-1]);
                std::swap(H_mat[j], H_mat[j-1]);
            } else {
                break;
            }
        }
        std::cout << count_vec[i] << std::endl;
    }
}

//Since the H matrix is sparse, this function searches the non-zero values and stores them with index, compressing the H matrix storage
void ldpc_bp::H_mat_comp_form() {
    int numRows = getNumRows(), numCols = getNumCols();
    H_comp.resize(numRows);
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            if (H_mat[i][j] > 0) {
                comp_form cf;
                cf.col = j;
                cf.val = H_mat[i][j];
                H_comp[i].push_back(cf);
            }
        }
    }
    /*
    for (int i = 0; i < H_comp.size(); i++) {
        std::cout << "|| ";
        for (int j = 0; j < H_comp[i].size(); j++) {
            std::cout << "(" << H_comp[i][j].col << "," << H_comp[i][j].val << ")";
        }
        std::cout << " ||\n";
    }
    std::cout << "\n";
    */
}


//Creates a generator matrix from H matrix after taking matrix invers of part of matrix
void ldpc_bp::gen_mat_from_H_mat_inv() {
    int numRows = getNumRows(), numCols = getNumCols(); 
    std::vector<std::vector<int> > A1(getNumRows(), std::vector<int> (getNumRows())), eye_mat(getNumRows(), std::vector<int> (getNumRows()));
    //Forming an I matrix where the matrix inverse will be stored
    
    //Shuffling the H matrix
    std::vector<std::vector<int> > H_mat = getAdjMat();
    std::vector<int> shuffled;
    for (int i = 0; i < getNumCols(); i++) {
        shuffled.push_back(i);
    }
    shuffle_vec(shuffled);
    int temp;
    for (int i = 0; i < getNumCols(); i++) {
        for (int j = 0; j < getNumRows(); j++) {
            temp = H_mat[j][i];
            H_mat[j][i] = H_mat[j][shuffled[i]];
            H_mat[j][shuffled[i]] = temp;
        }
    }
    

    for (int i = 0; i < A1.size(); i++) {
        for (int j = 0; j < A1[0].size(); j++) {
            A1[i][j] = H_mat[j][i];
            if (i == j) {
                eye_mat[i][j] = 1;
            }
        }
    }
    std::cout << "Before inv...\n";
    print_matrix(A1);

    //Taking matrix inverse
    for (int i = 0; i < numRows; i++) {
        //Checking diagonal value
        //If diagonal value is 0, then check subsequent rows
        if (A1[i][i] == 0) {
            for (int j = i+1; j < numRows; j++) {
                if (A1[j][i] > 0) {
                    std::swap(A1[j], A1[i]);
                    std::swap(eye_mat[j], eye_mat[i]);
                    break;
                }
            }
        }
        //If no 1 present in subsequent rows, then check column
        if (A1[i][i] == 0) {
            for (int j = i+1; j < numRows; j++) {
                if (A1[i][j] > 0) {
                    int temp;
                    for (int jj = 0; jj < numRows; jj++) {
                        temp = A1[jj][j];
                        A1[jj][j] = A1[jj][i];
                        A1[jj][i] = temp;
                        temp = eye_mat[jj][j];
                        eye_mat[jj][j] = eye_mat[jj][i];
                        eye_mat[jj][i] = temp;
                    }
                    break;
                }
            }
        }
        //If no 1 present in row or column then the matrix is singular
        if (A1[i][i] == 0) {
            std::cout << "Matrix is singular...\n";
            break;
        }
    //}

   // for (int i = 0; i < numRows; i++) {
        //Elimination of lower triangular values
        for (int j = i+1; j < numRows; j++) {
            if (A1[j][i] > 0) {
                for (int jj = 0; jj < numRows; jj++) {
                    A1[j][jj] = (A1[j][jj] + A1[i][jj]) % 2;
                    eye_mat[j][jj] = (eye_mat[j][jj] + eye_mat[i][jj]) % 2;
                }
            }
        }
    }
    print_matrix(A1);
    print_matrix(eye_mat);
}

//Creates a generator matrix from the parity check matrix
void ldpc_bp::gen_mat_from_H_mat() {
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
   // print_matrix(H_temp);
}

void ldpc_bp::print_matrix(std::vector<std::vector<int> > vec) {
    for (int i = 0; i < vec.size(); i++) {
        std::cout << "|| ";
        for (int j = 0; j < vec[i].size(); j++) {
            std::cout << vec[i][j];
        }
        std::cout << " ||\n";
    }
    std::cout << "\n";
}

//Converts the G and H matrices to standard form
void ldpc_bp::standard_form() {
    int temp, flag = 0;;
    if (G_mat.size() == 0) {
        std::cout << "Standard form not available. First create generator matrix...\n";
        return;
    }
    //Shifting rows to make sure that G is in [I:P] form, and if H matrix is constructed then shifting H matrix as well
    for (int i = 0; i < G_mat.size(); i++) {
        for (int j = 0; j < G_mat[0].size(); j++) {
            flag = 0;
            if (G_mat[i][j] > 0) {
                for (int ii = 0; ii < G_mat.size(); ii++) {
                    if (G_mat[ii][j] > 0 and ii != i) {
                        flag = 1;
                        break;
                    }
                }
                if (flag == 0) {
                   // std::cout << "Here...\n";
                    if (i != j) {
                        for (int jj = 0; jj < G_mat.size(); jj++) {
                            temp = G_mat[jj][j];
                            G_mat[jj][j] = G_mat[jj][i];
                            G_mat[jj][i] = temp;
                        }
                        if (H_mat.size() > 0) {
                            for (int jj = 0; jj < H_mat.size(); jj++) {
                                temp = H_mat[jj][j];
                                H_mat[jj][j] = H_mat[jj][i];
                                H_mat[jj][i] = temp;
                            }
                        }
                    }
                }
            }
        }
    }
}

void ldpc_bp::print_matrices() {
  //  std::vector<std::vector<int> > H_mat = getAdjMat();
    if (H_mat.size() > 0) {
        std::cout << "Parity check matrix:\n";
        print_matrix(H_mat);
    } else {
        std::cout << "Parity check matrix not created...\n";
    }

    if (G_mat.size() > 0) {
        std::cout << "Generator matrix:\n";
        print_matrix(G_mat);
    } else {
        std::cout << "Generator matrix not created...\n";
    }
}

//Checks whether the matrices generated are correct
int ldpc_bp::check_matrices() {
    int temp;
    for (int i = 0; i < G_mat.size(); i++) {
        for (int j = 0; j < H_mat.size(); j++) {
            temp = 0;
            for (int jj = 0; jj < n; jj++) {
                temp = (temp + (G_mat[i][jj] * H_mat[j][jj])) % 2;
            }
            if (temp > 0) {
                std::cout << "Matrices not correct...\n";
                return -1;
            }
        }
    }
    std::cout << "Matrices correct...\n";
    return 0;
}

//Creates an adjacency list from the compressed form of H matrix
void ldpc_bp::create_list_from_mat() {
    if (H_comp.size() == 0) {
        H_mat_comp_form();
    }
    var.resize(n);
    check.resize(H_comp.size());
   // Conn conn_var, conn_check;
    for (int i = 0; i < H_comp.size(); i++) {
        //std::cout << i << std::endl;
        for (int j = 0; j < H_comp[i].size(); j++) {
           // std::cout << j << " ";
            var[H_comp[i][j].col].list.push_back(&check[i]);
            var[H_comp[i][j].col].node_val = INF_VAL;
            check[i].list.push_back(&var[H_comp[i][j].col]);
            check[i].node_val = INF_VAL;
        }
    }
}

//Multiplies the input vector with the generator matrix and encodes it.
//Pass reference to input and output vectors (no need to specify input and output vector length)
void ldpc_bp::encode_using_G_mat(std::vector<int> &in, std::vector<int> &out) {
    if (G_mat.size() == 0) {
        if (H_mat.size() == 0) {
            std::cout << "Create matrices first...\n";
            return;
        } else {
            gen_mat_from_H_mat();
        }
    }
    long int len = in.size();
    if (fmod((float)len/(float)n, 1.0) != 0) {
        for (int i = 0; i < fmod((float)len/(float)n, 1.0); i++) {
            in.push_back(0);
        }
    }
    len = in.size();
    int num_msg_bits = G_mat.size();
    //Resizing the output vector
    if (out.size() != ceil(len*(n/G_mat.size()))) {
        out.resize(ceil(len*(n/G_mat.size())));
    }

    for (int i = 0; i < ceil(len/(float)n); i++) {
        memcpy(&out[i * n], &in[i * n], num_msg_bits*sizeof(int));
        for (int j = n - num_msg_bits; j < n; j++) {
            out[j + i*n] = 0;

            for (int jj = 0; jj < G_mat.size(); jj++) {
                out[j + i*n] = (out[j + i*n] + (in[jj] * G_mat[j + i*n][jj])) % 2; 
            }
        }
    }


}

//Set the n, m or k values
void ldpc_bp::setNMK(int _n = 0, int _m = 0, int _k = 0) {
    if (_n > 0) {
        n = _n;
    }
    if (_m > 0) {
        m = _m;
    }
    if (_k > 0) {
        k = _k;
    }
}

//Sets the rate and also resizes the generator matrix
void ldpc_bp::setRateAndPuncGenMat(int _rate) {
    if ((_rate * n) > G_mat.size()) {
        std::cout << "Rate not possible...\n";
    } else {
        rate = _rate;
        std::cout << "Rate set to " << rate << std::endl;
        std::cout << "Puncturing generator matrix...\n";
        G_mat.resize(ceil(rate * n));
    }
}

//Creates H and G matrix from given rate (incomplete)
void ldpc_bp::create_H_mat_based_on_rate(float _rate, int _n) {
    n = _n;

    float temp = n*(1-_rate);
    if ((int)floor(temp)%2 == 0) {
        temp = floor(temp);
    } else {
        temp = ceil(temp);
    }

    for (int i = 2; i < n; i++) {
        if ((int)temp % i == 0 and fmod((float)i/(float)(1 - _rate), 1.0) < 0.1) {
            m = i;
            break;
        }
    }
    k = ceil((float)m/(float)(1 - _rate));
    setMatSize(ceil((n * m)/k), n);
    create_adj_mat();
}

float ldpc_bp::getGenMatRate() {
    if (G_mat.size() > 0) {
        rate = (float)G_mat.size()/(float)G_mat[0].size();
        return rate; 
    } else {
        std::cout << "Rate cannot be determined from Generator matrix...\n";
        return -1;
    }
}

//Returns the rate based on the matrices created or given
float ldpc_bp::getRate() {
    if (rate > 0) {
        return rate;
    } else if (m > 0 and k > 0) {
        rate = 1-((float)m/(float)k);
        return rate;
    }
    else {
        std::cout << "Rate cannot be determined...\n";
        return -1;
    }
}