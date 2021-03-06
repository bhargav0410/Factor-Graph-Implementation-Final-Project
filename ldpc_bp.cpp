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
    std::ifstream outfile;
    outfile.open(H_file.c_str(), std::ofstream::in);
    H_mat.resize(rows);
    for (int i = 0; i < rows; i++) {
        H_mat[i].resize(cols);
        for (int j = 0; j < cols; j++) {
            outfile >> H_mat[i][j];
        }
    } 
    outfile.close();
    n = cols;
    rate = (float)(cols - rows)/(float)cols;
}

//Gets the generator matrix from file
void ldpc_bp::set_G_mat_from_file(std::string G_file, int rows, int cols) {
    std::ifstream outfile;
    outfile.open(G_file.c_str(), std::ofstream::in);
    G_mat.resize(rows);
    for (int i = 0; i < rows; i++) {
        G_mat[i].resize(cols);
        for (int j = 0; j < cols; j++) {
            outfile >> G_mat[i][j];
        }
    } 
    outfile.close();
    n = cols;
    rate = (float)(cols)/(float)rows;
}

//Creates the parity check martix using Gallager's construction method in this case, but is used to create the adjacency matrix of a graph.
void ldpc_bp::create_adj_mat() {
    //setMatSize(ceil(n/k) * m, n);
    
    //Creating the first sub-matrix of the H matrix
   // std::cout << "Creating first sub-matrix...\n";
    for (int i = 0; i < ceil(n/k); i++) {
        for (int j = 0; j < getNumCols(); j++) {
            if (j >= i*k && j < (i+1)*k) {
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

//Matrix to rref form
void mat_to_rref_form(std::vector<std::vector<int>> &H_mat) {
    int numRows = H_mat.size(), numCols = H_mat[0].size(), c;
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
}

//Copies part of SALT form in matrices in compressed form
std::vector<std::vector<comp_form>> H_salt_comp_form_copy(std::vector<std::vector<int>> &H_salt, int row_start, int row_end, int col_start, int col_end) {
    std::vector<std::vector<comp_form>> temp(row_end - row_start);
    for (int row = row_start; row < row_end; row++) {
        for (int col = col_start; col < col_end; col++) {    
            if (H_salt[row][col] > 0) {
               // printf("Copying into E matrix...\n");
                comp_form cf;
                cf.col = col - col_start;
                cf.val = H_salt[row][col];
                temp[row - row_start].push_back(cf);
            }
        }
    }
    return temp;
}

//Copies part of SALT form in matrices
void H_salt_copy(std::vector<std::vector<int>> &H_salt, std::vector<std::vector<int>> &temp, int row_start, int row_end, int col_start, int col_end) {
    temp.resize(row_end - row_start);
    for (int row = row_start; row < row_end; row++) {
        for (int col = col_start; col < col_end; col++) {    
            temp[row - row_start].push_back(H_salt[row][col]);
        }
    }
}

//Converts H matrix to approximate lower triangle form using row and column permutations
//Very useful for linear time encoding
void ldpc_bp::H_mat_to_salt_form() {
    //Initializing variables
    H_salt = this->H_mat;
    std::vector<std::vector<int>> D;
    int g_size, col_idx, row_idx, num_ones = 0, temp = 0;

    printf("Starting salt form conversion...\n");
    //Starting from last column to get to alt form
    for (int i = H_salt[0].size() - 1; i >= 0; i--) {
        //Checking number of ones per column
        if (i == H_salt[0].size() - 1) {
           // printf("For last column...\n");
            for (int col = i; col >= 0; col--) {
                temp = 0;
             //   printf("checking each row of each column...\n");
                for (int row = 0; row < H_salt.size(); row++) {
                    if (H_salt[row][col] > 0) {
                        temp += 1;
                    }
                }
                //Getting column index
                if (col == i || temp < num_ones) {
                    num_ones = temp;
                    col_idx = col;
                    if (num_ones == 1) {
                        break;
                    }
                }
            }
        } else {
            num_ones = H_salt.size();
            //For any column other than the last one, we only check for ones on and above the diagonal of matrix T
            for (int col = i; col >= 0; col--) {
                temp = 0;
                for (int row = 0; row < row_idx; row++) {
                    if (H_salt[row][col] > 0) {
                        temp += 1;
                    }
                }
                if (temp == 0) {
                    continue;
                }
                //Getting column index
                if (col == i || temp < num_ones) {
                    num_ones = temp;
                    col_idx = col;
                    if (num_ones == 1) {
                        break;
                    }
                }
            }
        }
     //   printf("Found column %d with num ones = %d...\n", col_idx, num_ones);

        //Swapping columns to get column with minimum number of ones at index i
        if (col_idx != i) {
            temp = 0;
            for (int row = 0; row < H_salt.size(); row++) {
                temp = H_salt[row][col_idx];
                H_salt[row][col_idx] = H_salt[row][i];
                H_salt[row][i] = temp;
            }
        }
    //    printf("Columns %d and %d swapped...\n", col_idx, i);

        //Getting the ones below diagonal of T matrix
        if (i == H_salt[0].size() - 1) {
           // printf("Starting row swap...\n");
            //For the last column only permute rows to bring ones to the last rows of the column
            int temp_idx = (int)H_salt.size() - 1;
            //Updating row index
            row_idx = (int)H_salt.size() - num_ones;
            for (int row = 0; row < row_idx; row++) {
                if (H_salt[row][i] > 0) {
                    for (int j = row_idx; j < H_salt.size(); j++) {
                        if (H_salt[j][i] == 0) {
                            H_salt[j].swap(H_salt[row]);
                        }
                    }
                }
            }
        } else {
            //For any column other than the last one
            while (num_ones > 0) {
                if (num_ones == 1) {
                    for (int row = 0; row < row_idx; row++) {
                        if (H_salt[row][i] > 0 && row < row_idx - 1) {
                            H_salt[row].swap(H_salt[row_idx - 1]);
                        }
                    }
                    num_ones -= 1;
                } else {
                    for (int row = 0; row < row_idx; row++) {
                        if (H_salt[row][i] > 0) {
                            for (int jj = (int)H_salt.size() - 1; jj > row; jj--) {
                                H_salt[row].swap(H_salt[jj]);
                            }
                            break;
                        }
                    }
                    num_ones -= 1;
                    row_idx -= 1;
                }
            }
            //Updating row index
            row_idx -= 1;
        }
    //    printf("Rows permuted...row index = %d\n", row_idx);
    //    std::cin.get();
        if (row_idx < 1) {
            break;
        }
    }

    //Gaussian elmination for making E part of salt form to 0

    //Finding row and column index for T matrix
    for (int row = 0; row < H_salt.size(); row++) {
        if (H_salt[row][(int)H_salt[0].size() - 1] > 0) {
            row_idx = row;
            break;
        }
    }
    for (int col = H_salt[0].size() - 1; col >= 0; col--) {
        if (H_salt[0][col] > 0) {
            col_idx = col;
            break;
        }
    }
    //Copying part fo salt form in A, B, C, D, E, T matrices
    T = H_salt_comp_form_copy(H_salt, 0, row_idx + 1, col_idx, (int)H_salt[0].size());
    E = H_salt_comp_form_copy(H_salt, row_idx + 1, (int)H_salt.size(), col_idx, (int)H_salt[0].size());
    std::vector<std::vector<int>> T_mat, E_mat, D_mat, temp_mat;
    H_salt_copy(H_salt, T_mat, 0, row_idx + 1, col_idx, (int)H_salt[0].size());
    H_salt_copy(H_salt, E_mat, row_idx + 1, (int)H_salt.size(), col_idx, (int)H_salt[0].size());

    //Gaussian elimination on temporary matrix to eliminate E part of SALT form
    temp_mat.resize(E_mat.size(), std::vector<int>(T_mat[0].size() + E_mat.size(), 0));
    std::vector<std::vector<int>> H_salt_temp = H_salt;
    for (int i = col_idx; i < H_salt[0].size(); i++) {
        for (int j = i - col_idx + 1; j < H_salt.size(); j++) {
            if (H_salt[j][i] > 0) {
                for (int jj = 0; jj < H_salt_temp[0].size(); jj++) {
                    H_salt_temp[j][jj] = (H_salt_temp[j][jj] + H_salt_temp[i - col_idx][jj]) % 2;
                }
            }
        }
    }
    for (int i = row_idx + 1; i < H_salt.size(); i++) {
        for (int j = 0; j < H_salt[0].size(); j++) {
            H_salt[i][j] = H_salt_temp[i][j];
        }
    }

    //Getting D part of salt form
    g_size = H_salt.size() - row_idx - 1;
    printf("g size: %d\n", g_size);
    H_salt_copy(H_salt, D_mat, row_idx + 1, (int)H_salt.size(), col_idx - g_size, col_idx);

    print_matrix(H_salt);
    print_matrix(T_mat);
    print_matrix(E_mat);
    print_matrix(D_mat);
    print_matrix(temp_mat);
    print_matrix(H_salt_temp);

    for (int i = 0; i < T.size(); i++) {
        std::cout << "|| ";
        for (int j = 0; j < T[i].size(); j++) {
            std::cout << "(" << T[i][j].col << "," << T[i][j].val << ")";
        }
        std::cout << " ||\n";
    }
    std::cout << "\n";
    for (int i = 0; i < E.size(); i++) {
        std::cout << "|| ";
        for (int j = 0; j < E[i].size(); j++) {
            std::cout << "(" << E[i][j].col << "," << E[i][j].val << ")";
        }
        std::cout << " ||\n";
    }
    std::cout << "\n";
    salt_form = 1;
}



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
    int numRows = H_mat.size(), numCols = H_mat[0].size();
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
    H_rref.clear();
    G_mat.clear();
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
    for (int i = 0; i < H_mat.size(); i++) {
     //   std::cout << "Checking diagonal value...\n";
        //Checking if the diagonal value of the I part of the parity check matrix is 0, and swapping with a row which has value 1
        for (int j = 0; j < H_mat[0].size(); j++) {
            if (H_temp[j][i] > 0) {
                c = j;
                break;
            }
        }
        

     //   std::cout << "Elimination of lower triangular values...\n";
        //Forward elimination step. Eliminates any non-zero elements in the lower triangular part to find the null of H matrix.
        for (int j = c + 1; j < H_mat[0].size(); j++) {
            if (H_temp[j][i] > 0) {
                for (int jj = 0; jj < H_temp[i].size(); jj++) {
                    H_temp[j][jj] = (H_temp[j][jj] + H_temp[c][jj]) % 2;
                }
            }
        }

        if (H_temp[i][i] == 0) {
            while (c < H_mat[0].size()) {
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
    int mm = H_mat.size();
    for (int i = 0; i < H_mat.size(); i++) {
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
   // print_matrix(H_rref);
    //print_matrix(H_temp);
}

//Print matrix
template<typename type>
void ldpc_bp::print_matrix(std::vector<std::vector<type> > &vec) {
    for (int i = 0; i < vec.size(); i++) {
        std::cout << "|| ";
        for (int j = 0; j < vec[i].size(); j++) {
            std::cout << vec[i][j];
        }
        std::cout << " ||\n";
    }
    std::cout << "\n";
}


//Functions to store created H and G matrices in file
void ldpc_bp::store_H_mat_in_file(std::string file) {
    std::ofstream outfile;
    outfile.open(file.c_str(), std::ofstream::out | std::ofstream::trunc);
    for (int i = 0; i < H_mat.size(); i++) {
        for (int j = 0; j < H_mat[i].size(); j++) {
            outfile << H_mat[i][j];
        }
    } 
    outfile.close();
}
void ldpc_bp::store_G_mat_in_file(std::string file) {
    std::ofstream outfile;
    outfile.open(file.c_str(), std::ofstream::out | std::ofstream::trunc);
    for (int i = 0; i < G_mat.size(); i++) {
        for (int j = 0; j < G_mat[i].size(); j++) {
            outfile << G_mat[i][j];
        }
    }
    outfile.close();
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
            //Checking each column of ith row
            if (G_mat[i][j] > 0) {
                for (int ii = 0; ii < G_mat.size(); ii++) {
                    if ((G_mat[ii][j] > 0) && (ii != i)) {
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
    standard_form_var = 1;
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


//Checks the encoded vector by multiplying it with the H matrix i.e. checks if H * vec = 0
int ldpc_bp::check_vector(std::vector<int> &vec) {
    int temp;
    if (H_comp.size() == 0) {
        H_mat_comp_form();
        //std::cout << "H mat comp form created...\n";
    }
    int flag = 0;
 //   std::cout << "n val: " << n << "\n";
    for (int i = 0; i < ceil((float)vec.size()/(float)n); i++) {
        #pragma omp parallel for num_threads(NUM_THREADS)
        for (int j = 0; j < H_comp.size(); j++) {
            temp = 0;
            for (int jj = 0; jj < H_comp[j].size(); jj++) {
                temp = (temp + (vec[H_comp[j][jj].col + i*n] * H_comp[j][jj].val)) % 2;
            }
            if (temp > 0 && flag == 0) {
                //std::cout << "Vector not correct...\n";
                flag = 1;
            }
        }
        #pragma omp barrier
        if (flag > 0) {
            break;
        }
    }
    if (flag > 0) {
        return -1;
    }
    //std::cout << "Vector correct...\n";
    return 0;
}

//Creates an adjacency list from the compressed form of H matrix.
//Two lists are created, one for variable nodes nad one for check nodes.
void ldpc_bp::create_list_from_mat() {
    if (H_comp.size() == 0) {
        H_mat_comp_form();
    }
    var.resize(H_mat[0].size());
    check.resize(H_comp.size());
    llr.extrin_llr.resize(H_mat.size());
    llr.intrin_llr.resize(H_mat.size());
    llr.llr.resize(H_mat[0].size());
   // Conn conn_var, conn_check;
    for (int i = 0; i < H_comp.size(); i++) {
        llr.extrin_llr[i].resize(H_mat[i].size());
        llr.intrin_llr[i].resize(H_mat[i].size());
        //std::cout << i << std::endl;
        for (int j = 0; j < H_comp[i].size(); j++) {
            var[H_comp[i][j].col].list.push_back(&check[i]);
            var[H_comp[i][j].col].node_val = 0;
            var[H_comp[i][j].col].vertex = H_comp[i][j].col;
            var[H_comp[i][j].col].conn_vertex.push_back(i);
            check[i].list.push_back(&var[H_comp[i][j].col]);
            check[i].node_val = 0;
            check[i].conn_vertex.push_back(H_comp[i][j].col);
            check[i].vertex = i;
        }
    }
}

//****************** Encoding part *****************

int ldpc_bp::check_standard_form() {
    for (int i = 0; i < G_mat.size(); i++) {
        for (int j = 0; j < G_mat[i].size(); j++) {
            if (i == j) {
                if (G_mat[i][j] == 0) {
                    return -1;
                }
            } else {
                if (G_mat[i][j] != 0) {
                    return -1;
                }
            }
        }
    }
    return 0;
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
    if (standard_form_var == 0) {
        standard_form();
    }
    int len = in.size();
    int num_msg_bits = G_mat.size();

    //std::cout << "Num msg bits: " << num_msg_bits << "\n";
    //std::cout << "N: " << n << "\n";

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
        #pragma omp parallel for num_threads(NUM_THREADS)
        for (int j = num_msg_bits; j < n; j++) {
            out[j + i*n] = 0;
            for (int jj = 0; jj < num_msg_bits; jj++) {
                out[j + i*n] = (out[j + i*n] + (in[jj + i*num_msg_bits] * G_mat[jj][j])) % 2; 
            }
        }
    }
}
//********************************************************


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
        if ((int)temp % i == 0 && fmod((float)i/(float)(1 - _rate), 1.0) < 0.1) {
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

//Gives the possible number of input symbols that can be encoded based on the final rate. The input vector size should be a multiple of this
int ldpc_bp::get_num_input_syms() {
    return (int)G_mat.size();
}

//Returns the rate based on the matrices created or given
float ldpc_bp::getRate() {
    if (rate > 0) {
        return rate;
    } else if (m > 0 && k > 0) {
        rate = 1-((float)m/(float)k);
        return rate;
    }
    else {
        std::cout << "Rate cannot be determined...\n";
        return -1;
    }
}


//******* Decoding part starts here *****************

//Decoding of signal added with AWGN noise
void ldpc_bp::sum_product_decode(std::vector<float> &in_vec, std::vector<int> &out_vec, int iter, float snr) {
    if (G_mat.size() == 0) {
        std::cout << "Need actual rate for decoding...Create generator matrix\n";
        return;
    }
    //Making sure input and output vectors have some values and are of correct sizes
    int in_vec_size = in_vec.size();
    if (in_vec_size == 0) {
        std::cout << "Input vector does not have input...\n";
        return;
    } else if (in_vec_size % n != 0) {
        std::cout << "Encoded vector size incorrect...\n";
        return;
    }
    if (out_vec.size() != G_mat.size()*(int)((float)in_vec_size/(float)n)) {
        out_vec.resize(G_mat.size()*(int)((float)in_vec_size/(float)n));
    }
    std::vector<float> in_temp(n);
    
    for (int i = 0; i < in_vec_size/n; i++) {
        //Adding input vector to adjacency list
        std::copy(in_vec.begin() + i*n, in_vec.begin() + (i+1)*n, in_temp.begin());
        add_input_to_list(in_temp);
      //  print_vector(in_temp);
        belief_propagation(iter, snr);
       // print_vector(in_temp);
        std::vector<int> temp_vec = get_output_from_list();
        std::copy(temp_vec.begin(), temp_vec.begin() + G_mat.size(), out_vec.begin() + i*G_mat.size());
    }
}

//Adds vector values to adjacency list
void ldpc_bp::add_input_to_list(std::vector<float> &vec) {
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int i = 0; i < vec.size(); i++) {
        var[i].node_val = vec[i];
    }
    /*
    if (vec.size() > G_mat.size()) {
        #pragma omp parallel for num_threads(NUM_THREADS)
        for (int i = G_mat.size(); i < var.size(); i++) {
            var[i].node_val = (float)vec[i];
        }
    }
    */
}

//Gets output vector from list
std::vector<int> ldpc_bp::get_output_from_list() {
    std::vector<int> out_vec(n);
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int j = 0; j < n; j++) {
        if (llr.llr[j] < 0) {
            out_vec[j] = 0;
        } else {
            out_vec[j] = 1;        
        }
    }
    return out_vec;
}



//Sum product (belief propagation) decoding and stores output in output vector provided
void ldpc_bp::belief_propagation(int iter, float snr) {
    //printf("Check node size: %d, Var node size: %d\n", (int)check.size(), (int)var.size()); 
    //Initial LLR values
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int i = 0; i < var.size(); i++) {
        //float prob = std::max((float)0.0, std::min((float)1.0, (float)(var[i].node_val + 1)/(float)2.0));
        //The initial r value
      //  var[i].node_val = 2 * var[i].node_val * pow(10, snr/(float)10);
        //The L value
        llr.llr[i] = var[i].node_val;
        //Updating the M value
        for (int j = 0; j < var[i].conn_vertex.size(); j++) {
            llr.intrin_llr[var[i].conn_vertex[j]][i] = var[i].node_val;
        }
    }
    //printf("Calculated initial M, L and r values...\n");

    for (int it = 0; it < iter; it++) {
        //Checking the updated L value with the H matrix
    //    std::vector<int> check_vec = get_output_from_list();
       // print_vector(check_vec);
    //    if (check_vector(check_vec) == 0) {
    //        return;
    //    }
        //Horizontal step
        //Each check node calculates the extrinsic LLR based on the LLRs of the variable nodes
        //The Extrinsic LLR value of each check node is updated.
        #pragma omp parallel for num_threads(NUM_THREADS)
        for (int i = 0; i < check.size(); i++) {
           // std::cout << "LLR size: " << check[i].llr.size() << "\n";
            for (int j = 0; j < check[i].conn_vertex.size(); j++) {
                llr.extrin_llr[i][check[i].conn_vertex[j]] = 1;
                for (int k = 0; k < check[i].conn_vertex.size(); k++) {
                    if (check[i].conn_vertex[k] != check[i].conn_vertex[j])
                        llr.extrin_llr[i][check[i].conn_vertex[j]] *= tanh(llr.intrin_llr[i][check[i].conn_vertex[k]]/2.0);
                }
                llr.extrin_llr[i][check[i].conn_vertex[j]] = log((1 + llr.extrin_llr[i][check[i].conn_vertex[j]])/(1 - llr.extrin_llr[i][check[i].conn_vertex[j]]));
            }
        }

        //Vertical step
        //Each variable node updates its own LLR based on the LLR of the check nodes
        #pragma omp parallel for num_threads(NUM_THREADS)
        for (int i = 0; i < var.size(); i++) {
            //Calculating value of L and M for each variable node
            if (it < iter - 1) {
                for (int j = 0; j < var[i].conn_vertex.size(); j++) {
                    llr.intrin_llr[var[i].conn_vertex[j]][i] = var[i].node_val;
                    for (int k = 0; k < var[i].conn_vertex.size(); k++) {
                        if (var[i].conn_vertex[j] != var[i].conn_vertex[k])
                            llr.intrin_llr[var[i].conn_vertex[j]][i] += llr.extrin_llr[var[i].conn_vertex[k]][i];
                    }
                }
            } else {
                llr.llr[i] = var[i].node_val;
                for (int j = 0; j < var[i].conn_vertex.size(); j++) {
                    //Updating output LLR values 
                    llr.llr[i] += llr.extrin_llr[var[i].conn_vertex[j]][i];
                }
            }
        }
    }
}

void ldpc_bp::belief_propagation_min_sum(int iter, float snr) {
    //printf("Check node size: %d, Var node size: %d\n", (int)check.size(), (int)var.size()); 
    //Initial LLR values
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int i = 0; i < var.size(); i++) {
        //float prob = std::max((float)0.0, std::min((float)1.0, (float)(var[i].node_val + 1)/(float)2.0));
        //The initial r value
      //  var[i].node_val = 2 * var[i].node_val * pow(10, snr/(float)10);
        //The L value
        llr.llr[i] = var[i].node_val;
        //Updating the M value
        for (int j = 0; j < var[i].conn_vertex.size(); j++) {
            llr.intrin_llr[var[i].conn_vertex[j]][i] = var[i].node_val;
        }
    }
    //printf("Calculated initial M, L and r values...\n");

    for (int it = 0; it < iter; it++) {
        //Checking the updated L value with the H matrix
    //    std::vector<int> check_vec = get_output_from_list();
       // print_vector(check_vec);
    //    if (check_vector(check_vec) == 0) {
    //        return;
    //    }
        //Horizontal step
        //Each check node calculates the extrinsic LLR based on the LLRs of the variable nodes
        //The Extrinsic LLR value of each check node is updated.
        #pragma omp parallel for num_threads(NUM_THREADS)
        for (int i = 0; i < check.size(); i++) {
           // std::cout << "LLR size: " << check[i].llr.size() << "\n";
            for (int j = 0; j < check[i].conn_vertex.size(); j++) {
                llr.extrin_llr[i][check[i].conn_vertex[j]] = 1;
                float temp_val = 1000000.0;
                for (int k = 0; k < check[i].conn_vertex.size(); k++) {
                    if (check[i].conn_vertex[k] != check[i].conn_vertex[j]) {
                        float val = llr.intrin_llr[i][check[i].conn_vertex[k]];
                        llr.extrin_llr[i][check[i].conn_vertex[j]] *= ((val >= 0) ? 1 : -1);
                        temp_val = (temp_val >= fabs(val)) ? fabs(val) : temp_val;
                    }
                }
                llr.extrin_llr[i][check[i].conn_vertex[j]] = llr.extrin_llr[i][check[i].conn_vertex[j]]*temp_val;
            }
        }

        //Vertical step
        //Each variable node updates its own LLR based on the LLR of the check nodes
        #pragma omp parallel for num_threads(NUM_THREADS)
        for (int i = 0; i < var.size(); i++) {
            //Calculating value of L and M for each variable node
            if (it < iter - 1) {
                for (int j = 0; j < var[i].conn_vertex.size(); j++) {
                    llr.intrin_llr[var[i].conn_vertex[j]][i] = var[i].node_val;
                    for (int k = 0; k < var[i].conn_vertex.size(); k++) {
                        if (var[i].conn_vertex[j] != var[i].conn_vertex[k])
                            llr.intrin_llr[var[i].conn_vertex[j]][i] += llr.extrin_llr[var[i].conn_vertex[k]][i];
                    }
                }
            } else {
                llr.llr[i] = var[i].node_val;
                for (int j = 0; j < var[i].conn_vertex.size(); j++) {
                    //Updating output LLR values 
                    llr.llr[i] += llr.extrin_llr[var[i].conn_vertex[j]][i];
                }
            }
        }
    }
}
