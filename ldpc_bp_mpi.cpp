#include "ldpc_bp_mpi.h"

ldpc_bp_mpi::ldpc_bp_mpi(int grank_, int gsize_) : grank(grank_), gsize(gsize_) {}

ldpc_bp_mpi::~ldpc_bp_mpi() {}

//Boradcasts H matrix from root process to others
void ldpc_bp_mpi::bcast_H_mat(int root) {
    //printf("Bcasting H mat to all procs...\n");
    int row_size, col_size;
    //For H mat
    if (grank == root) {
        row_size = H_mat.size();
        col_size = H_mat[0].size();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&row_size, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&col_size, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    //printf("Row size: %d, Col size: %d, Proc: %d\n", row_size, col_size, grank);
    if (grank != root) {
        H_mat.resize(row_size);
        for (int i = 0; i < H_mat.size(); i++) {
            H_mat[i].resize(col_size);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < H_mat.size(); i++) {
        MPI_Bcast(&H_mat[i][0], (int)H_mat[i].size(), MPI_INT, root, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Bcast(&n, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
   // printf("NMK vals %d %d %d proc %d\n", n, m, k, grank);
    /*
    //If compressed form of H matrix exists
    if (H_comp.size() > 0) {
        for (int i = 0; i < H_comp.size(); i++) {
            MPI_Bcast(&H_comp[i][0], (int)H_comp[i].size(), MPI_2INT, root, MPI_COMM_WORLD);
        }
    }
    */
   //printf("Have H mat at proc %d\n", grank);
}

//Broadcasts G matrix from root process to others
void ldpc_bp_mpi::bcast_G_mat(int root) {
    int row_size, col_size;
    //For G mat
    if (grank == root) {
        row_size = G_mat.size();
        col_size = G_mat[0].size();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&row_size, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&col_size, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    //printf("Row size: %d, Col size: %d, Proc: %d\n", row_size, col_size, grank);
    if (grank != root) {
        G_mat.resize(row_size);
        for (int i = 0; i < G_mat.size(); i++) {
            G_mat[i].resize(col_size);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < G_mat.size(); i++) {
        MPI_Bcast(&G_mat[i][0], (int)G_mat[i].size(), MPI_INT, root, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Bcast(&standard_form_var, 1, MPI_INT, root, MPI_COMM_WORLD);
    //printf("Std form var val is %d at proc %d\n", standard_form_var, grank);
    MPI_Barrier(MPI_COMM_WORLD);
    //printf("Have G mat at proc %d\n", grank);
}

//Checks if vector is error free by multiplying it with H matrix
int ldpc_bp_mpi::check_vector_mpi(std::vector<int> &vec) {
    int temp, flag = 0;
    if (H_comp.size() == 0) {
        H_mat_comp_form();
    }
    if (ceil((float)vec.size()/(float)n) < gsize) {
        for (int i = 0; i < ceil((float)vec.size()/(float)n); i++) {
            for (int j = 0; j < H_comp.size(); j += gsize) {
                if (j + grank >= H_comp.size()) {
                    break;
                }
                temp = 0;
                for (int jj = 0; jj < H_comp[j].size(); jj++) {
                    temp = (temp + (vec[H_comp[j][jj].col + i*n] * H_comp[j][jj].val)) % 2;
                }
                if (temp > 0) {
                    flag = 1;
                    break;
                }
            }
        }
    }
    //Using reduce and broadcast to check the flag value and send it to all processors so that they have a common return value
    if (grank == 0) {
        MPI_Reduce(MPI_IN_PLACE, &flag, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(&flag, &flag, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (flag != 0) {
        return -1;
    }
    return 0;
}

//Converts H matrix to rref form using mpi
void ldpc_bp_mpi::H_mat_to_rref_form_mpi() {
    std::vector<std::vector<int> > H_mat = this->H_mat;

    int numCols = H_mat[0].size(), numRows = H_mat.size(), c;

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
        for (int j = 0; j < numRows; j += gsize) {
            if (j + grank >= numRows) {
                break;
            }
            if (j + grank == i) {
                continue;
            }
            if (H_mat[j + grank][c] > 0) {
                for (int jj = 0; jj < numCols; jj++) {
                    H_mat[j + grank][jj] = (H_mat[j + grank][jj] + H_mat[i][jj]) % 2;
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        //printf("Copying values proc %d\n", grank);
        //Each procesor sends updated rows to all others
        //int loop_size;
        for (int i = 0; i < numRows; i++) {
            /*if ((int)numRows - i >= gsize) {
                loop_size = gsize;
            } else {
                loop_size = (int)numRows - i;
            }
            for (int proc = 0; proc < loop_size; proc++) {
                */
                MPI_Bcast((void *)&H_mat[i][0], (int)H_mat[i].size(), MPI_INT, i % gsize, MPI_COMM_WORLD);  
            //}
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    //printf("Finished conversion to rref form proc %d\n", grank);
    MPI_Barrier(MPI_COMM_WORLD);
    //print_matrix(H_mat);
    

    this->H_rref = H_mat;
    
}

void ldpc_bp_mpi::gen_mat_from_H_mat_mpi() {
    //printf("Converting to rref form proc %d\n", grank);
    //Converting the H matrix to reduced row echelon form
    H_mat_to_rref_form_mpi();
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

    //printf("Creating null space proc %d\n", grank);


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
        for (int j = c + 1; j < H_mat[0].size(); j += gsize) {
            if (j + grank >= H_mat[0].size()) {
                break;
            }
            if (H_temp[j + grank][i] > 0) {
                for (int jj = 0; jj < H_temp[i].size(); jj++) {
                    H_temp[j + grank][jj] = (H_temp[j + grank][jj] + H_temp[c][jj]) % 2;
                }
            }
        }
        for (int j = c + 1; j < H_mat[0].size(); j++) {
            MPI_Bcast((void *)&H_temp[j][0], (int)H_temp[j].size(), MPI_INT, (j - c - 1) % gsize, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);

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
    //printf("Tank from null space proc %d\n", grank);

    //Finding where the null space starts
    int mm = H_mat.size();
    for (int i = 0; i < H_mat.size(); i++) {
        if (H_temp[i][i] == 0) {
            mm = i;
            break;
        }
    }
    //printf("Copying into generator matrix proc %d\n", grank);

    //Copying the null of the H matrix into the generator matrix
    //printf("H mat size proc %d: %d\n", grank, n);
    int in = n - mm;
    G_mat.resize(in);
    //printf("G mat resized proc %d\n", grank);
   // std::cout << in << "\n";
    for (int i = 0; i < in; i++) {
        G_mat[i].resize(n);
        for (int j = 0; j < n; j++) {
            G_mat[i][j] = H_temp[i + (H_temp.size()) - in][j + (H_temp[0].size()) - n];
        }
    }
   // print_matrix(H_temp);
}

//Checks if generator matrix is in standard form
/*
int ldpc_bp_mpi::check_standard_form_mpi() {
    int elems_per_proc = ceil((float)G_mat.size()/(float)gsize);
    int flag = 0;
    for (int i = grank*elems_per_proc; i < std::min((int)G_mat.size(), (grank+1)*elems_per_proc); i++) {
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
    //Using reduce and braodcast to check the flag value and sned it to all processors so that they have a common return value
    MPI_Reduce(&flag, &flag, 1, MPI_INT, MPI_LOR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (flag != 0) {
        return -1;
    }
    return 0;
}
*/

//Encodes the input symbols using the generator matrix
void ldpc_bp_mpi::encode_using_G_mat_mpi(std::vector<int> &in, std::vector<int> &out) {
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
    

    //If the number of symbols to encode is less than the number of processors
    //Goes to block procesing only if the input length is very large
    if (ceil((float)len/(float)num_msg_bits) < 10*gsize) {
        
        //printf("proc > syms\n");
        //Size of vector that each processor takes
        int *size_of_proc_data, *displ;//, *size_of_copy_data;
        size_of_proc_data = (int *)malloc(gsize*sizeof(*size_of_proc_data));
        displ = (int *)malloc(gsize*sizeof(*displ));
        //Getting the elements per processor
        float elems_per_proc = ((float)(n - num_msg_bits)/(float)gsize);
        //printf("Elems per proc %d: %f\n", grank, elems_per_proc);
        int total_elems = 0;
        for (int i = 0; i < gsize - 1; i++) {
            if (i % 2 == 0) {
                size_of_proc_data[i] = (int)floor(elems_per_proc);
                displ[i] = total_elems;
                total_elems += (int)floor(elems_per_proc);
            } else {
                size_of_proc_data[i] = (int)ceil(elems_per_proc);
                displ[i] = total_elems;
                total_elems += (int)ceil(elems_per_proc);
            }
            
            //size_of_copy_data[i] = (std::min(num_msg_bits, (i+1)*elems_per_proc_copy) - i*elems_per_proc_copy);
        }
        if ((n - num_msg_bits) - total_elems >= 0) {
            size_of_proc_data[gsize - 1] = (n - num_msg_bits) - total_elems;
            displ[gsize - 1] = total_elems;
            total_elems += (n - num_msg_bits) - total_elems;
        } else {
            size_of_proc_data[gsize - 1] = 0;
            displ[gsize - 1] = total_elems;
            //total_elems += (n - num_msg_bits) - total_elems;
        }
        //if (grank == 0) {
        //    printf("Size per proc %d: %d\n", gsize - 1, size_of_proc_data[gsize - 1]);
        //}
        //printf("Encoding starting for %d proc\n", grank);
        for (int i = 0; i < ceil((float)len/(float)num_msg_bits); i++) {
            std::copy(in.begin() + i*num_msg_bits, in.begin() + (i+1)*num_msg_bits, out.begin() + i*n);

            for (int j = num_msg_bits + displ[grank]; j < std::min(n, num_msg_bits + displ[grank] + size_of_proc_data[grank]); j++) {
                out[j + i*n] = 0;
                for (int jj = 0; jj < num_msg_bits; jj++) {
                    out[j + i*n] = (out[j + i*n] + (in[jj + i*num_msg_bits] * G_mat[jj][j])) % 2; 
                }
            }
            //Gathering all outputs from each procesor
            MPI_Barrier(MPI_COMM_WORLD);
            if (grank == 0) {
                MPI_Gatherv(MPI_IN_PLACE, size_of_proc_data[grank], MPI_INT, (void *)&out[i*n + displ[grank] + num_msg_bits], size_of_proc_data, displ, MPI_INT, 0, MPI_COMM_WORLD);
            } else {
                MPI_Gatherv((void *)&out[i*n + displ[grank] + num_msg_bits], size_of_proc_data[grank], MPI_INT, (void *)&out[i*n + displ[grank] + num_msg_bits], size_of_proc_data, displ, MPI_INT, 0, MPI_COMM_WORLD);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast((void *)&out[0], (int)out.size(), MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        //free(size_of_copy_data);
        free(size_of_proc_data);
        free(displ);
        
        /*
        for (int i = 0; i < ceil((float)len/(float)num_msg_bits); i++) {
            //std::copy(in.begin() + i*num_msg_bits + grank*elems_per_proc_copy, in.begin() + i*num_msg_bits + std::min(num_msg_bits, (grank+1)*elems_per_proc_copy), out.begin() + i*n + grank*elems_per_proc_copy);

            for (int j = 0; j < n; j += gsize) {
                if (j + grank >= n) {
                    break;
                }
                out[j + grank + i*n] = 0;
                if ((j + grank) < num_msg_bits) {
                    out[j + i*n + grank] = in[i*num_msg_bits + j + grank];
                } else {
                    for (int jj = 0; jj < num_msg_bits; jj++) {
                        out[j + grank + i*n] = (out[j + grank + i*n] + (in[jj + i*num_msg_bits] * G_mat[jj][j + grank])) % 2; 
                    }
                }
            }
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        //std::vector out_temp(out_vec.size());
        for (int i = 0; i < (int)out.size(); i += gsize) {
            if (i + grank >= (int)out.size()) {
                //MPI_Gather((void *)&out[(i-gsize) + grank], 1, MPI_INT, (void *)&out[(i-gsize) + grank], 1, MPI_INT, 0, MPI_COMM_WORLD);
                break;
            }
            //Gathering all outputs from each processor
            if (grank == 0) {
                MPI_Gather(MPI_IN_PLACE, 1, MPI_INT, (void *)&out[i + grank], 1, MPI_INT, 0, MPI_COMM_WORLD);
            } else {
                MPI_Gather((void *)&out[i + grank], 1, MPI_INT, (void *)&out[i + grank], 1, MPI_INT, 0, MPI_COMM_WORLD);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        //Broadcasting final collected output to all processes
        MPI_Bcast((void *)&out[0], (int)out.size(), MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        */
        //If the number of symbols is much more than number of processors
    
    } else {
        int *size_of_proc_data, *displ;
        size_of_proc_data = (int *)malloc(gsize*sizeof(*size_of_proc_data));
        displ = (int *)malloc(gsize*sizeof(*displ));
        float elems_per_proc = ceil((float)len/(float)num_msg_bits)/(float)gsize;
        int total_elems = 0;
        for (int i = 0; i < gsize - 1; i++) {
            if (i % 2 == 0) {
                size_of_proc_data[i] = (int)floor(elems_per_proc);
                displ[i] = total_elems;
                total_elems += (int)floor(elems_per_proc);
            } else {
                size_of_proc_data[i] = (int)ceil(elems_per_proc);
                displ[i] = total_elems;
                total_elems += (int)ceil(elems_per_proc);
            }
            //size_of_copy_data[i] = (std::min(num_msg_bits, (i+1)*elems_per_proc_copy) - i*elems_per_proc_copy);
        }
        if (ceil((float)len/(float)num_msg_bits) - total_elems >= 0) {
            size_of_proc_data[gsize - 1] = (ceil((float)len/(float)num_msg_bits) - total_elems);
            displ[gsize - 1] = total_elems;
            //total_elems += in_vec_size/n - total_elems;
        } else {
            size_of_proc_data[gsize - 1] = 0;
            displ[gsize - 1] = total_elems;
            //total_elems += (n - num_msg_bits) - total_elems;
        }

        for (int i = displ[grank]; i < displ[grank] + size_of_proc_data[grank]; i++) {
            std::copy(in.begin() + i*num_msg_bits, in.begin() + (i+1)*num_msg_bits, out.begin() + i*n);

            for (int j = num_msg_bits; j < n; j++) {
                out[j + i*n] = 0;
                for (int jj = 0; jj < num_msg_bits; jj++) {
                    out[j + i*n] = (out[j + i*n] + (in[jj + i*num_msg_bits] * G_mat[jj][j])) % 2; 
                }
            }
        }
        for (int i = 0; i < gsize; i++) {
            size_of_proc_data[i] = size_of_proc_data[i]*n;
            displ[i] = displ[i]*n;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (grank == 0) {
            MPI_Gatherv(MPI_IN_PLACE, size_of_proc_data[grank], MPI_INT, (void *)&out[displ[grank]], size_of_proc_data, displ, MPI_INT, 0, MPI_COMM_WORLD);
        } else {
            MPI_Gatherv((void *)&out[displ[grank]], size_of_proc_data[grank], MPI_INT, (void *)&out[displ[grank]], size_of_proc_data, displ, MPI_INT, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast((void *)&out[0], (int)out.size(), MPI_INT, 0, MPI_COMM_WORLD);
        //MPI_Barrier(MPI_COMM_WORLD);
        free(size_of_proc_data);
        free(displ);
        
        /*
        for (int i = 0; i < (int)ceil((float)len/(float)num_msg_bits); i += gsize) {
            if (i + grank >= (int)ceil((float)len/(float)num_msg_bits)) {
                //printf("Breaking...%d\n", grank);
                break;
            }
            std::copy(in.begin() + (i+grank)*num_msg_bits, in.begin() + (i+1+grank)*num_msg_bits, out.begin() + (i+grank)*n);

            for (int j = num_msg_bits; j < n; j++) {
                out[j + (i+grank)*n] = 0;
                for (int jj = 0; jj < num_msg_bits; jj++) {
                    out[j + (i+grank)*n] = (out[j + (i+grank)*n] + (in[jj + (i+grank)*num_msg_bits] * G_mat[jj][j])) % 2; 
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        //printf("Num syms: %d\n", (int)ceil((float)len/(float)num_msg_bits));
        for (int i = 0; i < (int)ceil((float)len/(float)num_msg_bits); i += gsize) {
            //printf("Here...%d\n", grank);
            if (i + grank >= (int)ceil((float)len/(float)num_msg_bits)) {
                //printf("Breaking...%d\n", grank);
                //MPI_Gather((void *)&out[(i - gsize + grank)*n], n, MPI_INT, (void *)&out[(i - gsize + grank)*n], n, MPI_INT, 0, MPI_COMM_WORLD);
                break;
            }
            //Gathering all outputs from each processor
            if (grank == 0) {
                MPI_Gather(MPI_IN_PLACE, n, MPI_INT, (void *)&out[(i + grank)*n], n, MPI_INT, 0, MPI_COMM_WORLD);
                //printf("Copied...%d\n", grank);
            } else {
                MPI_Gather((void *)&out[(i + grank)*n], n, MPI_INT, (void *)&out[(i + grank)*n], n, MPI_INT, 0, MPI_COMM_WORLD);
                //printf("Copied...%d\n", grank);
            }
            
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast((void *)&out[0], (int)out.size(), MPI_INT, 0, MPI_COMM_WORLD);
        */
    } 
}

void ldpc_bp_mpi::sum_product_decode_mpi(std::vector<float> &in_vec, std::vector<int> &out_vec, int iter, float snr) {
    printf("N proc %d: %d\n", grank, n);
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
    //printf("Starting decoding...\n");
    for (int i = 0; i < in_vec_size/n; i++) {
        //Adding input vector to adjacency list
        std::copy(in_vec.begin() + i*n, in_vec.begin() + (i+1)*n, in_temp.begin());
        add_input_to_list(in_temp);
      //  print_vector(in_temp);
        belief_propagation_mpi(iter, snr);
       // print_vector(in_temp);
        std::vector<int> temp_vec = get_output_from_list();
        std::copy(temp_vec.begin(), temp_vec.begin() + G_mat.size(), out_vec.begin() + i*G_mat.size());
    }
}

//Blocked sum product decoding
void ldpc_bp_mpi::sum_product_decode_mpi_block(std::vector<float> &in_vec, std::vector<int> &out_vec, int iter, float snr) {
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
    int num_msg_bits = G_mat.size();
    std::vector<float> in_temp(n);
    //printf("Starting decoding...\n");
    if (in_vec_size/n < gsize) {
        for (int i = 0; i < in_vec_size/n; i++) {
            //Adding input vector to adjacency list
            std::copy(in_vec.begin() + i*n, in_vec.begin() + (i+1)*n, in_temp.begin());
            add_input_to_list(in_temp);
        //  print_vector(in_temp);
            belief_propagation_mpi(iter, snr);
        // print_vector(in_temp);
            std::vector<int> temp_vec = get_output_from_list();
            std::copy(temp_vec.begin(), temp_vec.begin() + G_mat.size(), out_vec.begin() + i*G_mat.size());
        }
    } else {
        int *size_of_proc_data, *displ;
        size_of_proc_data = (int *)malloc(gsize*sizeof(*size_of_proc_data));
        displ = (int *)malloc(gsize*sizeof(*displ));
        float elems_per_proc = (float)(in_vec_size/n)/(float)gsize;
        //printf("Elems per proc: %f\n", elems_per_proc);
        int total_elems = 0;
        for (int i = 0; i < gsize - 1; i++) {
            if (i % 2 == 0) {
                size_of_proc_data[i] = (int)floor(elems_per_proc);
                displ[i] = total_elems;
                total_elems += (int)floor(elems_per_proc);
            } else {
                size_of_proc_data[i] = (int)ceil(elems_per_proc);
                displ[i] = total_elems;
                total_elems += (int)ceil(elems_per_proc);
            }
        }
        if (in_vec_size/n - total_elems >= 0) {
            size_of_proc_data[gsize - 1] = (in_vec_size/n - total_elems);
            displ[gsize - 1] = total_elems;
            //total_elems += in_vec_size/n - total_elems;
        } else {
            size_of_proc_data[gsize - 1] = 0;
            displ[gsize - 1] = total_elems;
            //total_elems += (n - num_msg_bits) - total_elems;
        }
        for (int i = displ[grank]; i < displ[grank] + size_of_proc_data[grank]; i++) {
            //Adding input vector to adjacency list
            std::copy(in_vec.begin() + i*n, in_vec.begin() + (i+1)*n, in_temp.begin());
            add_input_to_list(in_temp);
        //  print_vector(in_temp);
            belief_propagation(iter, snr);
        // print_vector(in_temp);
            std::vector<int> temp_vec = get_output_from_list();
            std::copy(temp_vec.begin(), temp_vec.begin() + G_mat.size(), out_vec.begin() + i*G_mat.size());
        }
        for (int i = 0; i < gsize; i++) {
            size_of_proc_data[i] = size_of_proc_data[i]*G_mat.size();
            displ[i] = displ[i]*G_mat.size();
            //if (grank == 0)
           //     printf("Proc %d: %d\n", i, displ[i]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (grank == 0) {
            MPI_Gatherv(MPI_IN_PLACE, size_of_proc_data[grank], MPI_INT, (void *)&out_vec[displ[grank]], size_of_proc_data, displ, MPI_INT, 0, MPI_COMM_WORLD);
        } else {
            MPI_Gatherv((void *)&out_vec[displ[grank]], size_of_proc_data[grank], MPI_INT, (void *)&out_vec[displ[grank]], size_of_proc_data, displ, MPI_INT, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast((void *)&out_vec[0], (int)out_vec.size(), MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        free(size_of_proc_data);
        free(displ);
    }
}

//Sum product decoding
void ldpc_bp_mpi::belief_propagation_mpi(int iter, float snr) {
    //Initial LLR values
    for (int i = 0; i < var.size(); i++) {
        //The initial r value
        var[i].node_val = 2 * var[i].node_val * pow(10, snr/(float)10);
        //The initial LLR value value
        llr.llr[i] = var[i].node_val;
        //Updating the intrinsic LLR value of variable nodes
        for (int j = 0; j < var[i].conn_vertex.size(); j++) {
            llr.intrin_llr[var[i].conn_vertex[j]][i] = var[i].node_val;
        }
    }
    //printf("Calculated initial M, L and r values...\n");

    for (int it = 0; it < iter; it++) {
        //Checking the updated L value with the H matrix
        //std::vector<int> check_vec = get_output_from_list();
        //print_vector(check_vec);
      //  if (check_vector_mpi(check_vec) == 0) {
      //      return;
      //  }
        //Horizontal step
        //Each check node calculates the extrinsic LLR based on the LLRs of the variable nodes
        for (int i = 0; i < check.size(); i += gsize) {
            if (i + grank >= check.size()) {
                break;
            }
            for (int j = 0; j < check[i + grank].conn_vertex.size(); j++) {
                llr.extrin_llr[i + grank][check[i + grank].conn_vertex[j]] = 1;
                for (int k = 0; k < check[i + grank].conn_vertex.size(); k++) {
                    if (check[i + grank].conn_vertex[k] != check[i + grank].conn_vertex[j])
                        llr.extrin_llr[i + grank][check[i + grank].conn_vertex[j]] *= tanh(llr.intrin_llr[i + grank][check[i + grank].conn_vertex[k]]/2.0);
                }
                llr.extrin_llr[i + grank][check[i + grank].conn_vertex[j]] = log((1 + llr.extrin_llr[i + grank][check[i + grank].conn_vertex[j]])/(1 - llr.extrin_llr[i + grank][check[i + grank].conn_vertex[j]]));
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        //Extrinsic LLR updation
        int loop_size;
        for (int i = 0; i < check.size(); i += gsize) {
            if ((int)check.size() - i >= gsize) {
                loop_size = gsize;
            } else {
                loop_size = (int)check.size() - i;
            }
            for (int proc = 0; proc < loop_size; proc++) {
                //MPI_Bcast((void *)&llr.extrin_llr[i + proc][0], (int)llr.extrin_llr[i+proc].size(), MPI_FLOAT, proc, MPI_COMM_WORLD);
                
                for (int j = 0; j < check[i + proc].conn_vertex.size(); j++) {
                    if (proc != grank) {
                        if (grank == check[i + proc].conn_vertex[j]%gsize) {
                            MPI_Recv(&llr.extrin_llr[i  + proc][check[i + proc].conn_vertex[j]], 1, MPI_FLOAT, proc, proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        }
                    } else {
                        int receiver_rank = check[i + proc].conn_vertex[j]%gsize;
                        if (proc != receiver_rank) {
                            MPI_Send(&llr.extrin_llr[i + proc][check[i + proc].conn_vertex[j]], 1, MPI_FLOAT, receiver_rank, proc, MPI_COMM_WORLD);
                        }
                    }
                }
                
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        //Vertical step
        //Each variable node updates its own intrinsic LLR based on the LLR of the check nodes
       // if (it < iter - 1) {
            for (int i = 0; i < var.size(); i += gsize) {
                if (i + grank >= var.size()) {
                    break;
                }
                if (it < iter - 1) {
                    //Calculating value of LLR for each variable node
                    for (int j = 0; j < var[i + grank].conn_vertex.size(); j++) {
                        llr.intrin_llr[var[i + grank].conn_vertex[j]][i + grank] = var[i + grank].node_val;
                        for (int k = 0; k < var[i + grank].conn_vertex.size(); k++) {
                            if (var[i + grank].conn_vertex[j] != var[i + grank].conn_vertex[k])
                                llr.intrin_llr[var[i + grank].conn_vertex[j]][i + grank] += llr.extrin_llr[var[i + grank].conn_vertex[k]][i + grank];
                        }
                    }
                } else {
                    //Updating output LLR values
                    llr.llr[i + grank] = var[i + grank].node_val;
                    for (int j = 0; j < var[i + grank].conn_vertex.size(); j++) {
                        llr.llr[i + grank] += llr.extrin_llr[var[i + grank].conn_vertex[j]][i + grank];
                    }
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
       // }
        //Updating LLR values as all variable nodes
        for (int i = 0; i < var.size(); i += gsize) {
            if ((int)var.size() - i >= gsize) {
                loop_size = gsize;
            } else {
                loop_size = (int)var.size() - i;
            }
            for (int proc = 0; proc < loop_size; proc++) {
                //Final LLR values broadcased only for last iteration
                if (it == iter - 1) {
                    MPI_Bcast((void *)&llr.llr[i + proc], 1, MPI_INT, proc, MPI_COMM_WORLD);
                } else {
                    //Procs send the variable node LLR values
                    for (int j = 0; j < var[i + proc].conn_vertex.size(); j++) {
                        if (proc != grank) {
                            if (grank == var[i + proc].conn_vertex[j]%gsize) {
                                MPI_Recv(&llr.intrin_llr[var[i + proc].conn_vertex[j]][i + proc], 1, MPI_FLOAT, proc, proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            }
                        } else {
                            int receiver_rank = var[i + proc].conn_vertex[j]%gsize;
                            if (proc != receiver_rank) {
                                MPI_Send(&llr.intrin_llr[var[i + proc].conn_vertex[j]][i + proc], 1, MPI_FLOAT, receiver_rank, proc, MPI_COMM_WORLD);
                            }
                        }
                    }
                }
            }
        }
    }
}

/*
void ldpc_bp_mpi::belief_propagation_nonblock_mpi(int iter, float snr) {
    //Initial LLR values
    for (int i = 0; i < var.size(); i++) {
        //The initial r value
        var[i].node_val = 2 * var[i].node_val * pow(10, snr/(float)10);
        //The L value
        llr.llr[i] = var[i].node_val;
        //Updating the M value
        for (int j = 0; j < var[i].conn_vertex.size(); j++) {
            llr.intrin_llr[var[i].conn_vertex[j]][i] = var[i].node_val;
        }
    }
    //printf("Calculated initial M, L and r values...\n");

    //Initializing threads
    //std::thread t_send, t_recv;
    MPI_Request *request;
    for (int it = 0; it < iter; it++) {
        //Checking the updated L value with the H matrix
        //std::vector<int> check_vec = get_output_from_list();
        //print_vector(check_vec);
      //  if (check_vector_mpi(check_vec) == 0) {
      //      return;
      //  }
        //Horizontal step
        //Each check node calculates the extrinsic LLR based on the LLRs of the variable nodes
        for (int i = 0; i < check.size(); i += gsize) {
            if (i + grank >= check.size()) {
                break;
            }
            for (int j = 0; j < check[i + grank].conn_vertex.size(); j++) {
                llr.extrin_llr[i + grank][check[i + grank].conn_vertex[j]] = 1;
                for (int k = 0; k < check[i + grank].conn_vertex.size(); k++) {
                    if (check[i + grank].conn_vertex[k] != check[i + grank].conn_vertex[j])
                        llr.extrin_llr[i + grank][check[i + grank].conn_vertex[j]] *= tanh(llr.intrin_llr[i + grank][check[i + grank].conn_vertex[k]]/2.0);
                }
                llr.extrin_llr[i + grank][check[i + grank].conn_vertex[j]] = log((1 + llr.extrin_llr[i + grank][check[i + grank].conn_vertex[j]])/(1 - llr.extrin_llr[i + grank][check[i + grank].conn_vertex[j]]));
                //Sending output to other processes
                int receiver_rank = check[i + grank].conn_vertex[j]%gsize;
                MPI_Isend(&llr.extrin_llr[i + grank][check[i + grank].conn_vertex[j]], 1, MPI_FLOAT, receiver_rank, grank, MPI_COMM_WORLD, request);
                for (int proc = 0; proc < gsize; proc++) {
                    if (proc != grank) {
                        MPI_Recv(&llr.extrin_llr[i  + proc][check[i + grank].conn_vertex[j]], 1, MPI_FLOAT, proc, proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        //Extrinsic LLR updation
        int loop_size;
        for (int i = 0; i < check.size(); i += gsize) {
            if ((int)check.size() - i >= gsize) {
                loop_size = gsize;
            } else {
                loop_size = (int)check.size() - i;
            }
            for (int proc = 0; proc < loop_size; proc++) {
                //MPI_Bcast((void *)&llr.extrin_llr[i + proc][0], (int)llr.extrin_llr[i+proc].size(), MPI_FLOAT, proc, MPI_COMM_WORLD);
                
                for (int j = 0; j < check[i + proc].conn_vertex.size(); j++) {
                    if (proc != grank) {
                        if (grank == check[i + proc].conn_vertex[j]%gsize) {
                            MPI_Recv(&llr.extrin_llr[i  + proc][check[i + proc].conn_vertex[j]], 1, MPI_FLOAT, proc, proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        }
                    } else {
                        int receiver_rank = check[i + proc].conn_vertex[j]%gsize;
                        if (proc != receiver_rank) {
                            MPI_Send(&llr.extrin_llr[i + proc][check[i + proc].conn_vertex[j]], 1, MPI_FLOAT, receiver_rank, proc, MPI_COMM_WORLD);
                        }
                    }
                }
                
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        //Vertical step
        //Each variable node updates its own LLR based on the LLR of the check nodes
       // if (it < iter - 1) {
            for (int i = 0; i < var.size(); i += gsize) {
                if (i + grank >= var.size()) {
                    break;
                }
                if (it < iter - 1) {
                    //Calculating value of LLR for each variable node
                    for (int j = 0; j < var[i + grank].conn_vertex.size(); j++) {
                        llr.intrin_llr[var[i + grank].conn_vertex[j]][i + grank] = var[i + grank].node_val;
                        for (int k = 0; k < var[i + grank].conn_vertex.size(); k++) {
                            if (var[i + grank].conn_vertex[j] != var[i + grank].conn_vertex[k])
                                llr.intrin_llr[var[i + grank].conn_vertex[j]][i + grank] += llr.extrin_llr[var[i + grank].conn_vertex[k]][i + grank];
                        }
                    }
                } else {
                    //Updating output LLR values
                    llr.llr[i + grank] = var[i + grank].node_val;
                    for (int j = 0; j < var[i + grank].conn_vertex.size(); j++) {
                        llr.llr[i + grank] += llr.extrin_llr[var[i + grank].conn_vertex[j]][i + grank];
                    }
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
       // }
        //Updating LLR values as all variable nodes
        for (int i = 0; i < var.size(); i += gsize) {
            if ((int)var.size() - i >= gsize) {
                loop_size = gsize;
            } else {
                loop_size = (int)var.size() - i;
            }
            for (int proc = 0; proc < loop_size; proc++) {
                //Final LLR values broadcased only for last iteration
                if (it == iter - 1) {
                    MPI_Bcast((void *)&llr.llr[i + proc], 1, MPI_INT, proc, MPI_COMM_WORLD);
                } else {
                    //Procs send the variable node LLR values
                    for (int j = 0; j < var[i + proc].conn_vertex.size(); j++) {
                        if (proc != grank) {
                            if (grank == var[i + proc].conn_vertex[j]%gsize) {
                                MPI_Recv(&llr.intrin_llr[var[i + proc].conn_vertex[j]][i + proc], 1, MPI_FLOAT, proc, proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            }
                        } else {
                            int receiver_rank = var[i + proc].conn_vertex[j]%gsize;
                            if (proc != receiver_rank) {
                                MPI_Send(&llr.intrin_llr[var[i + proc].conn_vertex[j]][i + proc], 1, MPI_FLOAT, receiver_rank, proc, MPI_COMM_WORLD);
                            }
                        }
                    }
                }
            }
        }
    }
}
*/
