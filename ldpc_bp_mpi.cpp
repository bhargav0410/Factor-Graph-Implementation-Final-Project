#include "ldpc_bp_mpi.h"

ldpc_bp_mpi::ldpc_bp_mpi(int grank_, int gsize_) : grank(grank_),gsize(gsize_) {}

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
    //Using reduce and braodcast to check the flag value and send it to all processors so that they have a common return value
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
    if (ceil((float)len/(float)num_msg_bits) < gsize) {
        //printf("proc > syms\n");
        
        //printf("proc > syms\n");
        //Size of vector that each processor takes
        int *size_of_proc_data, *displ;//, *size_of_copy_data;
        size_of_proc_data = (int *)malloc(gsize*sizeof(*size_of_proc_data));
        //size_of_copy_data = (int *)malloc(gsize*sizeof(*size_of_proc_data));
        displ = (int *)malloc(gsize*sizeof(*displ));
        //Getting the elements per processor
        int elems_per_proc = (int)ceil((float)(n - num_msg_bits)/(float)gsize);
        //printf("Elems per proc: %d\n", elems_per_proc);
        //int elems_per_proc_copy = (int)ceil((float)num_msg_bits/(float)gsize);
        //Since each processor takes different size vectors
        for (int i = 0; i < gsize; i++) {
            size_of_proc_data[i] = (std::min(n - num_msg_bits, (i+1)*elems_per_proc) - i*elems_per_proc);
            //printf("Size per proc: %d\n", size_of_proc_data[i]);
            //size_of_copy_data[i] = (std::min(num_msg_bits, (i+1)*elems_per_proc_copy) - i*elems_per_proc_copy);
            displ[i] = i*elems_per_proc;
        }
        //printf("Encoding starting for %d proc\n", grank);
        for (int i = 0; i < ceil((float)len/(float)num_msg_bits); i++) {
            std::copy(in.begin() + i*num_msg_bits, in.begin() + (i+1)*num_msg_bits, out.begin() + i*n);

            for (int j = num_msg_bits + grank*elems_per_proc; j < std::min(n, num_msg_bits + (grank+1)*elems_per_proc); j++) {
                out[j + i*n] = 0;
                for (int jj = 0; jj < num_msg_bits; jj++) {
                    out[j + i*n] = (out[j + i*n] + (in[jj + i*num_msg_bits] * G_mat[jj][j])) % 2; 
                }
            }
            //Gathering all outputs from each procesor
            MPI_Barrier(MPI_COMM_WORLD);
            if (grank == 0) {
                MPI_Gatherv(MPI_IN_PLACE, size_of_proc_data[grank], MPI_INT, (void *)&out[i*n + grank*elems_per_proc + num_msg_bits], size_of_proc_data, displ, MPI_INT, 0, MPI_COMM_WORLD);
            } else {
                MPI_Gatherv((void *)&out[i*n + grank*elems_per_proc + num_msg_bits], size_of_proc_data[grank], MPI_INT, (void *)&out[i*n + grank*elems_per_proc + num_msg_bits], size_of_proc_data, displ, MPI_INT, 0, MPI_COMM_WORLD);
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
        
        //printf("proc < syms\n");
        //Size of the output vector that each procesor takes
        int *size_of_proc_data, *displ;
        size_of_proc_data = (int *)malloc(gsize*sizeof(*size_of_proc_data));
        displ = (int *)malloc(gsize*sizeof(*displ));
        //Getting the elements per processor
        int elems_per_proc = (int)ceil((float)ceil((float)len/(float)num_msg_bits)/(float)gsize);
        //Since each processor takes different size vectors
        for (int i = 0; i < gsize; i++) {
            size_of_proc_data[i] = (std::min((int)ceil((float)len/(float)num_msg_bits), (i+1)*elems_per_proc) - i*elems_per_proc)*n;
            displ[i] = i*elems_per_proc*n;
        }

        for (int i = grank*elems_per_proc; i < std::min((int)ceil((float)len/(float)num_msg_bits), (grank+1)*elems_per_proc); i++) {
            std::copy(in.begin() + i*num_msg_bits, in.begin() + (i+1)*num_msg_bits, out.begin() + i*n);

            for (int j = num_msg_bits; j < n; j++) {
                out[j + i*n] = 0;
                for (int jj = 0; jj < num_msg_bits; jj++) {
                    out[j + i*n] = (out[j + i*n] + (in[jj + i*num_msg_bits] * G_mat[jj][j])) % 2; 
                }
            }
        }
        if (grank == 0) {
            MPI_Gatherv(MPI_IN_PLACE, size_of_proc_data[grank], MPI_INT, (void *)&out[grank*elems_per_proc*n], size_of_proc_data, displ, MPI_INT, 0, MPI_COMM_WORLD);
        } else {
            MPI_Gatherv((void *)&out[grank*elems_per_proc*n], size_of_proc_data[grank], MPI_INT, (void *)&out[grank*elems_per_proc*n], size_of_proc_data, displ, MPI_INT, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast((void *)&out[0], (int)out.size(), MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
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

void ldpc_bp_mpi::add_input_to_list_mpi(std::vector<float> &) {}
std::vector<int> ldpc_bp_mpi::get_output_from_list_mpi() {}

//Sum product decoding
void ldpc_bp_mpi::belief_propagation_mpi(int iter, float snr) {
    //printf("Check node size: %d, Var node size: %d\n", (int)check.size(), (int)var.size()); 
    //Initial LLR values
    for (int i = 0; i < var.size(); i++) {
        //float prob = std::max((float)0.0, std::min((float)1.0, (float)(var[i].node_val + 1)/(float)2.0));
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

    for (int it = 0; it < iter; it++) {
        //Checking the updated L value with the H matrix
        //std::vector<int> check_vec = get_output_from_list();
        //print_vector(check_vec);
      //  if (check_vector_mpi(check_vec) == 0) {
      //      return;
      //  }
        //Horizontal step
        //Each check node calculates the extrinsic LLR based on the LLRs of the variable nodes
        //The Extrinsic LLR value of each check node is updated.
        for (int i = 0; i < check.size(); i += gsize) {
            if (i + grank >= check.size()) {
                break;
            }
           // std::cout << "LLR size: " << check[i].llr.size() << "\n";
            for (int j = 0; j < check[i + grank].conn_vertex.size(); j++) {
                llr.extrin_llr[i + grank][check[i + grank].conn_vertex[j]] = 0;
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

        //Vertical step
        //Each variable node updates its own LLR based on the LLR of the check nodes
       // if (it < iter - 1) {
            for (int i = 0; i < var.size(); i += gsize) {
                if (i + grank >= var.size()) {
                    break;
                }
                //Calculating value of L and M for each variable node
                llr.llr[i + grank] = var[i + grank].node_val;
                for (int j = 0; j < var[i + grank].conn_vertex.size(); j++) {
                    llr.intrin_llr[var[i + grank].conn_vertex[j]][i + grank] = var[i + grank].node_val;
                    for (int k = 0; k < var[i + grank].conn_vertex.size(); k++) {
                        if (var[i + grank].conn_vertex[j] != var[i + grank].conn_vertex[k])
                            llr.intrin_llr[var[i + grank].conn_vertex[j]][i + grank] += llr.extrin_llr[var[i + grank].conn_vertex[k]][i + grank];
                    }
                    //Updating output LLR values 
                    llr.llr[i + grank] += llr.extrin_llr[var[i + grank].conn_vertex[j]][i + grank];
                }
                for (int proc = 0; proc < gsize; proc++) {
                    
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
       // }
        //Updating LLR values as all variable nodes
        //int loop_size;
        for (int i = 0; i < var.size(); i += gsize) {
            if ((int)var.size() - i >= gsize) {
                loop_size = gsize;
            } else {
                loop_size = (int)var.size() - i;
            }
            for (int proc = 0; proc < loop_size; proc++) {
                MPI_Bcast((void *)&llr.llr[i + proc], 1, MPI_INT, proc, MPI_COMM_WORLD);
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