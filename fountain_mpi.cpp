#include "fountain_mpi.h"

fountain_mpi::fountain_mpi(int _grank, int _gsize) {
    grank = _grank;
    gsize = _gsize;
}

fountain_mpi::fountain_mpi(int _grank, int _gsize, int m_, int n_) : num_msg(m_), num_coded(n_) {
    grank = _grank;
    gsize = _gsize;
    create_encoding_mat(m_, n_);
}

fountain_mpi::~fountain_mpi() {}

//Creates an encoding matrix which will be used to create a bipartite graph
void fountain_mpi::create_encoding_mat(int m_, int n_) {
    this->n = n_;
    num_msg = m_;
    num_coded = n_;
    H_mat.resize(num_coded, std::vector<int> (num_msg));
    int degree;
    //Random degree distribution based matrix creation
    for (int i = 0; i < H_mat.size(); i += gsize) {
        if ((i + grank) >= H_mat.size()) {
            break;
        }
        if (i + grank < num_msg) {
            //degree = 1;
            H_mat[i + grank][i + grank] = 1;
        }
        else {
            degree = (rand() % (num_msg - 1)) + 1;
            for (int j = 0; j < degree; j++) {
                H_mat[i + grank][rand() % num_msg] = 1;
            }
        }
        
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //Broadcasting the created matrix to all processes
    for (int i = 0; i < H_mat.size(); i++) {
        MPI_Bcast((void *)&H_mat[i][0], (int)H_mat[i].size(), MPI_INT, i % gsize, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //Creating generator matrix using transpose of H matrix
    G_mat.resize(H_mat[0].size(), std::vector<int> (H_mat.size()));
    for (int i = 0; i < G_mat.size(); i++) {
        for (int j = 0; j < G_mat[i].size(); j++) {
            G_mat[i][j] = H_mat[j][i];
        }
    }
    standard_form_var = 1;
    MPI_Barrier(MPI_COMM_WORLD);
}

//Encoding using fountain codes
void fountain_mpi::encode_using_H_mat(std::vector<int> &in, std::vector<int> &out) {
    if (H_comp.size() == 0) {
        H_mat_comp_form();
    }
    int len = in.size();
    int num_msg_bits = num_msg;

    //std::cout << "Num msg bits: " << num_msg_bits << "\n";
    //std::cout << "N: " << n << "\n";

    //Padding zeros to input vector (changes the length)
    while (fmod((float)len/(float)num_msg_bits, 1.0) != 0.0) {
      //  std::cout << "Len: " << len << std::endl;
        //for (int i = 0; i < fmod((float)len/(float)num_msg_bits, 1.0); i++) {
        in.push_back(0);
        len = in.size();
    }
    if (out.size() != ceil(num_coded*((float)len/(float)num_msg_bits))) {
        out.resize(ceil(num_coded*((float)len/(float)num_msg_bits)));
    }
    for (int num_vecs = 0; num_vecs < len/num_msg_bits; num_vecs++) {
        for (int i = 0; i < H_comp.size(); i += gsize) {
            if ((i + grank) >= H_comp.size()) {
                break;
            }
            out[i + grank + num_vecs*num_coded] = 0;
            for (int j = 0; j < H_comp[i].size(); j++) {
                out[i + grank + num_vecs*num_coded] = (out[i + grank + num_vecs*num_coded] + in[H_comp[i + grank][j].col + num_vecs*num_msg]) % 2;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        //Broadcasting the created matrix to all processes
        for (int i = 0; i < H_comp.size(); i++) {
            MPI_Bcast((void *)&out[i + num_coded*num_vecs], 1, MPI_INT, i % gsize, MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

//Encodes the input symbols using the generator matrix
void fountain_mpi::encode_using_G_mat_mpi_f(std::vector<int> &in, std::vector<int> &out) {
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
    } 
}

//Sum product decoding for fountain codes
void fountain_mpi::sum_product_decode_mpi_f(std::vector<float> &in_vec, std::vector<int> &out_vec, int iter, float snr) {
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
        belief_propagation_mpi_f(iter, snr);
       // print_vector(in_temp);
        std::vector<int> temp_vec = get_output_from_list();
        std::copy(temp_vec.begin(), temp_vec.begin() + G_mat.size(), out_vec.begin() + i*G_mat.size());
    }
}

//Blocked sum product decoding
void fountain_mpi::sum_product_decode_mpi_block_f(std::vector<float> &in_vec, std::vector<int> &out_vec, int iter, float snr) {
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
            belief_propagation_mpi_f(iter, snr);
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
           /* if (grank == 0)
                printf("Proc %d: %d\n", i, displ[i]);

            if (grank == 0)
                printf("Proc %d: %d\n", i, size_of_proc_data[i]);
            */
            //size_of_copy_data[i] = (std::min(num_msg_bits, (i+1)*elems_per_proc_copy) - i*elems_per_proc_copy);
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
        /*if (grank == 0)
            printf("Proc %d: %d\n", gsize - 1, displ[gsize - 1]);

        if (grank == 0)
            printf("Proc %d: %d\n", gsize - 1, size_of_proc_data[gsize - 1]);
*/
        for (int i = displ[grank]; i < displ[grank] + size_of_proc_data[grank]; i++) {
            //Adding input vector to adjacency list
            std::copy(in_vec.begin() + i*n, in_vec.begin() + (i+1)*n, in_temp.begin());
            add_input_to_list(in_temp);
            belief_propagation(iter, snr);
            std::vector<int> temp_vec = get_output_from_list();
            std::copy(temp_vec.begin(), temp_vec.begin() + G_mat.size(), out_vec.begin() + i*G_mat.size());
        }
        for (int i = 0; i < gsize; i++) {
            size_of_proc_data[i] = size_of_proc_data[i]*G_mat.size();
            displ[i] = displ[i]*G_mat.size();
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
void fountain_mpi::belief_propagation_mpi_f(int iter, float snr) {
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
        //Each variable node updates its own LLR based on the LLR of the check nodes
       // if (it < iter - 1) {
            for (int i = 0; i < var.size(); i += gsize) {
                if (i + grank >= var.size()) {
                    break;
                }
                if (it < iter - 1) {
                    //Calculating value of L and M for each variable node
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
        //int loop_size;
        for (int i = 0; i < var.size(); i += gsize) {
            if ((int)var.size() - i >= gsize) {
                loop_size = gsize;
            } else {
                loop_size = (int)var.size() - i;
            }
            for (int proc = 0; proc < loop_size; proc++) {
                if (it == iter - 1) {
                    MPI_Bcast((void *)&llr.llr[i + proc], 1, MPI_INT, proc, MPI_COMM_WORLD);
                } else {
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