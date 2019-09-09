#include "qam_llr_mpi.h"

qam_llr_mpi::qam_llr_mpi() {}

qam_llr_mpi::qam_llr_mpi(int _rank, int _size) {
    qrank = _rank;
    qsize = _size;
}

qam_llr_mpi::~qam_llr_mpi() {}

//Distributed the tasks amongst workers in a balanced fashion
void qam_llr_mpi::load_balancing_mpi(int *size_of_proc_data, int *displ, int num_procs, int len) {
    int displ_of_proc = 0;
    for (int i = 0; i < num_procs; i++) {
        displ[i] = displ_of_proc;
        size_of_proc_data[i] = (int)floor((float)len/(float)num_procs);
        if (i < (int)len % num_procs) {
            size_of_proc_data[i] += 1;
        }
        displ_of_proc += size_of_proc_data[i];
    }
}

//Gives gray code output based on minimum distance match between the input noisy qam symbol and the constellation points
void qam_llr_mpi::qam_to_gray_mpi(std::vector<std::complex<float>> &in, std::vector<int> &out, int) {

}

//Divides the binary string and gives equivalent qam output based on the qam size and input binary string (assumes that it is a string of gray coded binary values)
void qam_llr_mpi::gray_to_qam_mpi(std::vector<int> &in, std::vector<std::complex<float>> &out) {
    int prev_len = in.size();
    //Checking the input vector length and appending zeros to make it a multiple of bits per symbol of the qam constellation
    while (in.size() % bits_per_sym != 0) {
        in.push_back(0);
    }
    out.resize((int)in.size()/bits_per_sym);

    //Load balancing
    int *size_of_proc_data, *displ;
    size_of_proc_data = (int *)malloc(qsize*sizeof(*size_of_proc_data));
    displ = (int *)malloc(qsize*sizeof(*displ));
    int out_size = out.size();
    load_balancing_mpi(size_of_proc_data, displ, qsize, out_size);
   // printf("Data size for proc %d: %d, Displacement for proc %d: %d\n", qrank, size_of_proc_data[qrank], qrank, displ[qrank]);

    //qam output for gray code input
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int i = displ[qrank]; i < displ[qrank] + size_of_proc_data[qrank]; i++) {
     //   printf("Proc %d\n", qrank);
        std::vector<int> temp(bits_per_sym);
        std::copy(in.begin() + i*bits_per_sym, in.begin() + (i+1)*bits_per_sym, temp.begin());
        int idx = 0;
        for (int k = 0; k < bits_per_sym; k++) {
            idx += temp[k] * pow(2, bits_per_sym - 1 - k);
        }
       // printf("Index: %d\n", idx);
        if (bits_per_sym <= 1) {
            out[i] = 2*idx - 1;
        } else {
            out[i] = std::complex<float>(std::real(constellation[(int)floor((float)idx/(float)constellation.size())][idx % (int)constellation.size()].const_place), std::imag(constellation[(int)floor((float)idx/(float)constellation.size())][idx % (int)constellation.size()].const_place));
        }
        /*
        int flag = 0;
        for (int j = 0; j < qam_size; j++) {
            for (int k = 0; k < qam_size; k++) {
                if (constellation[j][k].gray_str == temp) {
                    //Setting output value here
                    out[i] = constellation[j][k].const_place;
                    flag = 1;
                    break;
                }
            }
            if (flag == 1) {
                break;
            }
        }
        */
    }


    //Broadcasting data to all other processes 
    MPI_Barrier(MPI_COMM_WORLD);
    if (qrank == 0) {
        MPI_Gatherv(MPI_IN_PLACE, size_of_proc_data[qrank], MPI_COMPLEX, (void *)&out[displ[qrank]], size_of_proc_data, displ, MPI_COMPLEX, 0, MPI_COMM_WORLD);
    } else {
        MPI_Gatherv((void *)&out[displ[qrank]], size_of_proc_data[qrank], MPI_COMPLEX, (void *)&out[displ[qrank]], size_of_proc_data, displ, MPI_COMPLEX, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast((void *)&out[0], (int)out.size(), MPI_COMPLEX, 0, MPI_COMM_WORLD);
    //Freeeing allocated memory and resizing input vector
    in.resize(prev_len);
    free(size_of_proc_data);
    free(displ);
}

//Gives LLR values for qam input
void qam_llr_mpi::get_llr_mpi(std::vector<std::complex<float>> &in, std::vector<float> &out, int out_len, float noise) {
    out.resize(out_len);
 //   printf("In size: %d\n", (int)in.size());
    //Load balancing
    int *size_of_proc_data, *displ;
    size_of_proc_data = (int *)malloc(qsize*sizeof(*size_of_proc_data));
    displ = (int *)malloc(qsize*sizeof(*displ));
    //Distributing jobs amongst all workers
    load_balancing_mpi(size_of_proc_data, displ, qsize, (int)in.size());
    //LLR output for qam input
    float llr_for_zero, llr_for_one;
   // int num_threads_for_llr = std::min(NUM_THREADS, size_of_proc_data[qrank]/2);

    //Using OpenMP for local parallelism
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int i = displ[qrank]; i < displ[qrank] + size_of_proc_data[qrank]; i++) {
        for (int bit = 0; bit < bits_per_sym; bit++) {
            float llr_for_zero = 0;
            float llr_for_one = 0;
            for (int j = 0; j < constellation.size(); j++) {
                for (int k = 0; k < constellation[j].size(); k++) {
                    if (constellation[j][k].gray_str[bit] == 0) {
                        llr_for_zero += exp((float)pow(abs(constellation[j][k].const_place - in[i]), 2)/(float)(-2*noise*noise));
                    } else {
                        llr_for_one += exp((float)pow(abs(constellation[j][k].const_place - in[i]), 2)/(float)(-2*noise*noise));
                    } 
                }
            }
            std::cout << "(" << llr_for_zero << "," << llr_for_one << ") ";
            out[i*bits_per_sym + bit] = std::min((float)1000.0, std::max((float)-1000.0, log(llr_for_one/llr_for_zero)));
        }
    }
    std::cout << "\n";

    //Resizing the elements per process for data transfer of out vector
    for (int i = 0; i < qsize; i++) {
        size_of_proc_data[i] *= bits_per_sym;
        displ[i] *= bits_per_sym;
    }

    //Broadcasting data to all other processes 
    MPI_Barrier(MPI_COMM_WORLD);
    if (qrank == 0) {
        MPI_Gatherv(MPI_IN_PLACE, size_of_proc_data[qrank], MPI_FLOAT, (void *)&out[displ[qrank]], size_of_proc_data, displ, MPI_FLOAT, 0, MPI_COMM_WORLD);
    } else {
        MPI_Gatherv((void *)&out[displ[qrank]], size_of_proc_data[qrank], MPI_FLOAT, (void *)&out[displ[qrank]], size_of_proc_data, displ, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast((void *)&out[0], (int)out.size(), MPI_FLOAT, 0, MPI_COMM_WORLD);
    //Freeing allocated memory and resizing input vector
    free(size_of_proc_data);
    free(displ);
}

//Sets the constellation values in a 2D vector of complex floats (only works if bits per symbol is even, Eg: 16, 64, 256, 1024, etc.)
void qam_llr_mpi::set_contellation(int _qam_size) {
    qam_size = _qam_size;
    printf("QAM size: %d\n", qam_size);
    bits_per_sym = ceil(log2(_qam_size));
    printf("Bits per sym: %d\n", bits_per_sym);
    int points_per_side = (int)ceil(sqrt(_qam_size));
    printf("Points per side: %d\n", points_per_side);

    //For BPSK modulation
    if (_qam_size == 2) {
        printf("QAM size is 2...\n");
        constellation.resize(1);
        constellation[0].resize(2);
        constellation[0][0].const_place = std::complex<float>(-1,0);
        constellation[0][0].gray_str.push_back(0);
        constellation[0][1].const_place = std::complex<float>(1,0);
        constellation[0][1].gray_str.push_back(1);
        printf("Here...\n");
    } else {
        printf("Here with more than 2...\n");
        //Returns if qam size cannot be created
        if (fmod(sqrt((float)_qam_size), 1) != 0) {
            printf("Cannot create qam constellation if bits per symbol is odd...\n");
            return;
        }

        //For higher than 2 order modulation

        //Resizing qam constellation
        constellation.resize(points_per_side);
        for (int i = 0; i < constellation.size(); i++) {
            constellation[i].resize(points_per_side);
        }
        //Starting qam constellation creation
        float row_val = -1*(points_per_side - 1);
        int gray_val;
        for (int i = 0; i < constellation.size(); i++) {
            float col_val = -1*(points_per_side - 1);
            
            for (int j = 0; j < constellation[i].size(); j++) {
                //Setting constellation index based on decimal conversion of gray values
                int temp = i*constellation.size() + j, idx, prev;
                idx = temp ^ (temp >> 1);
                /*
                for (int k = 0; k < bits_per_sym; k++) {
                    gray_val = temp % 2;
                    temp = (int)floor((float)temp/(float)2);
                    if (k == 0) {
                        idx = gray_val * pow(2, bits_per_sym - 1 - k);
                        prev = gray_val;
                    } else {
                        idx += ((gray_val + prev) % 2) * pow(2, bits_per_sym - 1 - k);
                        prev = (gray_val + prev) % 2;
                    }
                }
                */
                //Setting index values
                int ii = (int)floor((float)idx/(float)points_per_side), jj = idx % points_per_side;
                //Converting from decimal to gray codes for storage in consellation
                constellation[ii][jj].gray_str.resize(bits_per_sym);
                //int prev;
                temp = idx;
                for (int k = 0; k < bits_per_sym; k++) {
                    gray_val = temp % 2;
                    temp = (int)floor((float)temp/(float)2);
                //  if (k == 0) {
                        constellation[ii][jj].gray_str[bits_per_sym - 1 - k] = gray_val;
                //  } else {
                //      constellation[ii][jj].gray_str[k] = (prev + gray_val) % 2;
                //  }
                //  prev = gray_val;
                    //temp = (int)floor((float)temp/(float)2);
                }

                //Setting value of constellation point
                constellation[ii][jj].const_place = std::complex<float>(row_val, col_val);
                col_val += 2;

                //Normalizing the power of qam constellation points
                constellation[ii][jj].const_place *= (float)1/(float)(sqrt(2)*(points_per_side - 1));
        //      printf("Index: %d, Row: %d, Cols: %d, Row val: %f, Col val: %f\n", idx, ii, jj, std::real(constellation[ii][jj].const_place), std::imag(constellation[ii][jj].const_place));
        //      for (int k = 0; k < bits_per_sym; k++) {
        //          printf("%d ", constellation[ii][jj].gray_str[k]);
        //      }
        //      printf("\n");
                for (int z = 0 ;z < constellation[ii][jj].gray_str.size(); z++) {
                    std::cout << constellation[ii][jj].gray_str[z] << " ";
                }
                std::cout << "\n";
                std::cout << constellation[ii][jj].const_place << " ";

            }
            std::cout << std::endl;
            row_val += 2;
        }
    }
}

//Gives qam constellation values
std::vector<std::vector<std::complex<float>>> qam_llr_mpi::get_constellation_vals() {
    std::vector<std::vector<std::complex<float>>> vals;
    vals.resize(constellation.size());
    for (int i = 0; i < vals.size(); i++) {
        vals[i].resize(constellation[i].size());
        for (int j = 0; j < vals[i].size(); i++) {
            vals[i][j] = constellation[i][j].const_place;
        }
    }
    return vals;
}