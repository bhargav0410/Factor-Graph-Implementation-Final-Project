#include "ofdm_mpi.h"

ofdm_mpi::ofdm_mpi(int _rank, int _size) : orank(_rank), osize(_size) {
    fft_size = 64;
    prefix_size = 16;
    num_ants = 1;
    std::vector<std::complex<float>> fft_in(fft_size, std::complex<float> (1, 0)), ifft_in(fft_size, std::complex<float> (1, 0));
    int num_threads = NUM_THREADS;
    fft_omp_mpi(&fft_in[0], &fft_in[0], num_threads);
    ifft_omp_mpi(&ifft_in[0], &ifft_in[0], num_threads);
}

ofdm_mpi::ofdm_mpi(int _rank, int _size, int _fft_size, int _prefix_size, int _num_ants) : orank(_rank), osize(_size), fft_size(_fft_size), prefix_size(_prefix_size), num_ants(_num_ants) {
    std::vector<std::complex<float>> fft_in(fft_size, std::complex<float> (1, 0)), ifft_in(fft_size, std::complex<float> (1, 0));
    int num_threads = NUM_THREADS;
    fft_omp_mpi(&fft_in[0], &fft_in[0], num_threads);
    ifft_omp_mpi(&ifft_in[0], &ifft_in[0], num_threads);
}

ofdm_mpi::~ofdm_mpi() {}

//Distributed the tasks amongst workers in a balanced fashion
void ofdm_mpi::load_balancing_mpi(int *size_of_proc_data, int *displ, int num_procs, int len) {
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

// ************************ Only using MPI *************************
/*
//FFT of one row
void ofdm_mpi::fft_mpi(std::vector<std::complex<float>> *fft_vec) {
    fftwf_plan plan;
    ptrdiff_t alloc_local, local_ni, local_i_start, local_no, local_o_start;
    fftwf_mpi_init();
    //Getting local data size
    alloc_local = fftwf_mpi_local_size_1d((ptrdiff_t)fft_size, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE, &local_ni, &local_i_start, &local_no, &local_o_start);
    //Planning and executing FFT
    plan = fftwf_mpi_plan_dft_1d((ptrdiff_t)fft_size, (fftwf_complex *)&fft_vec[(int)local_i_start], (fftwf_complex *)&fft_vec[(int)local_o_start], MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE);
    fftwf_execute(plan);
    //Destroying plan
    fftwf_destroy_plan(plan);
}

//IFFT of one row
void ofdm_mpi::ifft_mpi(std::vector<std::complex<float>> *fft_vec) {
    fftwf_plan plan;
    ptrdiff_t alloc_local, local_ni, local_i_start, local_no, local_o_start;
    fftwf_mpi_init();
    //Getting local data size
    alloc_local = fftwf_mpi_local_size_1d((ptrdiff_t)fft_size, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE, &local_ni, &local_i_start, &local_no, &local_o_start);
    //Planning and executing IFFT
    plan = fftwf_mpi_plan_dft_1d((ptrdiff_t)fft_size, (fftwf_complex *)&fft_vec[(int)local_i_start], (fftwf_complex *)&fft_vec[(int)local_o_start], MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE);
    fftwf_execute(plan);
    //Destroying plan
    fftwf_destroy_plan(plan);
}
*/
// ****************************************************************


//FFT of one row using OpenMP + MPI
void ofdm_mpi::fft_omp_mpi(std::complex<float> *fft_in, std::complex<float> *fft_out, int num_omp_threads) {
    if (fftwf_init_threads() == 0) {
        printf("Error...\n");
        return;
    }
    fftwf_plan_with_nthreads(num_omp_threads);
    fftwf_plan plan;
    
    //Spreading the FFT to be performed over multiple threads using OpenMP
        //Setting up plan to execute
    plan = fftwf_plan_dft_1d(fft_size, (fftwf_complex *)fft_in, (fftwf_complex *)fft_out, FFTW_FORWARD, FFTW_MEASURE /*FFTW_ESTIMATE*/);
    fftwf_execute(plan);
    //Destroying plan
    fftwf_destroy_plan(plan);
}

//IFFT of one row using OpenMP + MPI
void ofdm_mpi::ifft_omp_mpi(std::complex<float> *fft_in, std::complex<float> *fft_out, int num_omp_threads) {
    if (fftwf_init_threads() == 0) {
        printf("Error...\n");
        return;
    }
    //Spreading the FFT to be performed over multiple threads using OpenMP
    fftwf_plan_with_nthreads(num_omp_threads);
    //Setting up plan to execute
    fftwf_plan plan;
    plan = fftwf_plan_dft_1d(fft_size, (fftwf_complex *)fft_in, (fftwf_complex *)fft_out, FFTW_BACKWARD, FFTW_MEASURE /*FFTW_ESTIMATE*/);
    fftwf_execute(plan);
    //Destroying plan
    fftwf_destroy_plan(plan);
}

//Finding maximum absolute value within vector
float ofdm_mpi::find_max_val(std::complex<float> *in_vec, int len, int threads) {
    std::vector<float> abs_vec(len);
    float temp_ret;

    //Getting absolute value of complex number
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < len; i++) {
        abs_vec[i] = std::abs(in_vec[i]);
    }

    //Finding max value of absolute value vector
    for (int step = 1; step < len; step *= 2) {

        #pragma omp parallel for num_threads(threads)
        for (int i = 0; i < len; i += step*2) {
            if (i + step < len) {
                abs_vec[i] = std::max(abs_vec[i], abs_vec[i+step]);
            }
        }
    }

    return abs_vec[0];
}

//Swaps the [0:FFT/2-1] and [FFT/2:FFT-1] halves of OFDM symbols
void ofdm_mpi::swap_halves(std::complex<float> *vec, int num_threads) {
    std::vector<std::complex<float>> temp((int)ceil((float)fft_size/(float)2));
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < fft_size/2; i++) {
        temp[i] = vec[i];
        vec[i] = vec[i + fft_size/2];
        vec[i + fft_size/2] = temp[i];
    }
}

//Performs element by element division of complex vectors and stores answer in numerator
void ofdm_mpi::divide_by_vec_omp(std::complex<float> *numer, std::complex<float> *denom, std::complex<float> *out, int len, int _num_threads) {
    #pragma omp parallel for num_threads(_num_threads)
	for (int i = 0; i < len; i++) {
        out[i] = numer[i]/denom[i];
    }
}

//Performs element by element division of complex vectors and stores answer in numerator
void ofdm_mpi::mult_by_conj_omp(std::complex<float> *in_vec, std::complex<float> *conj_vec, std::complex<float> *out, int len, int _num_threads) {
    #pragma omp parallel for num_threads(_num_threads)
	for (int i = 0; i < len; i++) {
        out[i] = in_vec[i] * std::conj(conj_vec[i]);
    }
}


//Modulates a pilot vector (currently a vector of ones are sent)
//If pilot vector is already created but IFFT is to be performed, this function can be used
void ofdm_mpi::mod_pilot_vector_mpi() {
    int num_tx_ants = num_ants;
    std::vector<int> size_of_proc_data(osize), displ(osize);
    if (pilot.size() != num_tx_ants) {
        pilot.resize(num_tx_ants);
    }
    //Distributing pilot creation among servers
    load_balancing_mpi(&size_of_proc_data[0], &displ[0], osize, (int)pilot.size());
    //Defining threads to be given for each task
    int threads_for_task = NUM_THREADS; //(int)ceil((float)NUM_THREADS/(float)size_of_proc_data[orank]);

    //#pragma omp parallel for num_threads(size_of_proc_data[orank])
    for (int i = displ[orank]; i < displ[orank] + size_of_proc_data[orank]; i++) {
        if (pilot[i].size() == 0) {
            pilot[i].resize(fft_size + prefix_size);
            for (int j = prefix_size; j < fft_size + prefix_size; j++) {
                //Adding pilot symbols (except for DC sub-carrier)
                if (j - prefix_size != fft_size/2) {
                    pilot[i][j] = std::complex<float>(1, 0);
                } else {
                    pilot[i][j] = std::complex<float>(0, 0);
                }
            }
        } else if (pilot[i].size() < fft_size + prefix_size) {
            pilot[i].resize(fft_size + prefix_size);
            for (int j = fft_size - 1; j >= 0; j--) {
                pilot[i][j+prefix_size] = pilot[i][j];
            }
        }
        //Swapping vector halves
        swap_halves(&pilot[i][prefix_size], threads_for_task);

        //IFFT for conversion to time domain
        ifft_omp_mpi(&pilot[i][prefix_size], &pilot[i][prefix_size], threads_for_task);

        //Adding cyclic prefix
        for (int j = 0; j < prefix_size; j++) {
            pilot[i][j] = pilot[i][fft_size + j];
        }
    }
    pilot_vec_mod = true;
}

//Performs least squares channel estimation (assumes the pilot vector is in [0 FFT/2 : -FFT/2 -1] form, if pilot vector not in such form then conversion is necessary)
void ofdm_mpi::chan_est_update_mpi(std::vector<std::vector<std::complex<float>>> &chan_est_in, std::vector<std::complex<float>> &pilot) {
    std::vector<int> size_of_proc_data(osize), displ(osize);
    load_balancing_mpi(&size_of_proc_data[0], &displ[0], osize, num_ants);

    chan_est.clear();
    chan_est.resize(chan_est_in.size());
    chan_est_abs_sqrd.clear();
    chan_est_abs_sqrd.resize(fft_size - 1, std::complex<float> (0,0));
    int threads_for_task = NUM_THREADS; //(int)ceil((float)NUM_THREADS/(float)size_of_proc_data[orank]);
    //#pragma omp parallel for num_threads(size_of_proc_data[orank])
    for (int i = displ[orank]; i < displ[orank] + size_of_proc_data[orank]; i++) {
        chan_est[i].resize(fft_size);
        //Removing cyclic prefix
       // for (int j = 0; j < fft_size; j++) {
        //    chan_est[i][j] = chan_est[i][prefix_size + j];
       // }
        //Performing FFT for conversion to frequency domain
        fft_omp_mpi(&chan_est_in[i][prefix_size], &chan_est[i][0], threads_for_task);
        
        //Swapping halves
        swap_halves(&chan_est[i][0], threads_for_task);
        
        //Removing DC sub-carrier and resizing
        //std::rotate(chan_est[i].begin(), chan_est[i].begin() + 1, chan_est[i].end());
        //chan_est[i].resize(fft_size - 1);
        chan_est[i].erase(chan_est[i].begin() + (fft_size/2));
        
        //Dividing by pilot vector
        divide_by_vec_omp(&chan_est[i][0], &pilot[0], &chan_est[i][0], fft_size - 1, threads_for_task);

        //Creating the conjugate vector
        std::vector<std::complex<float>> temp_vec(fft_size - 1);
        mult_by_conj_omp(&chan_est[i][0], &chan_est[i][0], &temp_vec[0], fft_size - 1, threads_for_task);
        
        //Making the denominator to be used for MRC and MRT
        #pragma omp parallel for num_threads(NUM_THREADS)
        for (int j = 0; j < fft_size - 1; j++) {
            chan_est_abs_sqrd[j] += temp_vec[j];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //Combining output of all processors
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int s = 0; s < fft_size - 1; s++) {
        if (orank == 0) {
            MPI_Reduce(MPI_IN_PLACE, (void *)&chan_est_abs_sqrd[s], 1, MPI_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
        } else {
            MPI_Reduce((void *)&chan_est_abs_sqrd[s], (void *)&chan_est_abs_sqrd[s], 1, MPI_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
        }
    }
    MPI_Bcast((void *)&chan_est_abs_sqrd[0], (int)chan_est_abs_sqrd.size(), MPI_COMPLEX, 0, MPI_COMM_WORLD);
    
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int i = displ[orank]; i < displ[orank] + size_of_proc_data[orank]; i++) {
        if (orank == 0) {
            MPI_Gather(MPI_IN_PLACE, fft_size - 1, MPI_COMPLEX, (void *)&chan_est[i][0], fft_size - 1, MPI_COMPLEX, 0, MPI_COMM_WORLD);
        } else {
            MPI_Gather((void *)&chan_est[i][0], fft_size - 1, MPI_COMPLEX, (void *)&chan_est[i][0], fft_size - 1, MPI_COMPLEX, 0, MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int i = 0; i < num_ants; i++) {
        MPI_Bcast((void *)&chan_est[i][0], fft_size - 1, MPI_COMPLEX, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < num_ants; i++) {
        for (int j = 0; j < chan_est_in[i].size(); j++) {
            std::cout << chan_est_in[i][j] << " ";
        }
        std::cout << "\n";
    }

    for (int i = 0; i < num_ants; i++) {
        for (int j = 0; j < chan_est[i].size(); j++) {
            std::cout << chan_est[i][j] << " ";
        }
        std::cout << "\n";
    }

    for (int j = 0; j < chan_est_abs_sqrd.size(); j++) {
        std::cout << chan_est_abs_sqrd[j] << " ";
    }
    std::cout << "\n";

}



//MRC is performed to combine data of all receiving antennas at the MIMO system
//Data is combined and stored in processor with index 0
void ofdm_mpi::maximal_ratio_combining_mpi(std::vector<std::vector<std::complex<float>>> &in, std::vector<std::complex<float>> &out) {
    std::vector<int> size_of_proc_data(osize), displ(osize);
    load_balancing_mpi(&size_of_proc_data[0], &displ[0], osize, num_ants);

    if (in.size() != num_ants) {
        printf("Number of antennas set in object and input vector size do not match...\n");
    }

	std::vector<std::vector<std::complex<float>>> out1;
    out1.resize((int)in.size());
	out.resize(fft_size - 1);
    std::vector<std::complex<float>> out_temp(fft_size - 1, 0);
    int threads_for_task = NUM_THREADS; //(int)ceil((float)NUM_THREADS/(float)size_of_proc_data[orank]);
    //#pragma omp parallel for num_threads(size_of_proc_data[orank])
    for (int i = displ[orank]; i < displ[orank] + size_of_proc_data[orank]; i++) {
        out1[i].resize(fft_size);

        //Performing FFT for conversion to frequency domain
        fft_omp_mpi(&in[i][prefix_size], &out1[i][0], threads_for_task);
        
        //Swapping halves
        swap_halves(&out1[i][0], threads_for_task);

        for (int j = 0; j < out1[i].size(); j++) {
            std::cout << out1[i][j] << " ";
        }
        std::cout << "\n";
        
        //Removing DC sub-carrier and resizing
        //std::rotate(out1[i].begin(), out1[i].begin() + 1, out1[i].end());
        //out1[i].resize(fft_size - 1);
        out1[i].erase(out1[i].begin() + (fft_size/2));

        for (int j = 0; j < out1[i].size(); j++) {
            std::cout << out1[i][j] << " ";
        }
        std::cout << "\n";
        
        //Performing Y' = H^H * Y
        mult_by_conj_omp(&out1[i][0], &chan_est[i][0], &out1[i][0], fft_size - 1, threads_for_task);

        //Performing Y' / H^H * H
        divide_by_vec_omp(&out1[i][0], &chan_est_abs_sqrd[0], &out1[i][0], fft_size - 1, threads_for_task);

        //Combining local output
        #pragma omp parallel for num_threads(NUM_THREADS)
        for (int s = 0; s < fft_size - 1; s++) {
            out[s] += out1[i][s];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //Combining output of all processors
    //#pragma omp parallel for num_threads(NUM_THREADS)
    for (int s = 0; s < fft_size - 1; s++) {
        if (orank == 0) {
            MPI_Reduce(MPI_IN_PLACE, (void *)&out[s], 1, MPI_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
        } else {
            MPI_Reduce((void *)&out[s], (void *)&out[s], 1, MPI_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
        }
    }
    MPI_Bcast((void *)&out[0], (int)out.size(), MPI_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}

//MRT is performed on input vector to be transmitted via multiple antennas assuming channel estimate for receiver is known
//Assumes channel estimation vector is of the form [1 FFT/2 : -FFT/2 -1], also input vector should be of size (fft_size - 1)
void ofdm_mpi::maximal_ratio_transmission_mpi(std::vector<std::complex<float>> &in, std::vector<std::vector<std::complex<float>>> &out) {
    std::vector<int> size_of_proc_data(osize), displ(osize);
    //Distributing pilot creation among servers
    load_balancing_mpi(&size_of_proc_data[0], &displ[0], osize, num_ants);
    //std::cout << "Size for proc "<< orank << " data: " << size_of_proc_data[orank] << '\n';
    //std::cout << "From " << displ[orank] << " to " << size_of_proc_data[orank] + displ[orank] << "\n";
    //std::cout << "Input size for MRT: " << in.size() << '\n';
    //Defining threads to be given for each task
    int threads_for_task = NUM_THREADS; //(int)ceil((float)NUM_THREADS/(float)size_of_proc_data[orank]);
    //std::cout << "Threads for task: " << threads_for_task << "\n";
    out.resize(num_ants);

    //Checking input size
    if (in.size() != (fft_size - 1)) {
        printf("Incorrect input size %d...Input size should be %d\n", (int)in.size(), fft_size - 1);    
    }

    std::vector<std::complex<float>> in_temp = in;
    
    //Adding zero to DC sub-carrier of input vector
    in_temp.insert(in_temp.begin() + fft_size/2, std::complex<float>(0,0));

    //Swapping vector halves of input vector
    swap_halves(&in_temp[0], threads_for_task);

    
    std::cout << "Input size for MRT: " << in_temp.size() << '\n';
    for (int i = 0; i < in_temp.size(); i++) {
        std::cout << in_temp[i] << " ";
    }
    std::cout << "\n";
    

    //#pragma omp parallel for num_threads(size_of_proc_data[orank])
    for (int i = displ[orank]; i < displ[orank] + size_of_proc_data[orank]; i++) {
        //Resizing output vector
        out[i].resize(fft_size + prefix_size);

        //Adding 0 to DC sub-carrier of output vector
        //out[i][prefix_size] = 0;

        //Dividing by channel estimation vector
        //divide_by_vec_omp(&in_temp[1], &chan_est[i][0], &out[i][prefix_size + 1], fft_size - 1, threads_for_task);

        //Performing Y' = H^H * Y
        mult_by_conj_omp(&in_temp[1], &chan_est[i][0], &out[i][prefix_size + 1], fft_size - 1, threads_for_task);

        //Performing Y' / H^H * H
        divide_by_vec_omp(&out[i][prefix_size + 1], &chan_est_abs_sqrd[0], &out[i][prefix_size + 1], fft_size - 1, threads_for_task);

        //IFFT for conversion to time domain
        ifft_omp_mpi(&out[i][prefix_size], &out[i][prefix_size], threads_for_task);

        //Adding cyclic prefix
        #pragma omp parallel for num_threads(NUM_THREADS)
        for (int j = 0; j < prefix_size; j++) {
            out[i][j] = out[i][fft_size + j];
        }
        /*
        //Finding maximum absolute value
        float abs_val = find_max_val(&out[i][0], (int)out[i].size(), threads_for_task);

        //Dividing by max absolute value for each antenna
        #pragma omp parallel for num_threads(NUM_THREADS)
        for (int j = 0; j < out[i].size(); j++) {
            out[i][j] = out[i][j]/abs_val;
        }
        */
    }
}

//Creates frame to be used for chanenl sounding transmission
void ofdm_mpi::create_chan_sound_frame(std::vector<std::complex<float>> &pilot_in, std::vector<std::vector<std::complex<float>>> &out) {
    //Setting pilot vector for modulation
    this->set_pilot(pilot_in);
    this->mod_pilot_vector_mpi();

    //Resizing output vector and adding pilot vectors
    out.resize(num_ants);
    for (int i = 0; i < num_ants; i++) {
        out[i].resize((fft_size + prefix_size)*num_ants);
        for (int j = 0; j < fft_size + prefix_size; j++) {
            out[i][j + i*(fft_size + prefix_size)] = pilot[i][j];
        }
    }
}

//Creates frame with LS channel estimate vector for CSI feedback
void ofdm_mpi::give_csi_fb(std::vector<std::complex<float>> &in, std::vector<std::complex<float>> &out, std::vector<std::complex<float>> &pilot_in) {
    //Getting number of antennas by using FFT size and prefix size
    int temp_num_ants = in.size()/(fft_size + prefix_size);
    
    std::vector<int> size_of_proc_data(osize), displ(osize);
    //Distributing pilot creation among servers
    load_balancing_mpi(&size_of_proc_data[0], &displ[0], osize, temp_num_ants);
    //Defining threads to be given for each task
    int threads_for_task = (int)ceil((float)NUM_THREADS/(float)size_of_proc_data[orank]);

    //resizing output vector based on number of antennas the transmisster has
    out.resize((fft_size-1)*temp_num_ants);

    std::vector<std::vector<std::complex<float>>> temp_out;
    temp_out.resize(temp_num_ants);
    #pragma omp parallel for num_threads(size_of_proc_data[orank])
    for (int i = displ[orank]; i < displ[orank] + size_of_proc_data[orank]; i++) {
        temp_out[i].resize(fft_size);
        //Performing FFT for conversion to frequency domain
        fft_omp_mpi(&in[prefix_size + i*(prefix_size + fft_size)], &temp_out[i][0], threads_for_task);
        
        //Swapping halves
        swap_halves(&temp_out[i][0], threads_for_task);
        
        //Removing DC sub-carrier and resizing
        std::rotate(temp_out[i].begin(), temp_out[i].begin() + 1, temp_out[i].end());
        temp_out[i].resize(fft_size - 1);
        
        //Dividing by pilot vector
        divide_by_vec_omp(&temp_out[i][0], &pilot_in[0], &out[i*(fft_size-1)], fft_size - 1, threads_for_task);
    }
    //Gathering CSi from all processors (MPI based distributed procs) and boradcasting gathered data back to all
    if (orank == 0) {
        MPI_Gatherv(MPI_IN_PLACE, size_of_proc_data[orank], MPI_COMPLEX, (void *)&out[displ[orank]*(fft_size - 1)], &size_of_proc_data[0], &displ[0], MPI_COMPLEX, 0, MPI_COMM_WORLD);
    } else {
        MPI_Gatherv((void *)&out[displ[orank]*(fft_size - 1)], size_of_proc_data[orank], MPI_COMPLEX, (void *)&out[displ[orank]*(fft_size - 1)], &size_of_proc_data[0], &displ[0], MPI_COMPLEX, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast((void *)&out[0], (int)out.size(), MPI_COMPLEX, 0, MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);
}

//OFDM modulation of one frame (perform MRT for each symbol including the pilot vector) (assumes channel estimation vector is already set)
void ofdm_mpi::mod_one_frame_mpi(std::vector<std::complex<float>> &in, std::vector<std::vector<std::complex<float>>> &out, std::vector<std::complex<float>> &unmod_pilot_vec) {
	std::vector<int> size_of_proc_data(osize), displ(osize);
    //Distributing pilot creation among servers
    load_balancing_mpi(&size_of_proc_data[0], &displ[0], osize, num_ants);
    
    //Calculating input vector length and adding zeros if necessary
    int total_syms = (int)ceil(((float)in.size())/((float)(fft_size - 1)));
    int zeros_to_add = (fft_size - 1) - (((int)in.size()) % (fft_size - 1));
    in.resize(in.size() + zeros_to_add);
    std::cout << "Input vector size: " << in.size() << "\n";
    std::cout << "Total number of syms: " << total_syms << "\n";
    std::cout << "Pilot vector size: " << unmod_pilot_vec.size() << "\n";

    // Clearing output vector
    out.resize(num_ants);

    //std::vector<std::complex<float>> temp_vec(fft_size, std::complex<float> (1,0));
    //ifft_omp_mpi(&temp_vec[0], &temp_vec[0], NUM_THREADS);


    printf("Modulating pilot vector...\n");
    //Modulating pilot vector if unmodulated pilot vector is given as input
    if (unmod_pilot_vec.size() != 0) {
        if (unmod_pilot_vec.size() > (fft_size - 1)) {
            unmod_pilot_vec.resize(fft_size - 1);
        } else if (unmod_pilot_vec.size() < fft_size - 1) {
            std::cout << "Cannot modulate pilot vector. Please make pilot vector size (FFT size - 1)...\n";
            return;
        }
        //Performing MRT on pilot symbols
        std::cout << "Performing MRT on pilot vector...\n";
        std::vector<std::complex<float>> temp_pilot = unmod_pilot_vec;
        maximal_ratio_transmission_mpi(temp_pilot, pilot);
        this->pilot_vec_mod == true;
    } else {
        std::cout << "Pilot not given. Creating a pilot of all 1s...\n";
        //Modulating pilot vector into output if input pilot vector is not given
        std::vector<std::complex<float>> temp;
        if (pilot.size() != num_ants) {
            pilot.resize(num_ants);
            temp.resize(fft_size - 1);
            for (int j = 0; j < fft_size - 1; j++) {
                //Adding pilot symbols (except for DC sub-carrier)
                temp[j] = std::complex<float>(1, 0);
            }
        } else if (pilot.size() == num_ants) {
            temp.insert(temp.end(), pilot[0].begin(), pilot[0].end());
        }
        //Performing MRT on pilot symbols
        maximal_ratio_transmission_mpi(temp, pilot);
        this->pilot_vec_mod == true;
    }
    //Adding modulated pilot vector to output
    for (int i = 0; i < num_ants; i++) {
        for (int j = 0; j < pilot[i].size(); j++) {
            std::cout << pilot[i][j] << " ";
        }
        std::cout << "\n";

        out[i].insert(out[i].end(), pilot[i].begin(), pilot[i].end());
    }

    printf("Pilot vector modulation done.\n");

    //Modulating each symbol
    for (int syms = 0; syms < total_syms; syms++) {
        //Performing MRT on input symbols
        std::vector<std::complex<float>> temp_in(fft_size - 1);
        for (int i = 0; i < fft_size - 1; i++) {
            temp_in[i] = in[syms*(fft_size - 1) + i];
        }
        std::vector<std::vector<std::complex<float>>> temp_out;
        maximal_ratio_transmission_mpi(temp_in, temp_out);
        //Inserting MRT output in output vector
		for (int i = 0; i < num_ants; i++) {
			out[i].insert(out[i].end(), temp_out[i].begin(), temp_out[i].end());
		}
    }
    
    for (int i = displ[orank]; i < displ[orank] + size_of_proc_data[orank]; i++) {
        //Finding maimum absolute value
        float abs_val = /*fft_size - 1; */ find_max_val(&out[i][0], (int)out[i].size(), NUM_THREADS);

        //Dividing by max absolute value for each antenna
        #pragma omp parallel for num_threads(NUM_THREADS)
        for (int j = 0; j < out[i].size(); j++) {
            out[i][j] = out[i][j]/abs_val;
        }
    }

}

//OFDM demodulation of one frame including channel estimation
void ofdm_mpi::demod_one_frame_mpi(std::vector<std::vector<std::complex<float>>> &in, std::vector<std::complex<float>> &out, std::vector<std::complex<float>> &pilot_vec) {
    //Getting input length
    int num_syms = (int)in[0].size()/(fft_size+prefix_size);
    std::cout << "Num syms to demod: " << num_syms << "\n";
    
    //Performing channel sounding
    if (pilot_vec.size() != (fft_size - 1))  {
        std::cout << "Pilot vec size: " << pilot_vec.size() << "\n";
        std::cout << "Add pilot vector of size (FFT size - 1)...\n";
        return;
    }

    std::cout << "\n Pilot vec for demodulation:\n";
    for (int i = 0; i < pilot_vec.size(); i++) {
        std::cout << pilot_vec[i] << " ";
    }
    std::cout << "\n";

    //Demodulating frame
    std::vector<std::vector<std::complex<float>>> temp_in;
    temp_in.resize(num_ants);
    for (int i = 0; i < num_ants; i++) {
        temp_in[i].resize(fft_size + prefix_size);
    }
    std::vector<std::complex<float>> temp_out;
    //Performing MRC for each symbol
    for (int syms = 0; syms < num_syms; syms++) {
        for (int a = 0; a < num_ants; a++) {
            for (int i = 0; i < fft_size + prefix_size; i++) {
                temp_in[a][i] = in[a][syms*(fft_size + prefix_size) + i];
            }
        }
        if (syms == 0) {
            chan_est_update_mpi(temp_in, pilot_vec);
        } else {
            maximal_ratio_combining_mpi(temp_in, temp_out);
            out.insert(out.end(), temp_out.begin(), temp_out.end());
        }
        for (int i = 0; i < temp_out.size(); i++) {
            std::cout << temp_out[i] << " ";
        }
        std::cout << "\n";
    }
    
    std::cout << "Output vector size: " << out.size() << "\n";

    //Getting max value of output vector

}




// ******************** Get and set vectors for FFT size, Cyclc Prefix size, pilot vector, and channel estimation vector values *************
int ofdm_mpi::get_fft_size() {
    return fft_size;
}

void ofdm_mpi::set_fft_size(int _fft_size) {
    fft_size = _fft_size;
}
int ofdm_mpi::get_prefix_size() {
    return prefix_size;
}

void ofdm_mpi::set_prefix_size(int _prefix_size) {
    prefix_size = _prefix_size;
}

void ofdm_mpi::set_chan_est(std::vector<std::vector<std::complex<float>>> &in) {
    chan_est = in;
    chan_est_abs_sqrd.resize(chan_est[0].size());
    //Making the denominator to be used for MRC and MRT
    for (int i = 0; i < chan_est.size(); i++) {
        #pragma omp parallel for num_threads(NUM_THREADS)
        for (int j = 0; j < chan_est[i].size(); j++) {
            chan_est_abs_sqrd[j] += chan_est[i][j]*std::conj(chan_est[i][j]);
        }
    }
}
std::vector<std::vector<std::complex<float>>> ofdm_mpi::get_chan_est() {
    return chan_est;
}

//Set symbols of pilot vector (another function used for modulating it)
void ofdm_mpi::set_pilot(std::vector<std::complex<float>> &in) {
    if (in.size() == fft_size) {
        std::cout << "Removing last value as pilot size needs to be (FFT size - 1)...\n";
        in.resize(fft_size - 1);
    }
    if (in.size() > fft_size || in.size() < fft_size - 1) {
        std::cout << "Cannot be used for pilot vector. Unmodulated pilot vector size should be (FFT Size - 1)...\n";
        return;
    }
    pilot.resize(num_ants);
    for (int i = 0; i < num_ants; i++) {
        pilot[i].insert(pilot[i].end(), in.begin(), in.end());
    }
    pilot_vec_mod = false;
}
std::vector<std::vector<std::complex<float>>> ofdm_mpi::get_pilot() {
    return pilot;
}
//************************************************************************