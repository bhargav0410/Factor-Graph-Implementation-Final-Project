#include "ofdm_mpi.h"

ofdm_mpi::ofdm_mpi(int _rank, int _size) : orank(_rank), osize(_size) {}

ofdm_mpi::ofdm_mpi(int _rank, int _size, int _fft_size, int _prefix_size, int _num_ants) : orank(_rank), osize(_size), fft_size(_fft_size), prefix_size(_prefix_size), num_ants(_num_ants) {}

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
// ****************************************************************


//FFT of one row using OpenMP + MPI
void ofdm_mpi::ifft_omp_mpi(std::complex<float> *fft_in, std::complex<float> *fft_out, int num_omp_threads) {
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
    plan = fftwf_plan_dft_1d(fft_size, (fftwf_complex *)fft_in, (fftwf_complex *)fft_out, FFTW_REVERSE, FFTW_MEASURE /*FFTW_ESTIMATE*/);
    fftwf_execute(plan);
    //Destroying plan
    fftwf_destroy_plan(plan);
}

//Swaps the [0:FFT/2-1] and [FFT/2:FFT-1] halves of OFDM symbols
void ofdm_mpi::swap_halves(std::complex<float> *vec, int num_threads) {
    std::complex<float> temp;
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < fft_size/2; i++) {
        temp = vec[i];
        vec[i] = vec[i + fft_size/2];
        vec[i + fft_size/2] = temp;
    }
}

//Performs element by element division of complex vectors and stores answer in numerator
void ofdm_mpi::divide_by_vec_omp(std::complex<float> *numer, std::complex<float> *denom, std::complex<float> *out, int len, int num_threads) {
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < len, i++) {
        out[i] = numer[i]/denom[i];
    }

}

//Creates a pilot vector (currently a vector of ones are sent)
//If pilot vector is already created but IFFT is to be performed, this function can be used
void ofdm_mpi::create_pilot_vector_mpi(int num_tx_ants) {
    std::vector<int> size_of_proc_data(osize), displ(osize);
    if (pilot.size() != num_tx_ants) {
        pilot.resize(num_tx_ants);
    }
    //Distributing pilot creation among servers
    load_balancing_mpi(&size_of_proc_data[0], &displ[0], osize, (int)pilot.size());
    //Defining threads to be given for each task
    int threads_for_task = (int)ceil((float)NUM_THREADS/(float)size_of_proc_data[orank]);

    #pragma omp parallel for num_threads(size_of_proc_data[orank])
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
}

//Performs least squares channel estimation (assumes the pilot vector is in [0 FFT/2 : -FFT/2 -1] form, if pilot vector not in such form then conversion is necessary)
void ofdm_mpi::chan_est_update_mpi(std::vector<std::vector<std::complex<float>>> &chan_est_in) {
    std::vector<int> size_of_proc_data(osize), displ(osize);
    load_balancing_mpi(&size_of_proc_data[0], &displ[0], osize, num_ants);

    chan_est.resize(chan_est_in.size());
    int threads_for_task = (int)ceil((float)NUM_THREADS/(float)size_of_proc_data[orank]);
    #pragma omp parallel for num_threads(size_of_proc_data[orank])
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
        std::rotate(chan_est[i].begin(), chan_est[i].begin() + 1, chan_est[i].end());
        chan_est[i].resize(fft_size - 1);
        
        //Dividing by pilot vector
        divide_by_vec_omp(&chan_est[i][0], &pilot[i][prefix_size + 1], &chan_est[i][0], fft_size - 1, threads_for_task);
    }
}

//MRC is performed to combine data of all receiving antennas at the MIMO system
//Data is combined and stored in processor with index 0
void ofdm_mpi::maximal_ratio_combining_mpi(std::vector<std::vector<std::complex<float>>> &in, std::vector<std::complex<float>> &out) {
    std::vector<int> size_of_proc_data(osize), displ(osize);
    load_balancing_mpi(&size_of_proc_data[0], &displ[0], osize, num_ants);

    if (in.size() != num_ants) {
        printf("Number of antennas set in object and input vector size do not match...\n");
    }

    int threads_for_task = (int)ceil((float)NUM_THREADS/(float)size_of_proc_data[orank]);
    #pragma omp parallel for num_threads(size_of_proc_data[orank])
    for (int i = displ[orank]; i < displ[orank] + size_of_proc_data[orank]; i++) {
        out[0].resize(fft_size);
 
        //Performing FFT for conversion to frequency domain
        fft_omp_mpi(&in[i][prefix_size], &out[i][0], threads_for_task);
        
        //Swapping halves
        swap_halves(&out[i][0], threads_for_task);
        
        //Removing DC sub-carrier and resizing
        std::rotate(out[i].begin(), out[i].begin() + 1, out[i].end());
        out[i].resize(fft_size - 1);
        
        //Dividing by pilot vector
        divide_by_vec_omp(&out[i][0], &chan_est[i][0], &out[i][0], fft_size - 1, threads_for_task);

        //Combining local output
        #pragma omp parallel for num_threads(threads_for_task)
        for (int s = 0; s < fft_size - 1; s++) {
            #pragma omp atomic
            out[0][s] += out[i][s];
        }
    }
    //Combining output of all processors
    for (int s = 0; s < fft_size - 1; s++) {
        MPI_Reduce(&out[0][s], &out[0][s], 1, MPI_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
    }
}

//MRT is performed on input vector to be transmitted via multiple antennas assuming channel estimate for receiver is known
//Assumes channel estimation vector is of the form [1 FFT/2 : -FFT/2 -1], also input vector should be of size (fft_size - 1)
void ofdm_mpi::maximal_ratio_transmission_mpi(std::vector<std::complex<float>> &in, std::vector<std::vector<std::complex<float>>> &out) {
    std::vector<int> size_of_proc_data(osize), displ(osize);
    //Distributing pilot creation among servers
    load_balancing_mpi(&size_of_proc_data[0], &displ[0], osize, num_ants);
    //Defining threads to be given for each task
    int threads_for_task = (int)ceil((float)NUM_THREADS/(float)size_of_proc_data[orank]);
    out.resize(num_ants);

    //Checking input size
    if (in[0].size() != fft_size - 1) {
        printf("Incorrect input size...Input size should be %d\n", fft_size - 1);    
    }
    //Adding zero to DC sub-carrier of input vector
    in.insert(in[0].begin() + fft_size/2, std::complex<float>(0,0));

    //Swapping vector halves of input vector
    swap_halves(&in[0][prefix_size], threads_for_task);

    #pragma omp parallel for num_threads(size_of_proc_data[orank])
    for (int i = displ[orank]; i < displ[orank] + size_of_proc_data[orank]; i++) {
        //Resizing output vector
        out[i].resize(fft_size + prefix_size);

        //Adding 0 to DC sub-carrier of output vector
        out[i][prefix_size] = 0;

        //Dividing by channel estimation vector
        divide_by_vec_omp(&in[0][1], &chan_est[i][0], &out[i][prefix_size + 1], fft_size - 1, threads_for_task);

        //IFFT for conversion to time domain
        ifft_omp_mpi(&out[i][prefix_size], &out[i][prefix_size], threads_for_task);

        //Adding cyclic prefix
        for (int j = 0; j < prefix_size; j++) {
            out[i][j] = out[i][fft_size + j];
        }
    }
}

//void ofdm_mpi::mod_one_frame_mpi(std::vector<std::vector<std::complex<float>>> &, std::vector<std::vector<std::complex<float>>> &);

//void ofdm_mpi::demod_one_frame_mpi(std::vector<std::vector<std::complex<float>>> &, std::vector<std::vector<std::complex<float>>> &);

/*
int ofdm_mpi::get_fft_size();
void ofdm_mpi::set_fft_size(int);
int ofdm_mpi::get_prefix_size();
void ofdm_mpi::set_prefix_size(int);
void ofdm_mpi::set_chan_est(std::vector<std::vector<std::complex<float>>>);
std::vector<std::vector<std::complex<float>>> ofdm_mpi::get_chan_est();
void ofdm_mpi::set_pilot(std::vector<std::vector<std::complex<float>>>);
std::vector<std::vector<std::complex<float>>> ofdm_mpi::get_pilot();
*/