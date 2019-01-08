#include "ofdm_mpi.h"

ofdm_mpi::ofdm_mpi(int _rank, int _size) : orank(_rank), osize(_size) {}

ofdm_mpi::ofdm_mpi(int _rank, int _size, int _fft_size, int _prefix_size) : orank(_rank), osize(_size), fft_size(_fft_size), prefix_size(_prefix_size) {}

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

//Creates a pilot vector (currently a vector of ones are sent)
void ofdm_mpi::create_pilot_vector_mpi(std::vector<std::vector<std::complex<float>>> &pilot) {
    std::vector<int> size_of_proc_data(osize), displ(osize);
    load_balancing_mpi(&size_of_proc_data[0], &displ[0], osize, fft_size + prefix_size);
    for (int i = 0; i < pilot.size(); i++) {
        for (int j = prefix_size; j < fft_size + prefix_size; j++) {
            //Adding pilot symbols
            pilot[i][j] = std::complex<float>(1, 0);
        }
        fft_mpi(&pilot[i][prefix_size]);
        for (int j = 0; j < prefix_size; j++) {
            pilot[i][j] = pilot[i][fft_size + j];
        }
    }
}

void ofdm_mpi::chan_est_update_mpi(std::vector<std::vector<std::complex<float>>> &);

void ofdm_mpi::maximal_ratio_combining_mpi(std::vector<std::vector<std::complex<float>>> &, std::vector<std::complex<float>> &);

void ofdm_mpi::maximal_ratio_transmission_mpi(std::vector<std::complex<float>> &in, std::vector<std::vector<std::complex<float>>> &out) {

}

void ofdm_mpi::mod_one_frame_mpi(std::vector<std::vector<std::complex<float>>> &, std::vector<std::vector<std::complex<float>>> &);

void ofdm_mpi::demod_one_frame_mpi(std::vector<std::vector<std::complex<float>>> &, std::vector<std::vector<std::complex<float>>> &);

int ofdm_mpi::get_fft_size();
void ofdm_mpi::set_fft_size(int);
int ofdm_mpi::get_prefix_size();
void ofdm_mpi::set_prefix_size(int);