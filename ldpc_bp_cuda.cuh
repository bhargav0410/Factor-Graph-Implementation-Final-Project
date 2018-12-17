#ifndef LDPC_BP_CUDA_CUH
#define LDPC_BP_CUDA_CUH

#include "ldpc_bp.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cmath>
#include <algorithm>
#include <chrono>

using namespace std::chrono;

//Encodes the input vector in GPU using the generator matrix
__global__ void encode(int *G_mat, int *in, int *out, int len, int num_vecs, int num_msg_bits) {
    int tidx = threadIdx.x + blockDim.x*blockIdx.x;
    int tidy = threadIdx.y + blockDim.y*blockIdx.y;
    int tid = tidx + len*tidy;
    if (tidx < len) {
        if (tidx < num_msg_bits) {
            out[tid] = in[tidx + num_msg_bits*tidy];
        } else {
            int temp = 0;
            for (int i = 0; i < num_msg_bits; i++) {
                temp = (temp + (in[num_msg_bits*tidy + i]*G_mat[tidx + len*i])) % 2;
            }
            out[tid] = temp;
        }
    }
}

//Sum product i.e. belief propagation

__global__ void get_output(int *H_mat, float *extrin_llr, float *llr, int *out, int rows_, int cols_, int num_msg_bits) {
    int rows = rows_, cols = cols_;
    int tidx = threadIdx.x + blockDim.x*blockIdx.x, tidy = threadIdx.y + blockDim.y*blockIdx.y;
    int tid = tidx + cols*tidy;
    int tidz = rows*cols*blockIdx.z;
    float temp;
    if (tidx < cols && tidy < rows) {
        if (H_mat[tid] > 0) {
            //Variable node update step
            temp = llr[tidx + cols*blockIdx.z];
            for (int k = 0; k < rows; k++) {
                if (H_mat[tidx + k*cols] > 0) {
                    temp += extrin_llr[tidx + k*cols + tidz];
                }
            }
            //llr[tidx] = temp;
            //printf("(%f, %f, %d)\n", llr[tidx], temp, tidx);

            if (tidx < num_msg_bits) {
                if (temp >= 0) {
                    out[tidx + num_msg_bits*blockIdx.z] = 1;
                } else {
                    out[tidx + num_msg_bits*blockIdx.z] = 0;
                }
            }
        }
    }
}

//Variable nodes update their own partial LLR
__global__ void var_node_updation_step(int *H_mat, float *intrin_llr, float *extrin_llr, float *llr, int rows_, int cols_) {
    int rows = rows_, cols = cols_;
    int tidx = threadIdx.x + blockDim.x*blockIdx.x, tidy = threadIdx.y + blockDim.y*blockIdx.y;
    int tid = tidx + cols*tidy;
    int tidz = rows*cols*blockIdx.z;
    float temp;
    if (tidx < cols && tidy < rows) {
        if (H_mat[tid] > 0) {
            //Variable node update step
            temp = llr[tidx + cols*blockIdx.z];
            for (int k = 0; k < rows; k++) {
                if (H_mat[tidx + k*cols] > 0 && k != tidy) {
                    temp += extrin_llr[tidx + k*cols + tidz];
                }
            }
            intrin_llr[tid + tidz] = temp;
        }
    }
}

//Check nodes update their partial LLR
__global__ void check_node_updation_step(int *H_mat, float *intrin_llr, float *extrin_llr, int rows_, int cols_) {
    int rows = rows_, cols = cols_;
    int tidx = threadIdx.x + blockDim.x*blockIdx.x, tidy = threadIdx.y + blockDim.y*blockIdx.y;
    int tid = tidx + cols*tidy;
    int tidz = rows*cols*blockIdx.z;
    float temp;
    if (tidx < cols && tidy < rows) {
        if (H_mat[tid] > 0) {
            temp = 0;
            //Check node update (horizontal step)
            for (int k = 0; k < cols; k++) {
                if (H_mat[cols*tidy + k] > 0 && k != tidx) {
                    temp *= tanhf(intrin_llr[tidy*cols + k + tidz]/(float)2.0);
                }
            }
            temp = logf((1 + temp)/(1 - temp));
            extrin_llr[tid + tidz] = temp;
        }
    }
}

//Initializing belief propagation
__global__ void belief_propagation_init(int *H_mat, float *intrin_llr, float *llr, float *in, int rows_, int cols_, float snr_) {
    int rows = rows_, cols = cols_;
    float snr = snr_;
    int tidx = threadIdx.x + blockDim.x*blockIdx.x, tidy = threadIdx.y + blockDim.y*blockIdx.y;
    int tid = tidx + cols*tidy;
    int tidz = rows*cols*blockIdx.z;
    float temp;
    if (tidx < cols && tidy < rows) {
        //Getting LLR values from input vector
        temp = 2 * in[tidx + cols*blockIdx.z] * powf((float)10, (float)snr/(float)10);
        /*
        if (tidx == 0)
            printf("(%f, %f)\n", temp, snr);
            */
        if (H_mat[tid] > 0) {
            intrin_llr[tid + tidz] = temp;
        }
        if (tidy == 0) {
            llr[tidx + cols*blockIdx.z] = temp;
        }
    }
}

class ldpc_bp_cuda : public ldpc_bp {

public:

    double encode_using_G_mat_cuda(std::vector<int> &in, std::vector<int> &out) {
        duration<double> timediff;
        high_resolution_clock::time_point start, finish;
        cudaDeviceProp devProp;
        cudaGetDeviceProperties(&devProp, 0);
        int threads_per_block  = devProp.maxThreadsPerBlock;
        while (standard_form_var == 0) {
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

        //Resizing the output vector
        if (out.size() != ceil(n*((float)len/(float)num_msg_bits))) {
            out.resize(ceil(n*((float)len/(float)num_msg_bits)));
        }

        //Device vector allocation
        int *dev_G, *dev_in, *dev_out;
        cudaMalloc((void **)&dev_G, (int)G_mat.size()*(int)G_mat[0].size()*sizeof(*dev_G));
        cudaMalloc((void **)&dev_in, (int)in.size()*sizeof(*dev_in));
        cudaMalloc((void **)&dev_out, (int)out.size()*sizeof(*dev_out));
        //Copying from host to device
        for (int i = 0; i < G_mat.size(); i++) {
            cudaMemcpy(&dev_G[i*(int)G_mat[0].size()], &G_mat[i][0], (int)G_mat[0].size()*sizeof(*dev_G), cudaMemcpyHostToDevice);
        }

        //Number of threads and blocks
        int thread_x = std::min((int)floor(sqrt(threads_per_block)), (int)n);
        int thread_y = std::min((int)floor(sqrt(threads_per_block)), (int)(int)ceil((float)len/(float)num_msg_bits));//std::min((int)floor(sqrt(threads_per_block)), (int)G_mat.size());
        int block_x = (int)ceil((float)n/(float)thread_x);
        int block_y = (int)ceil((float)ceil((float)len/(float)num_msg_bits)/(float)thread_y);
        int block_z = 1;//(int)ceil((float)len/(float)num_msg_bits);
        dim3 gridDims(block_x, block_y, block_z), blockDims(thread_x, thread_y);

        start = high_resolution_clock::now();
        cudaMemcpy(dev_in, &in[0], (int)in.size()*sizeof(*dev_in), cudaMemcpyHostToDevice);
        encode << <gridDims, blockDims>> > (dev_G, dev_in, dev_out, n, (int)ceil((float)len/(float)num_msg_bits), num_msg_bits);
        cudaMemcpy(&out[0], dev_out, (int)out.size()*sizeof(*dev_out), cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        finish = high_resolution_clock::now();

        cudaFree(dev_G);
        cudaFree(dev_in);
        cudaFree(dev_out);

        return duration_cast<duration<double>>(finish - start).count();
    }

    double sum_product_decoding_cuda(std::vector<float> &in_vec, std::vector<int> &out_vec, float snr, int iter) {
        duration<double> timediff;
        high_resolution_clock::time_point start, finish;
        cudaDeviceProp devProp;
        cudaGetDeviceProperties(&devProp, 0);
        int threads_per_block  = devProp.maxThreadsPerBlock;
        if (G_mat.size() == 0) {
            std::cout << "Need actual rate for decoding...Create generator matrix\n";
            return -1;
        }
        //Making sure input and output vectors have some values and are of correct sizes
        int in_vec_size = in_vec.size();
        if (in_vec_size == 0) {
            std::cout << "Input vector does not have input...\n";
            return -1;
        } else if (in_vec_size % n != 0) {
            std::cout << "Encoded vector size incorrect...\n";
            return -1;
        }
        if (out_vec.size() != G_mat.size()*(int)((float)in_vec_size/(float)n)) {
            out_vec.resize(G_mat.size()*(int)((float)in_vec_size/(float)n));
        }
        int num_msg_bits = G_mat.size();

        //Device vector allocation
        int *dev_H, *dev_out;
        float *extrin_llr, *intrin_llr, *llr, *dev_in; 
        cudaMalloc((void **)&dev_H, (int)H_mat.size()*(int)H_mat[0].size()*sizeof(*dev_H));
        cudaMalloc((void **)&extrin_llr, in_vec_size/n*(int)H_mat.size()*(int)H_mat[0].size()*sizeof(*extrin_llr));
        cudaMalloc((void **)&intrin_llr, in_vec_size/n*(int)H_mat.size()*(int)H_mat[0].size()*sizeof(*intrin_llr));
        cudaMalloc((void **)&llr, in_vec_size/n*(int)H_mat[0].size()*sizeof(*llr));
        cudaMalloc((void **)&dev_in, (int)in_vec.size()*sizeof(*dev_in));
        cudaMalloc((void **)&dev_out, (int)out_vec.size()*sizeof(*dev_out));
        //Copying from host to device
        for (int i = 0; i < H_mat.size(); i++) {
            cudaMemcpy(&dev_H[i*(int)G_mat[0].size()], &H_mat[i][0], (int)H_mat[0].size()*sizeof(*dev_H), cudaMemcpyHostToDevice);
        }

        //Number of threads and blocks
        int thread_x = std::min((int)floor(sqrt(threads_per_block)), (int)H_mat[0].size());
        int thread_y = std::min((int)floor(sqrt(threads_per_block)), (int)H_mat.size());//std::min((int)floor(sqrt(threads_per_block)), (int)G_mat.size());
        int block_x = (int)ceil((float)H_mat[0].size()/(float)thread_x);
        int block_y = (int)ceil((float)H_mat.size()/(float)thread_y);
        int block_z = in_vec_size/n;//(int)ceil((float)len/(float)num_msg_bits);
        dim3 gridDims(block_x, block_y, block_z), blockDims(thread_x, thread_y);
        //Sum product decoding
        start = high_resolution_clock::now();
        cudaMemcpy(dev_in, &in_vec[0], in_vec.size()*sizeof(*dev_in), cudaMemcpyHostToDevice);
        //for (int i = 0; i < in_vec_size/n; i++) {
            belief_propagation_init <<<gridDims, blockDims>>> (dev_H, intrin_llr, llr, dev_in, (int)H_mat.size(), n, snr);
           //cudaDeviceSynchronize();
            for (int it = 0; it < iter; it++) {
                check_node_updation_step <<<gridDims, blockDims>>> (dev_H, intrin_llr, extrin_llr, (int)H_mat.size(), n);
                if (it < iter - 1)
                    var_node_updation_step <<<gridDims, blockDims>>> (dev_H, intrin_llr, extrin_llr, llr, (int)H_mat.size(), n);
            }
            get_output <<<gridDims, blockDims>>> (dev_H, extrin_llr, llr, dev_out, (int)H_mat.size(), n, num_msg_bits);
        //}
        cudaMemcpy(&out_vec[0], dev_out, out_vec.size()*sizeof(*dev_out), cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        finish = high_resolution_clock::now();
        cudaFree(dev_in);
        cudaFree(dev_out);
        cudaFree(dev_H);
        cudaFree(extrin_llr);
        cudaFree(intrin_llr);
        cudaFree(llr);
        cudaDeviceReset();
        return duration_cast<duration<double>>(finish - start).count();
    }
};

#endif