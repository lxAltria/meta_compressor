#include "sz_compress_pwr.hpp"
#include "sz_compress_3d.hpp"
#include "sz_compression_utils.hpp"
#include <cmath>
#include <cstring>
#include <limits>

template<typename T>
unsigned char *
sz_compress_2d_pwr(const T * data, size_t r1, size_t r2, double precision, size_t& compressed_size, int BSIZE, bool block_independant){
	size_t num_elements = r1 * r2;
	T * log_data = (T *) malloc(num_elements * sizeof(T));

	unsigned char * signs = (unsigned char *) malloc(num_elements);
	memset(signs, 0, num_elements);
	// preprocess
	T min_data = data[0];
	T max_data = data[0];
	for(int i=1; i<num_elements; i++){
		min_data = MIN(min_data, data[i]);
		max_data = MAX(max_data, data[i]);
	}
	T max_abs_log_data = 0;
    if(min_data == 0) max_abs_log_data = fabs(log2(fabs(max_data)));
    else if(max_data == 0) max_abs_log_data = fabs(log2(fabs(min_data)));
    else max_abs_log_data = fabs(log2(fabs(min_data))) > fabs(log2(fabs(max_data))) ? fabs(log2(fabs(min_data))) : fabs(log2(fabs(max_data)));
    T min_log_data = max_abs_log_data;
	for(size_t i=0; i<num_elements; i++){
		if(data[i] < 0){
			signs[i] = 1;
			log_data[i] = -data[i];
		}
		else
			log_data[i] = data[i];
		if(log_data[i] > 0){
			log_data[i] = log2(log_data[i]);
			if(log_data[i] > max_abs_log_data) max_abs_log_data = log_data[i];
			if(log_data[i] < min_log_data) min_log_data = log_data[i];
		}
	}
	if(fabs(min_log_data) > max_abs_log_data) max_abs_log_data = fabs(min_log_data);
	double log_precision = log2(1.0 + precision) - max_abs_log_data * std::numeric_limits<T>::epsilon();
	for(size_t i=0; i<num_elements; i++){
		if(data[i] == 0){
			log_data[i] = min_log_data - 2.0001*log_precision;
		}
	}
	size_t log_compressed_size = 0;
	unsigned char * log_compressed = sz_compress_2d(log_data, r1, r2, log_precision, log_compressed_size, BSIZE, block_independant);
    free(log_data);
	size_t sign_size = (num_elements-1)/8 + 1;
	T threshold = min_log_data - 1.0001*log_precision;
	compressed_size = sizeof(T) + sign_size + log_compressed_size;
    unsigned char * final_compressed = (unsigned char *) malloc(compressed_size);
    unsigned char * final_compressed_pos = final_compressed;
    memcpy(final_compressed_pos, &threshold, sizeof(T));
    final_compressed_pos += sizeof(T);
	convertIntArray2ByteArray_fast_1b_to_result_sz(signs, num_elements, final_compressed_pos);
	memcpy(final_compressed_pos, log_compressed, log_compressed_size);
	free(signs);
	free(log_compressed);
	return final_compressed;
}

template 
unsigned char * 
sz_compress_2d_pwr<float>(const float * data, size_t r1, size_t r2, double precision, size_t& compressed_size, int BSIZE, bool block_independant);

