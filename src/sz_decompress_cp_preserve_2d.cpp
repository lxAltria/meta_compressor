#include "sz_decompress_3d.hpp"
#include "sz_decompress_cp_preserve_2d.hpp"
#include "sz_decompress_block_processing.hpp"
#include <limits>

template<typename T>
void
sz_decompress_cp_preserve_2d_offline(const unsigned char * compressed, size_t r1, size_t r2, T *& U, T *& V){
	if(U) free(U);
	if(V) free(V);
	size_t num_elements = r1 * r2;
	const unsigned char * compressed_pos = compressed;
	int base = 0;
	read_variable_from_src(compressed_pos, base);
	printf("base = %d\n", base);
	size_t compressed_eb_size = 0;
	read_variable_from_src(compressed_pos, compressed_eb_size);
	size_t compressed_u_size = 0;
	read_variable_from_src(compressed_pos, compressed_u_size);
	size_t compressed_v_size = 0;
	read_variable_from_src(compressed_pos, compressed_v_size);
	printf("eb_size = %ld, u_size = %ld, v_size = %ld\n", compressed_eb_size, compressed_u_size, compressed_v_size);
	int * type = Huffman_decode_tree_and_data(2*256, 2*num_elements, compressed_pos);
	double * eb = (double *) malloc(num_elements*sizeof(double));
	const double threshold=std::numeric_limits<float>::epsilon();
	for(int i=0; i<num_elements; i++){
		if(type[i] == 0) eb[i] = 0;
		else eb[i] = pow(base, type[i]) * threshold;
	}
	U = sz_decompress_2d_with_eb<T>(compressed_pos, eb, r1, r2);
	compressed_pos += compressed_u_size;
	for(int i=0; i<num_elements; i++){
		if(type[num_elements + i] == 0) eb[i] = 0;
		else eb[i] = pow(base, type[num_elements + i]) * threshold;
	}
	V = sz_decompress_2d_with_eb<T>(compressed_pos, eb, r1, r2);
	free(eb);
}

template
void
sz_decompress_cp_preserve_2d_offline<float>(const unsigned char * compressed, size_t r1, size_t r2, float *& U, float *& V);

template
void
sz_decompress_cp_preserve_2d_offline<double>(const unsigned char * compressed, size_t r1, size_t r2, double *& U, double *& V);

template<typename T>
void
sz_decompress_cp_preserve_2d_offline_log(const unsigned char * compressed, size_t r1, size_t r2, T *& U, T *& V){
	if(U) free(U);
	if(V) free(V);
	size_t num_elements = r1 * r2;
	const unsigned char * compressed_pos = compressed;
	int base = 0;
	read_variable_from_src(compressed_pos, base);
	printf("base = %d\n", base);
	size_t compressed_eb_size = 0;
	read_variable_from_src(compressed_pos, compressed_eb_size);
	size_t compressed_u_size = 0;
	read_variable_from_src(compressed_pos, compressed_u_size);
	size_t compressed_v_size = 0;
	read_variable_from_src(compressed_pos, compressed_v_size);
	printf("eb_size = %ld, u_size = %ld, v_size = %ld\n", compressed_eb_size, compressed_u_size, compressed_v_size);
	int * type = Huffman_decode_tree_and_data(2*256, num_elements, compressed_pos);
	double * eb = (double *) malloc(num_elements*sizeof(double));
	const double threshold=std::numeric_limits<float>::epsilon();
	for(int i=0; i<num_elements; i++){
		if(type[i] == 0) eb[i] = 0;
		else eb[i] = pow(base, type[i]) * threshold;
	}
	size_t sign_map_size = (num_elements - 1)/8 + 1;
	unsigned char * sign_map_u = convertByteArray2IntArray_fast_1b_sz(num_elements, compressed_pos, sign_map_size);	
	unsigned char * sign_map_v = convertByteArray2IntArray_fast_1b_sz(num_elements, compressed_pos, sign_map_size);	
	U = sz_decompress_2d_with_eb<T>(compressed_pos, eb, r1, r2);
	compressed_pos += compressed_u_size;
	for(int i=0; i<num_elements; i++){
		U[i] = sign_map_u[i] ? exp2(U[i]) : -exp2(U[i]);
	}
	V = sz_decompress_2d_with_eb<T>(compressed_pos, eb, r1, r2);
	for(int i=0; i<num_elements; i++){
		V[i] = sign_map_v[i] ? exp2(V[i]) : -exp2(V[i]);
	}
	free(sign_map_u);
	free(sign_map_v);
	free(eb);
}

template
void
sz_decompress_cp_preserve_2d_offline_log<float>(const unsigned char * compressed, size_t r1, size_t r2, float *& U, float *& V);

template
void
sz_decompress_cp_preserve_2d_offline_log<double>(const unsigned char * compressed, size_t r1, size_t r2, double *& U, double *& V);

