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

template<typename T>
void
sz_decompress_cp_preserve_2d_online(const unsigned char * compressed, size_t r1, size_t r2, T *& U, T *& V){
	if(U) free(U);
	if(V) free(V);
	size_t num_elements = r1 * r2;
	const unsigned char * compressed_pos = compressed;
	int base = 0;
	read_variable_from_src(compressed_pos, base);
	printf("base = %d\n", base);
	int intv_radius = 0;
	read_variable_from_src(compressed_pos, intv_radius);
	const int capacity = (intv_radius << 1);
	size_t unpred_data_count = 0;
	read_variable_from_src(compressed_pos, unpred_data_count);
	const T * unpred_data_pos = (T *) compressed_pos;
	compressed_pos += unpred_data_count*sizeof(T);
	int * eb_quant_index = Huffman_decode_tree_and_data(2*256, 2*num_elements, compressed_pos);
	int * data_quant_index = Huffman_decode_tree_and_data(2*capacity, 2*num_elements, compressed_pos);
	U = (T *) malloc(num_elements*sizeof(T));
	V = (T *) malloc(num_elements*sizeof(T));
	T * U_pos = U;
	T * V_pos = V;
	int * eb_quant_index_pos = eb_quant_index;
	int * data_quant_index_pos = data_quant_index;
	const double threshold=std::numeric_limits<float>::epsilon();
	for(int i=0; i<r1; i++){
		for(int j=0; j<r2; j++){
			// get eb
			if(*eb_quant_index_pos == 0){
				*U_pos = *(unpred_data_pos ++);
				*V_pos = *(unpred_data_pos ++);
				eb_quant_index_pos += 2;
			}
			else{
				for(int k=0; k<2; k++){
					T * cur_data_pos = (k == 0) ? U_pos : V_pos;					
					double eb = pow(base, *eb_quant_index_pos ++) * threshold;
					T d0 = (i && j) ? cur_data_pos[-1 - r2] : 0;
					T d1 = (i) ? cur_data_pos[-r2] : 0;
					T d2 = (j) ? cur_data_pos[-1] : 0;
					T pred = d1 + d2 - d0;
					*cur_data_pos = pred + 2 * (data_quant_index_pos[k] - intv_radius) * eb;
				}
			}
			U_pos ++;
			V_pos ++;
			data_quant_index_pos += 2;
		}
	}
	free(eb_quant_index);
	free(data_quant_index);
}

template
void
sz_decompress_cp_preserve_2d_online<float>(const unsigned char * compressed, size_t r1, size_t r2, float *& U, float *& V);

template
void
sz_decompress_cp_preserve_2d_online<double>(const unsigned char * compressed, size_t r1, size_t r2, double *& U, double *& V);

