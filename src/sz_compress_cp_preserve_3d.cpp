#include "sz_cp_preserve_utils.hpp"
#include "sz_compress_3d.hpp"
#include "sz_compress_cp_preserve_3d.hpp"
#include "sz_def.hpp"
#include "sz_compression_utils.hpp"
#include <unordered_map>

// maximal error bound to keep the sign of A*(1 + e_1) + B*(1 + e_2) + C*(1+e_3) + D
template<typename T>
inline double max_eb_to_keep_sign_3d_online(const T A, const T B, const T C, const T D=0){
	// double positive = 0, negative = 0;
 // 	accumulate(A, positive, negative);
	// accumulate(B, positive, negative);
	// accumulate(C, positive, negative);
	// double sum = A + B + C + D;
	// return fabs(sum) / (positive + negative);
	return fabs(A + B + C + D) / (fabs(A) + fabs(B) + fabs(C));
}

template<typename T>
double 
max_eb_to_keep_position_and_type_3d_online(const T u0, const T u1, const T u2, const T u3, const T v0, const T v1, const T v2, const T v3,
	const T w0, const T w1, const T w2, const T w3, int n, std::unordered_map<int, vector<bool>>& flags){
	//det = -u2 v1 w0 + u3 v1 w0 + u1 v2 w0 - u3 v2 w0 - u1 v3 w0 + u2 v3 w0 + 
	//  u2 v0 w1 - u3 v0 w1 - u0 v2 w1 + u3 v2 w1 + u0 v3 w1 - u2 v3 w1 - 
	//  u1 v0 w2 + u3 v0 w2 + u0 v1 w2 - u3 v1 w2 - u0 v3 w2 + u1 v3 w2 + 
	//  u1 v0 w3 - u2 v0 w3 - u0 v1 w3 + u2 v1 w3 + u0 v2 w3 - u1 v2 w3
	//    = P0 + P1 + P2 + P3
	// mu0 = (u3 v2 w1 - u2 v3 w1 - u3 v1 w2 + u1 v3 w2 + u2 v1 w3 - u1 v2 w3) / det = P0/(P1 + P2 + P3 + P4)
	// mu1 = (-u3 v2 w0 + u2 v3 w0 + u3 v0 w2 - u0 v3 w2 - u2 v0 w3 + u0 v2 w3) / det = P1/(P1 + P2 + P3 + P4)
	// mu2 = (u3 v1 w0 - u1 v3 w0 - u3 v0 w1 + u0 v3 w1 + u1 v0 w3 - u0 v1 w3) / det = P2/(P1 + P2 + P3 + P4)
	// mu3 = (-u2 v1 w0 + u1 v2 w0 + u2 v0 w1 - u0 v2 w1 - u1 v0 w2 + u0 v1 w2) / det = P3/(P1 + P2 + P3 + P4)
	std::vector<bool> flag;
	if(flags.find(n) == flags.end()){
	    double M0 = -u1*v2*w3 + u1*v3*w2 - u2*v3*w1 + u2*v1*w3 - u3*v1*w2 + u3*v2*w1;
	    double M1 = -u0*v3*w2 + u0*v2*w3 - u2*v0*w3 + u2*v3*w0 - u3*v2*w0 + u3*v0*w2;
	    double M2 = -u0*v1*w3 + u0*v3*w1 - u1*v3*w0 + u1*v0*w3 - u3*v0*w1 + u3*v1*w0;
	    double M3 = u0*v1*w2 - u0*v2*w1 + u1*v2*w0 - u1*v0*w2 + u2*v0*w1 - u2*v1*w0;
	    double M = M0 + M1 + M2 + M3;
	    flag.push_back(M / M0 >= 1);
	    flag.push_back(M / M1 >= 1);
	    flag.push_back(M / M2 >= 1);
	    flag.push_back(M / M3 >= 1);
	    // flags.insert(make_pair(n, flag));
	}
	else{
		flag = flags[n];
	}
	double u3_0 = u3*v1*w2 + u3*v2*w1, u3_1 = - u3*v2*w0 + u3*v0*w2, u3_2 = - u3*v0*w1 + u3*v1*w0;
	double v3_0 = u1*v3*w2 - u2*v3*w1, v3_1 = u2*v3*w0 - u0*v3*w2, v3_2 = u0*v3*w1 - u1*v3*w0;
	double w3_0 = - u1*v2*w3 + u2*v1*w3, w3_1 = u0*v2*w3 - u2*v0*w3, w3_2 = - u0*v1*w3 + u1*v0*w3;
	double c_4 = u0*v1*w2 - u0*v2*w1 + u1*v2*w0 - u1*v0*w2 + u2*v0*w1 - u2*v1*w0;
	if(flag[0] && flag[1] && flag[2] && flag[3]){
		// cp found
		double eb = 1;
		double cur_eb = 0;
		cur_eb = MIN(max_eb_to_keep_sign_3d_online(u3_0, v3_0, w3_0), 
				max_eb_to_keep_sign_3d_online(u3_1 + u3_2, v3_1 + v3_2, w3_1 + w3_2, c_4));
		eb = MIN(eb, cur_eb);
		cur_eb = MIN(max_eb_to_keep_sign_3d_online(u3_1, v3_1, w3_1), 
				max_eb_to_keep_sign_3d_online(u3_0 + u3_2, v3_0 + v3_2, w3_0 + w3_2, c_4));
		eb = MIN(eb, cur_eb);
		cur_eb = MIN(max_eb_to_keep_sign_3d_online(u3_2, v3_2, w3_2), 
				max_eb_to_keep_sign_3d_online(u3_0 + u3_1, v3_0 + v3_1, w3_0 + w3_1, c_4));
		eb = MIN(eb, cur_eb);
		cur_eb = max_eb_to_keep_sign_3d_online(u3_0 + u3_1 + u3_2, v3_0 + v3_1 + v3_2, w3_0 + w3_1 + w3_2);
		eb = MIN(eb, cur_eb);
		return eb;
	}
	else{
		double eb = 0;
		double cur_eb = 0;
		if(!flag[0]){
			cur_eb = MIN(max_eb_to_keep_sign_3d_online(u3_0, v3_0, w3_0), 
					max_eb_to_keep_sign_3d_online(u3_1 + u3_2, v3_1 + v3_2, w3_1 + w3_2, c_4));
			eb = MAX(eb, cur_eb);
		}
		if(!flag[1]){
			cur_eb = MIN(max_eb_to_keep_sign_3d_online(u3_1, v3_1, w3_1), 
					max_eb_to_keep_sign_3d_online(u3_0 + u3_2, v3_0 + v3_2, w3_0 + w3_2, c_4));
			eb = MAX(eb, cur_eb);
		}
		if(!flag[2]){
			cur_eb = MIN(max_eb_to_keep_sign_3d_online(u3_2, v3_2, w3_2), 
					max_eb_to_keep_sign_3d_online(u3_0 + u3_1, v3_0 + v3_1, w3_0 + w3_1, c_4));
			eb = MAX(eb, cur_eb);
		}
		if(!flag[3]){
			cur_eb = max_eb_to_keep_sign_3d_online(u3_0 + u3_1 + u3_2, v3_0 + v3_1 + v3_2, w3_0 + w3_1 + w3_2);
			eb = MAX(eb, cur_eb);
		}
		return eb;
	}
}

template<typename T>
unsigned char *
sz_compress_cp_preserve_3d_online_log(const T * U, const T * V, const T * W, size_t r1, size_t r2, size_t r3, size_t& compressed_size, bool transpose, double max_pwr_eb){

	size_t num_elements = r1 * r2 * r3;
	size_t sign_map_size = (num_elements - 1)/8 + 1;
	unsigned char * sign_map_compressed = (unsigned char *) malloc(3*sign_map_size);
	unsigned char * sign_map_compressed_pos = sign_map_compressed;
	unsigned char * sign_map = (unsigned char *) malloc(num_elements*sizeof(unsigned char));
	// Note the convert function has address auto increment
	T * log_U = log_transform(U, sign_map, num_elements);
	convertIntArray2ByteArray_fast_1b_to_result_sz(sign_map, num_elements, sign_map_compressed_pos);
	T * log_V = log_transform(V, sign_map, num_elements);
	convertIntArray2ByteArray_fast_1b_to_result_sz(sign_map, num_elements, sign_map_compressed_pos);
	T * log_W = log_transform(W, sign_map, num_elements);
	convertIntArray2ByteArray_fast_1b_to_result_sz(sign_map, num_elements, sign_map_compressed_pos);
	free(sign_map);

	T * decompressed_U = (T *) malloc(num_elements*sizeof(T));
	memcpy(decompressed_U, U, num_elements*sizeof(T));
	T * decompressed_V = (T *) malloc(num_elements*sizeof(T));
	memcpy(decompressed_V, V, num_elements*sizeof(T));
	T * decompressed_W = (T *) malloc(num_elements*sizeof(T));
	memcpy(decompressed_W, W, num_elements*sizeof(T));

	int * eb_quant_index = (int *) malloc(num_elements*sizeof(int));
	int * data_quant_index = (int *) malloc(3*num_elements*sizeof(int));
	int * eb_quant_index_pos = eb_quant_index;
	int * data_quant_index_pos = data_quant_index;
	// next, row by row
	const int base = 2;
	const double log_of_base = log2(base);
	const int capacity = 65536;
	const int intv_radius = (capacity >> 1);
	// offsets to get 24 adjacent simplex indices
	// x -> z, high -> low
	// current data would always be the last index, i.e. x[i][3]
	const int coordinates[24][4][3] = {
		// offset = 0, 0, 0
		{
			{0, 0, 1},
			{0, 1, 1},
			{1, 1, 1},
			{0, 0, 0}
		},
		{
			{0, 1, 0},
			{0, 1, 1},
			{1, 1, 1},
			{0, 0, 0}
		},
		{
			{0, 0, 1},
			{1, 0, 1},
			{1, 1, 1},
			{0, 0, 0}
		},
		{
			{1, 0, 0},
			{1, 0, 1},
			{1, 1, 1},
			{0, 0, 0}
		},
		{
			{0, 1, 0},
			{1, 1, 0},
			{1, 1, 1},
			{0, 0, 0}
		},
		{
			{1, 0, 0},
			{1, 1, 0},
			{1, 1, 1},
			{0, 0, 0}
		},
		// offset = -1, 0, 0
		{
			{0, 0, 0},
			{1, 0, 1},
			{1, 1, 1},
			{1, 0, 0}
		},
		{
			{0, 0, 0},
			{1, 1, 0},
			{1, 1, 1},
			{1, 0, 0}
		},
		// offset = 0, -1, 0
		{
			{0, 0, 0},
			{0, 1, 1},
			{1, 1, 1},
			{0, 1, 0}
		},
		{
			{0, 0, 0},
			{1, 1, 0},
			{1, 1, 1},
			{1, 1, 0}
		},
		// offset = -1, -1, 0
		{
			{0, 0, 0},
			{0, 1, 0},
			{1, 1, 1},
			{1, 1, 0}
		},
		{
			{0, 0, 0},
			{1, 0, 1},
			{1, 1, 1},
			{0, 1, 0}
		},
		// offset = 0, 0, -1
		{
			{0, 0, 0},
			{0, 1, 1},
			{1, 1, 0},
			{0, 0, 1}
		},
		{
			{0, 0, 0},
			{1, 0, 1},
			{1, 1, 1},
			{0, 0, 1}
		},
		// offset = -1, 0, -1
		{
			{0, 0, 0},
			{0, 0, 1},
			{1, 1, 1},
			{1, 0, 1}
		},
		{
			{0, 0, 0},
			{1, 0, 0},
			{1, 1, 1},
			{1, 0, 1}
		},
		// offset = 0, -1, -1
		{
			{0, 0, 0},
			{0, 0, 1},
			{1, 1, 1},
			{0, 1, 1}
		},
		{
			{0, 0, 0},
			{0, 1, 0},
			{1, 1, 1},
			{0, 1, 1}
		},
		// offset = -1, -1, -1
		{
			{0, 0, 0},
			{0, 0, 1},
			{0, 1, 1},
			{1, 1, 1}
		},
		{
			{0, 0, 0},
			{0, 1, 0},
			{0, 1, 1},
			{1, 1, 1}
		},
		{
			{0, 0, 0},
			{0, 0, 1},
			{1, 0, 1},
			{1, 1, 1}
		},
		{
			{0, 0, 0},
			{1, 0, 0},
			{1, 0, 1},
			{1, 1, 1}
		},
		{
			{0, 0, 0},
			{0, 1, 0},
			{1, 1, 0},
			{1, 1, 1}
		},
		{
			{0, 0, 0},
			{1, 0, 0},
			{1, 1, 0},
			{1, 1, 1}
		}
	};
	T * cur_log_U_pos = log_U;
	T * cur_log_V_pos = log_V;
	T * cur_log_W_pos = log_W;
	T * cur_U_pos = decompressed_U;
	T * cur_V_pos = decompressed_V;
	T * cur_W_pos = decompressed_W;
	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;
	ptrdiff_t offset[24][3];
	for(int i=0; i<24; i++){
		for(int x=0; x<3; x++){
			offset[i][x] = (coordinates[i][x][0] - coordinates[i][3][0]) * dim0_offset + (coordinates[i][x][1] - coordinates[i][3][1]) * dim1_offset + (coordinates[i][x][2] - coordinates[i][3][2]);
		}
	}
	unpred_vec<T> unpred_data = unpred_vec<T>();
	ptrdiff_t max_pointer_pos = num_elements;
	std::unordered_map<int, vector<bool>> flags;
	for(int i=0; i<r1; i++){
		// printf("start %d row\n", i);
		for(int j=0; j<r2; j++){
			for(int k=0; k<r3; k++){
				double required_eb = max_pwr_eb;
				// derive eb given 24 adjacent simplex
				for(int n=0; n<24; n++){
					bool in_mesh = true;
					for(int p=0; p<3; p++){
						if(!pointer_in_range(cur_U_pos - decompressed_U + offset[n][p], max_pointer_pos)){
							in_mesh = false;
							break;
						}
					}
					if(in_mesh){
						int index = 0; // TODO: define index for each simplex
						required_eb = MIN(required_eb, max_eb_to_keep_position_and_type_3d_online(
							cur_U_pos[offset[n][0]], cur_U_pos[offset[n][1]], cur_U_pos[offset[n][2]], *cur_U_pos,
							cur_V_pos[offset[n][0]], cur_V_pos[offset[n][1]], cur_V_pos[offset[n][2]], *cur_V_pos,
							cur_W_pos[offset[n][0]], cur_W_pos[offset[n][1]], cur_W_pos[offset[n][2]], *cur_W_pos,
							index, flags));
					}
				}
				if(required_eb > 0){
					bool unpred_flag = false;
					T decompressed[3];
					double abs_eb = log(1 + required_eb);
					*eb_quant_index_pos = eb_exponential_quantize(abs_eb, base, log_of_base);
					if(*eb_quant_index_pos > 0){
						// compress vector fields
						T * log_data_pos[3] = {cur_log_U_pos, cur_log_V_pos, cur_log_W_pos};
						T * data_pos[3] = {cur_U_pos, cur_V_pos, cur_W_pos};
						for(int p=0; p<3; p++){
							T * cur_log_data_pos = log_data_pos[p];
							T cur_data = *cur_log_data_pos;
							// get adjacent data and perform Lorenzo
							/*
								d6	X
								d4	d5
								d2	d3
								d0	d1
							*/
							T d0 = (i && j && k) ? cur_log_data_pos[- dim0_offset - dim1_offset - 1] : 0;
							T d1 = (i && j) ? cur_log_data_pos[- dim0_offset - dim1_offset] : 0;
							T d2 = (i && k) ? cur_log_data_pos[- dim0_offset - 1] : 0;
							T d3 = (i) ? cur_log_data_pos[- dim0_offset] : 0;
							T d4 = (j && k) ? cur_log_data_pos[- dim1_offset - 1] : 0;
							T d5 = (j) ? cur_log_data_pos[- dim1_offset] : 0;
							T d6 = (k) ? cur_log_data_pos[- 1] : 0;
							T pred = d0 + d3 + d5 + d6 - d1 - d2 - d4;
							double diff = cur_data - pred;
							double quant_diff = fabs(diff) / abs_eb + 1;
							if(quant_diff < capacity){
								quant_diff = (diff > 0) ? quant_diff : -quant_diff;
								int quant_index = (int)(quant_diff/2) + intv_radius;
								data_quant_index_pos[p] = quant_index;
								decompressed[p] = pred + 2 * (quant_index - intv_radius) * abs_eb; 
								// check original data
								if(fabs(decompressed[p] - cur_data) >= abs_eb){
									unpred_flag = true;
									break;
								}
							}
							else{
								unpred_flag = true;
								break;
							}
						}
					}
					else unpred_flag = true;
					if(unpred_flag){
						*(eb_quant_index_pos ++) = 0;
						*(data_quant_index_pos ++) = intv_radius;
						*(data_quant_index_pos ++) = intv_radius;
						*(data_quant_index_pos ++) = intv_radius;
						unpred_data.push_back(*cur_U_pos);
						unpred_data.push_back(*cur_V_pos);
						unpred_data.push_back(*cur_W_pos);
					}
					else{
						eb_quant_index_pos ++;
						data_quant_index_pos += 3;
						*cur_log_U_pos = decompressed[0];
						*cur_log_V_pos = decompressed[1];
						*cur_log_W_pos = decompressed[2];
						*cur_U_pos = (*cur_U_pos > 0) ? exp2(*cur_log_U_pos) : -exp2(*cur_log_U_pos);
						*cur_V_pos = (*cur_V_pos > 0) ? exp2(*cur_log_V_pos) : -exp2(*cur_log_V_pos);
						*cur_W_pos = (*cur_W_pos > 0) ? exp2(*cur_log_W_pos) : -exp2(*cur_log_W_pos);
					}
				}
				else{
					// record as unpredictable data
					*(eb_quant_index_pos ++) = 0;
					*(data_quant_index_pos ++) = intv_radius;
					*(data_quant_index_pos ++) = intv_radius;
					*(data_quant_index_pos ++) = intv_radius;
					unpred_data.push_back(*cur_U_pos);
					unpred_data.push_back(*cur_V_pos);
					unpred_data.push_back(*cur_W_pos);
				}
				cur_log_U_pos ++, cur_log_V_pos ++, cur_log_W_pos ++;
				cur_U_pos ++, cur_V_pos ++, cur_W_pos ++;
			}
		}
	}
	free(log_U);
	free(log_V);
	free(log_W);
	free(decompressed_U);
	free(decompressed_V);
	free(decompressed_W);
	printf("offset eb_q, data_q, unpred: %ld %ld %ld\n", eb_quant_index_pos - eb_quant_index, data_quant_index_pos - data_quant_index, unpred_data.size());
	unsigned char * compressed = (unsigned char *) malloc(2*num_elements*sizeof(T));
	unsigned char * compressed_pos = compressed;
	write_variable_to_dst(compressed_pos, base);
	write_variable_to_dst(compressed_pos, intv_radius);
	write_array_to_dst(compressed_pos, sign_map_compressed, 3*sign_map_size);
	free(sign_map_compressed);
	printf("offset after sign: %ld\n", compressed_pos - compressed);
	size_t unpredictable_count = unpred_data.size();
	write_variable_to_dst(compressed_pos, unpredictable_count);
	write_array_to_dst(compressed_pos, (T *)&unpred_data[0], unpredictable_count);	
	printf("offset after unpred: %ld\n", compressed_pos - compressed);
	Huffman_encode_tree_and_data(2*256, eb_quant_index, num_elements, compressed_pos);
	free(eb_quant_index);
	Huffman_encode_tree_and_data(2*capacity, data_quant_index, 3*num_elements, compressed_pos);
	free(data_quant_index);
	printf("offset after huffman: %ld\n", compressed_pos - compressed);
	compressed_size = compressed_pos - compressed;
	return compressed;	
}

template
unsigned char *
sz_compress_cp_preserve_3d_online_log(const float * U, const float * V, const float * W, size_t r1, size_t r2, size_t r3, size_t& compressed_size, bool transpose, double max_pwr_eb);

template
unsigned char *
sz_compress_cp_preserve_3d_online_log(const double * U, const double * V, const double * W, size_t r1, size_t r2, size_t r3, size_t& compressed_size, bool transpose, double max_pwr_eb);
