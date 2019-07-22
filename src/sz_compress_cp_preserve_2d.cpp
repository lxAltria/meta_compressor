#include "sz_compress_3d.hpp"
#include "sz_compress_cp_preserve_2d.hpp"
#include "sz_def.hpp"
#include "sz_compression_utils.hpp"
#include <vector>
#include <limits>
#include <cmath>
#include <cstdio>
#include <cstdlib>

template<typename T>
using unpred_vec = std::vector<T>;

inline int eb_exponential_quantize(double& eb, const int base, const double log_of_base, const double threshold=std::numeric_limits<float>::epsilon()){
	if(eb <= threshold){
		eb = 0;
		return 0;
	}
	int id = log2(eb / threshold)/log_of_base;
	eb = pow(base, id) * threshold;
	return id;
}

inline int eb_linear_quantize(double& eb, double threshold=1e-5){
	int id = eb / threshold;
	eb = id * threshold;
	return id;
}

template<typename T>
inline void swap_buffer_pointer(T*& b1, T*& b2){
	T * tmp = b1;
	b1 = b2;
	b2 = tmp;
}

template<typename T>
inline void accumulate(const T value, double& positive, double& negative){
	if(value >= 0) positive += value;
	else negative += - value;
}

// // maximal error bound to keep the sign of W_i 
// // e.g W0 = u1v2 - u2v1
// // where u2, v2 is current data, u1, v1 is decompressed data
// template<typename T>
// inline double max_eb_to_keep_sign_online_W(const T u1v2, const T u2v1){
// 	if(u1v2 * u2v1 <= 0) return 1;
// 	double c = u1v2 / u2v1;
// 	return fabs(c - 1) / (c + 1);
// }

// // maximal error bound to keep the sign of W_i(W_(i+1)+W_(i+2))
// // which preserves the FP/FN for positions
// // W0 + W1 = u1v2 - u2v1 + u2v0 - u0v2
// // W1 + W2 = u2v0 - u0v2 + u0v1 - u1v0
// // W2 + W0 = u0v1 - u1v0 + u1v2 - u2v1
// // where u2, v2 is current data, u0, v0, u1, v1 is decompressed data

// // W0 + W2
// template<typename T>
// inline double max_eb_to_keep_sign_online_W1_W2(const T u0v1, const T u1v0, const T u2v0, const T u0v2){
// 	T c = (u1v0 - u0v1) / (u2v0 - u0v2);
// 	if(c < 0) return 1;
// 	return (u2v0 - u0v2 > 0) ? 1 - c : c - 1;
// }

// // W1 + W2
// template<typename T>
// inline double max_eb_to_keep_sign_online_W0_W2(const T u0v1, const T u1v0, const T u1v2, const T u2v1){
// 	T c = (u1v0 - u0v1) / (u1v2 - u2v1);
// 	if(c < 0) return 1;
// 	return (u1v2 - u2v1 > 0) ? 1 - c : c - 1;
// }

// // W0 + W1
// template<typename T>
// inline double max_eb_to_keep_sign_online_W0_W1(const T u1v2, const T u2v1, const T u2v0, const T u0v2){
// 	double positive = 0;
// 	double negative = 0;
// 	accumulate(u1v2, positive, negative);
// 	accumulate(-u2v1, positive, negative);
// 	accumulate(u2v0, positive, negative);
// 	accumulate(-u0v2, positive, negative);
// 	if(positive * negative <= 0) return 1;
// 	double c = positive / negative;
// 	return MIN(1, fabs(c - 1) / (c + 1));
// }

// template<typename T>
// inline double max_eb_to_keep_position_online(const T u0v1, const T u1v0, const T u1v2, const T u2v1, const T u2v0, const T u0v2){
// 	// for W_i
// 	double eb = MIN(max_eb_to_keep_sign_online_W(u1v2, u2v1), max_eb_to_keep_sign_online_W(u2v0, u0v2));
// 	// for W1 + W2, W2 + W0
// 	eb = MIN(eb, MIN(max_eb_to_keep_sign_online_W1_W2(u0v1, u1v0, u2v0, u0v2), max_eb_to_keep_sign_online_W0_W2(u0v1, u1v0, u1v2, u2v1)));
// 	// for W0 + W1
// 	eb = MIN(eb, max_eb_to_keep_sign_online_W0_W1(u1v2, u2v1, u2v0, u0v2));
// 	return eb;
// }

// int x,y;
// /*
// triangle mesh Dx, Data0, Data1, derive cp-preserving eb for Dx given Data0, Data1
// */
// template<typename T>
// double 
// derive_cp_eb_for_positions_online(const T * Data0, const T * Data1, const T * DecData0, const T * DecData1, const T * Dx, double max_pwr_eb){
// 	T u0 = Data0[0], u1 = Data1[0], u2 = Dx[0];
// 	T v0 = Data0[1], v1 = Data1[1], v2 = Dx[1];
// 	double u1v2 = u1*v2;
// 	double u2v1 = u2*v1;
// 	double u2v0 = u2*v0;
// 	double u0v2 = u0*v2;
// 	double u0v1 = u0*v1;
// 	double u1v0 = u1*v0;
// 	T du0 = DecData0[0], du1 = DecData1[0], u2 = Dx[0];
// 	T dv0 = DecData0[1], dv1 = DecData1[1], v2 = Dx[1];
// 	double du1v2 = du1*v2;
// 	double u2dv1 = u2*dv1;
// 	double u2dv0 = u2*dv0;
// 	double du0v2 = du0*v2;
// 	double du0v1 = du0*dv1;
// 	double du1dv0 = du1*dv0;
// 	double eb = 0;
//     double det = u0v1 - u1v0 + u1v2 - u2v1 + u2v0 - u0v2;
//     double d_det = du0dv1 - du1dv0 + du1v2 - u2dv1 + u2dv0 - du0v2;
//     if(det * d_det <= 0) return -1;
//     if(det != 0){
// 		bool f1 = (det / (u2v0 - u0v2) >= 1);
// 		bool f2 = (det / (u1v2 - u2v1) >= 1); 
// 		bool f3 = (det / (u0v1 - u1v0) >= 1); 
// 		if(f1 && f2 && f3){
// 			// eb = max_eb_to_keep_position_online(u0v1, u1v0, u1v2, u2v1, u2v0, u0v2);
// 			eb = max_eb_to_keep_position_online(du0dv1, du1dv0, du1v2, u2dv1, u2dv0, du0v2);
// 			if((x>=6 && x<506) && (y>=6 && y<506) && (eb == 0)) {
// 				printf("%d %d, eb = 0\n", x, y);
// 			}
// 		}
// 		else{
// 			eb = 0;
// 			if(!f1){
// 				// W1(W0 + W2)
// 				double cur_eb = MIN(max_eb_to_keep_sign_online_W(u2v0, u0v2), max_eb_to_keep_sign_online_W0_W2(u0v1, u1v0, u1v2, u2v1));
// 				// double cur_eb = MIN(max_eb_to_keep_sign_online_W(u2dv0, du0v2), max_eb_to_keep_sign_online_W0_W2(du0dv1, du1dv0, du1v2, u2dv1));
// 				eb = MAX(eb, cur_eb);
// 			}
// 			if(!f2){
// 				// W0(W1 + W2)
// 				double cur_eb = MIN(max_eb_to_keep_sign_online_W(u1v2, u2v1), max_eb_to_keep_sign_online_W1_W2(u0v1, u1v0, u2v0, u0v2));
// 				eb = MAX(eb, cur_eb);				
// 			}
// 			if(!f3){
// 				// W2(W0 + W1)
// 				double cur_eb = MIN(max_eb_to_keep_sign_online_W(u0v1, u1v0), max_eb_to_keep_sign_online_W0_W1(u1v2, u2v1, u2v0, u0v2));
// 				eb = MAX(eb, cur_eb);				
// 			}
// 			eb = MIN(eb, max_pwr_eb);
// 			if((x>=6 && x<506) && (y>=6 && y<506) && (eb == 0)) {
// 				printf("%d %d, eb = 0\n", x, y);
// 			}
// 		}
// 	}
// 	if((x>=6 && x<506) && (y>=6 && y<506) && (eb == 0)) {
// 		printf("%d %d, eb = 0\n", x, y);
// 	}
// 	return eb;
// }

// template<typename T>
// void
// copy_row(T * dst, const T * U, const T * V, size_t n){
// 	for(int i=0; i<n; i++){
// 		dst[2*i] = U[i];
// 		dst[2*i+1] = V[i];
// 	}
// }

template<typename T>
T *
log_transform(const T * data, unsigned char * sign, size_t n){
	T * log_data = (T *) malloc(n*sizeof(T));
	for(int i=0; i<n; i++){
		sign[i] = 0;
		if(data[i] != 0){
			sign[i] = (data[i] > 0);
			log_data[i] = (data[i] > 0) ? log2(data[i]) : log2(-data[i]); 
		}
		else{
			sign[i] = 0;
			log_data[i] = 0; //TODO???
		}
	}
	return log_data;
}

// template<typename T>
// unsigned char *
// sz_compress_cp_preserve_2d_log(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb){
// 	size_t num_elements = r1 * r2;
// 	unsigned char * compressed = (unsigned char *) malloc(num_elements * 4 * sizeof(T));
// 	unsigned char * compressed_pos = compressed;
// 	size_t encode_element = (r1 - 1) * r2;
// 	// size_t sign_map_size = (encode_element - 1)/8 + 1;
// 	unsigned char * sign_map = (unsigned char *) malloc(num_elements);
// 	// Note the convert function has address auto increment
// 	T * log_U = log_transform(U, sign_map, num_elements);
// 	convertIntArray2ByteArray_fast_1b_to_result_sz(sign_map + r2, encode_element, compressed_pos);
// 	T * log_V = log_transform(V, sign_map, num_elements);
// 	convertIntArray2ByteArray_fast_1b_to_result_sz(sign_map + r2, encode_element, compressed_pos);
// 	free(sign_map);
// 	// for(int i=0; i<10; i++){
// 	// 	printf("%.4g %4g %.4g %4g\n", U[i], V[i], log_U[i], log_V[i]);
// 	// }
// 	// for(int i=0; i<10; i++){
// 	// 	printf("%.4g %4g %.4g %4g\n", U[r2+i], V[r2+i], log_U[r2+i], log_V[r2+i]);
// 	// }
// 	// if(transpose){
// 	// 	transpose_2d(U);
// 	// 	transpose_2d(V);
// 	// }
// 	printf("num_elements = %ld, encode_element = %ld\n", num_elements, encode_element);
// 	printf("sign recorded: %ld\n", compressed_pos - compressed);
// 	T * log_U_pos = log_U;
// 	T * log_V_pos = log_V;
// 	// 1 eb for 2 data (u, v)
// 	int * eb_quant_index = (int *) malloc(encode_element*sizeof(int));
// 	int * data_quant_index = (int *) malloc(2*encode_element*sizeof(int));
// 	T * unpred_data = (T *) malloc(2*num_elements*sizeof(T));
// 	T * unpred_data_pos = unpred_data;
// 	int * eb_quant_index_pos = eb_quant_index;
// 	int * data_quant_index_pos = data_quant_index;
// 	T * buffer_top = (T *) malloc(2*r2*sizeof(T));
// 	T * buffer_bot = (T *) malloc(2*r2*sizeof(T));
// 	// copy first row as unpredicable data
// 	memcpy(unpred_data_pos, U, r2*sizeof(T));
// 	unpred_data_pos += r2;
// 	memcpy(unpred_data_pos, V, r2*sizeof(T));
// 	unpred_data_pos += r2;
// 	copy_row(buffer_top, log_U_pos, log_V_pos, r2);
// 	log_U_pos += r2, log_V_pos += r2;
// 	copy_row(buffer_bot, log_U_pos, log_V_pos, r2);
// 	log_U_pos += r2, log_V_pos += r2;
// 	// next, row by row
// 	int base = 4;
// 	double log_of_base = log2(base);
// 	int capacity = 65536;
// 	int intv_radius = (capacity >> 1);
// 	unpred_vec<T> unpred_data_top;
// 	unpred_vec<T> unpred_data_bot;
// 	// printf("start compression: %ld\n", compressed_pos - compressed);
// 	for(int i=1; i<r1; i++){
// 		// printf("i = %d, log_U_offset = %ld, %ld %ld %ld\n", i, log_U_pos - log_U, eb_quant_index_pos - eb_quant_index, data_quant_index_pos - data_quant_index, unpred_data_pos - unpred_data);
// 		// for(int j=0; j<r2; j++){
// 		// 	printf("%.4g ~ %.4g, %.4g ~ %.4g \t", buffer_top[2*j], log_U_pos[- 2*r2 + j], buffer_top[2*j + 1], log_V_pos[- 2*r2 + j]);
// 		// }
// 		// printf("\n");
// 		// for(int j=0; j<r2; j++){
// 		// 	printf("%.4g ~ %.4g, %.4g ~ %.4g \t", buffer_bot[2*j], log_U_pos[- r2 + j], buffer_bot[2*j + 1], log_V_pos[- r2 + j]);
// 		// }
// 		// cin.get();
// 		T * cur_top_pos = buffer_top;
// 		T * cur_bot_pos = buffer_bot;
// 		for(int j=0; j<r2; j++){
// 			// TODO: recover data to origin domain instead of log domain !!!!!!!
// 			double eb1 = (j==0) ? 1 : derive_cp_eb_for_positions_online(cur_bot_pos - 2, cur_top_pos, cur_bot_pos, max_pwr_eb);
// 			double eb2 = (j==r2-1) ? 1 : derive_cp_eb_for_positions_online(cur_top_pos, cur_top_pos + 2, cur_bot_pos, max_pwr_eb);
// 			double required_eb = MIN(eb1, eb2);
// 			required_eb = MIN(required_eb, max_pwr_eb);
// 			if(required_eb < 0){
// 				// roll back?
// 				// TODO
// 				printf("eb = %.4g\n", required_eb);
// 				exit(0);
// 			}
// 			double log_eb = log2(1 + required_eb);
// 			log_eb -= MAX(cur_bot_pos[0], cur_bot_pos[1]) * std::numeric_limits<T>::epsilon();
// 			log_eb = MAX(log_eb, 0);
// 			*(eb_quant_index_pos ++) = (log_eb) ? eb_exponential_quantize(log_eb, base, log_of_base) : 0;
// 			// *(eb_quant_index_pos ++) = eb_linear_quantize(log_eb, 1e-4);
// 			// compress the data
// 			if(log_eb > 0){
// 				bool unpred_flag = false;
// 				size_t offset = log_U_pos - log_U - r2 + j;
// 				T origin_data[2] = {U[offset], V[offset]};
// 				// printf("origin_data = %.4g, log_data = %.4g, cur_data = %.4g\n", origin_data[0], log_U_pos[j - r2], cur_bot_pos[0]);
// 				T log_decompressed[2];
// 				for(int k=0; k<2; k++){
// 					T pred = (j == 0) ? cur_top_pos[k] : cur_top_pos[k] + cur_bot_pos[-2 + k] - cur_top_pos[-2 + k];
// 					T cur_data = cur_bot_pos[k];
// 					double diff = cur_data - pred;
// 					double quant_diff = fabs(diff) / log_eb + 1;
// 					if(quant_diff < capacity){
// 						quant_diff = (diff > 0) ? quant_diff : -quant_diff;
// 						int quant_index = (int)(quant_diff/2) + intv_radius;
// 						*(data_quant_index_pos ++) = quant_index;
// 						log_decompressed[k] = pred + 2 * (quant_index - intv_radius) * log_eb; 
// 						// check original data
// 						{
// 							T decompressed = exp2(log_decompressed[k]);
// 							if(fabs((decompressed - fabs(origin_data[k])) / origin_data[k]) >= required_eb){
// 								unpred_flag = true;
// 								break;
// 							}
// 						}
// 					}
// 					else unpred_flag = true;
// 				}
// 				if(unpred_flag){
// 					// recover quant index
// 					data_quant_index_pos[-2] = 0;
// 					data_quant_index_pos[-1] = 0;
// 					unpred_data_bot.push_back(origin_data[0]);
// 					unpred_data_bot.push_back(origin_data[1]);
// 				}
// 				else{
// 					// assign decompressed data
// 					cur_bot_pos[0] = log_decompressed[0];
// 					cur_bot_pos[1] = log_decompressed[1];
// 				}
// 			}
// 			else{
// 				// record as unpredictable data
// 				// original row: (log_U_pos - log_U) / r2
// 				// original column: j
// 				size_t offset = log_U_pos - log_U - r2 + j;
// 				*(data_quant_index_pos ++) = 0;
// 				*(data_quant_index_pos ++) = 0;
// 				unpred_data_bot.push_back(U[offset]);
// 				unpred_data_bot.push_back(V[offset]);
// 				// log data remain unchanged
// 			}
// 			cur_top_pos += 2;
// 			cur_bot_pos += 2;
// 			// for(int j=0; j<r2; j++){
// 			// 	printf("%.4g ~ %.4g, %.4g ~ %.4g \t", buffer_bot[2*j], log_U_pos[- r2 + j], buffer_bot[2*j + 1], log_V_pos[- r2 + j]);
// 			// }
// 			// cin.get();
// 		}
// 		// commit top row and swap
// 		memcpy(unpred_data_pos, &unpred_data_top[0], unpred_data_top.size()*sizeof(T));
// 		unpred_data_pos += unpred_data_top.size();
// 		unpred_data_top = unpred_data_bot;
// 		unpred_data_bot.clear();
// 		swap_buffer_pointer(buffer_top, buffer_bot);
// 		if(i<r1 - 1){
// 			copy_row(buffer_bot, log_U_pos, log_V_pos, r2);
// 			log_U_pos += r2, log_V_pos += r2;
// 		}
// 	}
// 	free(log_U);
// 	free(log_V);
// 	// commit the remaining top row
// 	memcpy(unpred_data_pos, &unpred_data_top[0], unpred_data_top.size()*sizeof(T));
// 	unpred_data_pos += unpred_data_top.size();
// 	free(buffer_top);
// 	free(buffer_bot);
// 	printf("offsets eb_q, data_q, unpred: %ld %ld %ld\n", eb_quant_index_pos - eb_quant_index, data_quant_index_pos - data_quant_index, unpred_data_pos - unpred_data);
// 	write_variable_to_dst(compressed_pos, intv_radius);
// 	size_t unpredictable_count = unpred_data_pos - unpred_data;
// 	write_variable_to_dst(compressed_pos, unpredictable_count);
// 	write_array_to_dst(compressed_pos, unpred_data, unpredictable_count);	
// 	free(unpred_data);
// 	printf("before huffman: %ld\n", compressed_pos - compressed);
// 	Huffman_encode_tree_and_data(2*65536, eb_quant_index, encode_element, compressed_pos);
// 	printf("eb size: %ld\n", compressed_pos - compressed);
// 	free(eb_quant_index);
// 	Huffman_encode_tree_and_data(2*capacity, data_quant_index, 2*encode_element, compressed_pos);
// 	printf("data size: %ld\n", compressed_pos - compressed);
// 	free(data_quant_index);
// 	compressed_size = compressed_pos - compressed;
// 	printf("final: %ld\n", compressed_pos - compressed);
// 	return compressed;	
// }

// template
// unsigned char *
// sz_compress_cp_preserve_2d_log(const float * U, const float * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

// template<typename T>
// unsigned char *
// sz_compress_cp_preserve_2d_online(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb){
// 	size_t num_elements = r1 * r2;
// 	unsigned char * compressed = (unsigned char *) malloc(num_elements * 4 * sizeof(T));
// 	unsigned char * compressed_pos = compressed;
// 	size_t encode_element = (r1 - 1) * r2;
// 	printf("num_elements = %ld, encode_element = %ld\n", num_elements, encode_element);
// 	const T * U_pos = U;
// 	const T * V_pos = V;
// 	int * eb_quant_index = (int *) malloc(2*encode_element*sizeof(int));
// 	int * data_quant_index = (int *) malloc(2*encode_element*sizeof(int));
// 	T * unpred_data = (T *) malloc(2*num_elements*sizeof(T));
// 	T * unpred_data_pos = unpred_data;
// 	int * eb_quant_index_pos = eb_quant_index;
// 	int * data_quant_index_pos = data_quant_index;
// 	T * buffer_top = (T *) malloc(2*r2*sizeof(T));
// 	T * buffer_bot = (T *) malloc(2*r2*sizeof(T));
// 	// copy first row as unpredicable data
// 	memcpy(unpred_data_pos, U, r2*sizeof(T));
// 	unpred_data_pos += r2;
// 	memcpy(unpred_data_pos, V, r2*sizeof(T));
// 	unpred_data_pos += r2;
// 	copy_row(buffer_top, U_pos, V_pos, r2);
// 	U_pos += r2, V_pos += r2;
// 	copy_row(buffer_bot, U_pos, V_pos, r2);
// 	U_pos += r2, V_pos += r2;
// 	// next, row by row
// 	int base = 4;
// 	double log_of_base = log2(base);
// 	int capacity = 65536;
// 	int intv_radius = (capacity >> 1);
// 	unpred_vec<T> unpred_data_top;
// 	unpred_vec<T> unpred_data_bot;
// 	// printf("start compression: %ld\n", compressed_pos - compressed);
// 	T * decompressed = (T *) malloc(2*num_elements*sizeof(T));
// 	T * decompressed_pos = decompressed;
// 	memcpy(decompressed_pos, buffer_top, 2*r2*sizeof(T));
// 	decompressed_pos += 2*r2;
// 	for(int i=1; i<r1; i++){
// 		T * cur_top_pos = buffer_top;
// 		T * cur_bot_pos = buffer_bot;
// 		for(int j=0; j<r2; j++){
// 			double required_eb = 0;
// 			if((i>=6 && i<506) && (j>=6 && j<506)){
// 				if(i>=6 && j>=6){
// 					x = i, y = j;
// 				}
// 				double eb1 = (j==0) ? 1 : derive_cp_eb_for_positions_online(cur_bot_pos - 2, cur_top_pos, cur_bot_pos, max_pwr_eb);
// 				double eb2 = (j==r2-1) ? 1 : derive_cp_eb_for_positions_online(cur_top_pos, cur_top_pos + 2, cur_bot_pos, max_pwr_eb);
// 				required_eb = MIN(eb1, eb2);
// 				required_eb = MIN(required_eb, max_pwr_eb);
// 				// if((i>=6 && i<506) && (j>=6 && j<506) && (required_eb == 0)) {
// 				// 	printf("%d %d, eb = 0\n", i, j);
// 				// 	exit(0);
// 				// }
// 				if(required_eb < 0){
// 					// roll back?
// 					// TODO
// 					printf("eb = %.4g\n", required_eb);
// 					exit(0);
// 				}
// 			}
// 			// *(eb_quant_index_pos ++) = eb_linear_quantize(log_eb, 1e-4);
// 			// compress the data
// 			if(required_eb > 0){
// 				bool unpred_flag = false;
// 				size_t offset = - r2 + j;
// 				T origin_data[2] = {U_pos[offset], V_pos[offset]};
// 				// printf("origin_data = %.4g, log_data = %.4g, cur_data = %.4g\n", origin_data[0], U_pos[j - r2], cur_bot_pos[0]);
// 				T decompressed[2];
// 				for(int k=0; k<2; k++){
// 					double abs_eb = fabs(cur_bot_pos[k]) * required_eb;
// 					*(eb_quant_index_pos ++) = eb_exponential_quantize(abs_eb, base, log_of_base);
// 					T pred = (j == 0) ? cur_top_pos[k] : cur_top_pos[k] + cur_bot_pos[-2 + k] - cur_top_pos[-2 + k];
// 					T cur_data = cur_bot_pos[k];
// 					double diff = cur_data - pred;
// 					double quant_diff = fabs(diff) / abs_eb + 1;
// 					if(quant_diff < capacity){
// 						quant_diff = (diff > 0) ? quant_diff : -quant_diff;
// 						int quant_index = (int)(quant_diff/2) + intv_radius;
// 						*(data_quant_index_pos ++) = quant_index;
// 						decompressed[k] = pred + 2 * (quant_index - intv_radius) * abs_eb; 
// 						// check original data
// 						{
// 							if(fabs((decompressed[k] - origin_data[k]) / origin_data[k]) >= required_eb){
// 								unpred_flag = true;
// 								break;
// 							}
// 						}
// 					}
// 					else unpred_flag = true;
// 				}
// 				if(unpred_flag){
// 					// recover quant index
// 					data_quant_index_pos[-2] = 0;
// 					data_quant_index_pos[-1] = 0;
// 					unpred_data_bot.push_back(origin_data[0]);
// 					unpred_data_bot.push_back(origin_data[1]);
// 				}
// 				else{
// 					// assign decompressed data
// 					cur_bot_pos[0] = decompressed[0];
// 					cur_bot_pos[1] = decompressed[1];
// 				}
// 			}
// 			else{
// 				// record as unpredictable data
// 				// original row: (U_pos - U) / r2
// 				// original column: j
// 				size_t offset = - r2 + j;
// 				*(eb_quant_index_pos ++) = 0;
// 				*(eb_quant_index_pos ++) = 0;
// 				*(data_quant_index_pos ++) = 0;
// 				*(data_quant_index_pos ++) = 0;
// 				unpred_data_bot.push_back(U_pos[offset]);
// 				unpred_data_bot.push_back(V_pos[offset]);
// 			}
// 			cur_top_pos += 2;
// 			cur_bot_pos += 2;
// 		}
// 		// commit top row and swap
// 		memcpy(unpred_data_pos, &unpred_data_top[0], unpred_data_top.size()*sizeof(T));
// 		unpred_data_pos += unpred_data_top.size();
// 		unpred_data_top = unpred_data_bot;
// 		unpred_data_bot.clear();
// 		// for(int j=0; j<r2; j++){
// 		// 	printf("%.4g ~ %.4g, %.4g ~ %.4g \t", buffer_bot[2*j], U_pos[- r2 + j], buffer_bot[2*j + 1], V_pos[- r2 + j]);
// 		// }
// 		// cin.get();

// 		swap_buffer_pointer(buffer_top, buffer_bot);
// 		if(i<r1 - 1){
// 			copy_row(buffer_bot, U_pos, V_pos, r2);
// 			U_pos += r2, V_pos += r2;
// 		}
// 		memcpy(decompressed_pos, buffer_top, 2*r2*sizeof(T));
// 		decompressed_pos += 2*r2;
// 	}
// 	writefile("/Users/xin/github/ftk_xin/ftk/build/tmp.dat", decompressed, 2*num_elements);
// 	free(decompressed);
// 	// commit the remaining top row
// 	memcpy(unpred_data_pos, &unpred_data_top[0], unpred_data_top.size()*sizeof(T));
// 	unpred_data_pos += unpred_data_top.size();
// 	free(buffer_top);
// 	free(buffer_bot);
// 	printf("offsets eb_q, data_q, unpred: %ld %ld %ld\n", eb_quant_index_pos - eb_quant_index, data_quant_index_pos - data_quant_index, unpred_data_pos - unpred_data);
// 	write_variable_to_dst(compressed_pos, intv_radius);
// 	size_t unpredictable_count = unpred_data_pos - unpred_data;
// 	write_variable_to_dst(compressed_pos, unpredictable_count);
// 	write_array_to_dst(compressed_pos, unpred_data, unpredictable_count);	
// 	free(unpred_data);
// 	printf("before huffman: %ld\n", compressed_pos - compressed);
// 	Huffman_encode_tree_and_data(2*65536, eb_quant_index, 2*encode_element, compressed_pos);
// 	printf("eb size: %ld\n", compressed_pos - compressed);
// 	free(eb_quant_index);
// 	Huffman_encode_tree_and_data(2*capacity, data_quant_index, 2*encode_element, compressed_pos);
// 	printf("data size: %ld\n", compressed_pos - compressed);
// 	free(data_quant_index);
// 	compressed_size = compressed_pos - compressed;
// 	printf("final: %ld\n", compressed_pos - compressed);
// 	writefile("tmp.dat", compressed, compressed_pos - compressed);
// 	return compressed;	
// }

// template
// unsigned char *
// sz_compress_cp_preserve_2d_online(const float * U, const float * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

// maximal error bound to keep the sign of postive*(1+e)^d - negative*(1-e)^d
template<typename T>
inline double max_eb_to_keep_sign(const T positive, const T negative, int degree){
  if((negative < 0) || (positive < 0)){
    printf("%.4f, %.4f\n", negative, positive);
    exit(0);
  }
  if((negative < std::numeric_limits<double>::epsilon()) || (positive < std::numeric_limits<double>::epsilon())){
    return 1;
  }
  double c = 0;
  switch(degree){
    case 1:
      c = positive / negative;
      break;
    case 2:
      c = sqrt(positive / negative);
      break;
    case 3:
      c = cbrt(positive / negative);
      break;
    default:
      printf("Degree higher than 3 not supported yet\n");
      exit(0);
  }
  return MIN(1, fabs(c - 1) / (c + 1));  
}

// maximal error bound to keep the sign of u0v1 - u0v2 + u1v2 - u1v0 + u2v0 - u2v1
template<typename T>
inline double max_eb_to_keep_sign_det2x2(const T u0, const T u1, const T u2, const T v0, const T v1, const T v2){
  double positive = 0;
  double negative = 0;
  accumulate(u0*v1, positive, negative);
  accumulate(-u0*v2, positive, negative);
  accumulate(u1*v2, positive, negative);
  accumulate(-u1*v0, positive, negative);
  accumulate(u2*v0, positive, negative);
  accumulate(-u2*v1, positive, negative);
  return max_eb_to_keep_sign(positive, negative, 2);
}

// maximal error bound to keep the sign of u0v1 - u1v0
template<typename T>
inline double max_eb_to_keep_sign_2(const T u0, const T u1, const T v0, const T v1){
  double positive = 0;
  double negative = 0;
  accumulate(u0*v1, positive, negative);
  accumulate(-u1*v0, positive, negative);
  return max_eb_to_keep_sign(positive, negative, 2);
}

// maximal error bound to keep the sign of u0v1 - u1v0 + u1v2 - u2v1
template<typename T>
inline double max_eb_to_keep_sign_4(const T u0, const T u1, const T u2, const T v0, const T v1, const T v2){
  double positive = 0;
  double negative = 0;
  accumulate(u0*v1, positive, negative);
  accumulate(-u1*v0, positive, negative);
  accumulate(u1*v2, positive, negative);
  accumulate(-u2*v1, positive, negative);
  return max_eb_to_keep_sign(positive, negative, 2);
}

// det(c) = (x0 - x2)*(y1 - y2) - (x1 - x2)*(y0 - y2)
// c0 = (y1 - y2) / det(c)   c1 = -(y0 - y2) / det(c)
// c1 = -(x1 - x2) / det(c)  c3 = (x0 - x2) / det(c)
template<typename T>
inline void get_adjugate_matrix_for_position(const T x0, const T x1, const T x2, const T y0, const T y1, const T y2, T c[4]){
  T determinant = (x0 - x2)*(y1 - y2) - (x1 - x2)*(y0 - y2);
  c[0] = (y1 - y2) / determinant;
  c[1] = -(y0 - y2) / determinant;
  c[2] = -(x1 - x2) / determinant;
  c[3] = (x0 - x2) / determinant;
  // printf("%.4g, %.2g %.2g %.2g %.2g\n", determinant, c[0], c[1], c[2], c[3]);
  // exit(0);
}

// accumulate positive and negative in (a + b + c ...)^2
template<typename T>
inline void accumulate_in_square(const std::vector<T>& coeff, double& positive, double& negative){
  for(int i=0; i<coeff.size(); i++){
    for(int j=0; j<coeff.size(); j++){
      accumulate(coeff[i]*coeff[j], positive, negative);
    }
  }
}

// maximal error bound to keep the sign of B^2 - 4C
// where  B = - (c0 * (u0 - u2) - c1 * (u1 - u2) + c2 * (v0 - v2) - c3 * (v1 - v2))
//        C = det2x2 = u0v1 - u0v2 + u1v2 - u1v0 + u2v0 - u2v1
template<typename T>
inline double max_eb_to_keep_sign_eigen_delta_2(const T u0, const T u1, const T u2, const T v0, const T v1, const T v2,
  const T x0, const T x1, const T x2, const T y0, const T y1, const T y2){
  double eb = 1;
  T c[4] = {0};
  {
    get_adjugate_matrix_for_position(x0, x1, x2, y0, y1, y2, c);
    // keep sign for B
    double positive = 0;
    double negative = 0;
    accumulate(c[0]*u0, positive, negative);
    accumulate(c[1]*u1, positive, negative);
    accumulate(-(c[0] + c[1])*u2, positive, negative);
    accumulate(c[2]*v0, positive, negative);
    accumulate(c[3]*v1, positive, negative);
    accumulate(-(c[2] + c[3])*v2, positive, negative);
    eb = max_eb_to_keep_sign(positive, negative, 1);
    // keep sign for C
    eb = MIN(eb, max_eb_to_keep_sign_det2x2(u0, u1, u2, v0, v1, v2));
  }
  T m = c[1]*c[2] - c[0]*c[3];
  T C = (-m) * (u0*v1 - u0*v2 + u1*v2 - u1*v0 + u2*v0 - u2*v1);
  if(C <= 0) return eb;
  {
    std::vector<T> coeff(6);
    coeff[0] = c[0]*u0;
    coeff[1] = c[1]*u1;
    coeff[2] = - (c[1] + c[0])*u2;
    coeff[3] = c[2]*v0;
    coeff[4] = c[3]*v1;
    coeff[5] = - (c[3] + c[2])*v2;
    // keep sign for B^2 - 4*C
    double positive = 0;
    double negative = 0;
    accumulate_in_square(coeff, positive, negative);
    accumulate(-4*m*u1*v0, positive, negative);
    accumulate(4*m*u2*v0, positive, negative);
    accumulate(4*m*u0*v1, positive, negative);
    accumulate(-4*m*u2*v1, positive, negative);
    accumulate(-4*m*u0*v2, positive, negative);
    accumulate(4*m*u1*v2, positive, negative);
    eb = MIN(eb, max_eb_to_keep_sign(positive, negative, 2));
  }
  return eb;
}

template<typename T>
inline double max_eb_to_keep_position_and_type(const T u0, const T u1, const T u2, const T v0, const T v1, const T v2,
											const T x0, const T x1, const T x2, const T y0, const T y1, const T y2){
	double u0v1 = u0 * v1;
	double u1v0 = u1 * v0;
	double u0v2 = u0 * v2;
	double u2v0 = u2 * v0;
	double u1v2 = u1 * v2;
	double u2v1 = u2 * v1;
	double det = u0v1 - u1v0 + u1v2 - u2v1 + u2v0 - u0v2;
	double eb = 0;
	if(det != 0){
		bool f1 = (det / (u2v0 - u0v2) >= 1);
		bool f2 = (det / (u1v2 - u2v1) >= 1); 
		bool f3 = (det / (u0v1 - u1v0) >= 1); 
		if(f1 && f2 && f3){
			// critical point found
			eb = 1;
			eb = MIN(max_eb_to_keep_sign_2(u2, u0, v2, v0), eb);
			eb = MIN(max_eb_to_keep_sign_2(u1, u2, v1, v2), eb);
			eb = MIN(max_eb_to_keep_sign_2(u0, u1, v0, v1), eb);
			eb = MIN(max_eb_to_keep_sign_4(u0, u1, u2, v0, v1, v2), eb);
			eb = MIN(max_eb_to_keep_sign_4(u2, u0, u1, v2, v0, v1), eb);
			eb = MIN(max_eb_to_keep_sign_4(u1, u2, u0, v1, v2, v0), eb);
			eb = MIN(max_eb_to_keep_sign_eigen_delta_2(u0, u1, u2, v0, v1, v2, x0, x1, x2, y0, y1, y2), eb);
		}
		else{
			// no critical point
			eb = 0;
			if(!f1){
				double eb_cur = MIN(max_eb_to_keep_sign_2(u2, u0, v2, v0), max_eb_to_keep_sign_4(u0, u1, u2, v0, v1, v2));
				eb = MAX(eb, eb_cur);
			}
			if(!f2){
				double eb_cur = MIN(max_eb_to_keep_sign_2(u1, u2, v1, v2), max_eb_to_keep_sign_4(u2, u0, u1, v2, v0, v1));
				eb = MAX(eb, eb_cur);
			}
			if(!f3){
				double eb_cur = MIN(max_eb_to_keep_sign_2(u0, u1, v0, v1), max_eb_to_keep_sign_4(u1, u2, u0, v1, v2, v0));
				eb = MAX(eb, eb_cur);
			}
			eb = MIN(eb, DEFAULT_EB);
		}
	}
	return eb;
}

// compression with pre-computed error bounds
template<typename T>
unsigned char *
sz_compress_cp_preserve_2d_offline(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb){

	size_t num_elements = r1 * r2;
	double * eb = (double *) malloc(num_elements * sizeof(double));
	for(int i=0; i<num_elements; i++) eb[i] = 1;
	const T * U_pos = U;
	const T * V_pos = V;
	double * eb_pos = eb;
	// coordinates for triangle_coordinates
	const T X_upper[3][2] = {{0, 0}, {1, 0}, {1, 1}};
	const T X_lower[3][2] = {{0, 0}, {0, 1}, {1, 1}};
	const size_t offset_upper[3] = {0, r2, r2+1};
	const size_t offset_lower[3] = {0, 1, r2+1};
	printf("compute eb\n");
	for(int i=0; i<r1-1; i++){
		const T * U_row_pos = U_pos;
		const T * V_row_pos = V_pos;
		double * eb_row_pos = eb_pos;
		for(int j=0; j<r2-1; j++){
			for(int k=0; k<2; k++){
				auto X = (k == 0) ? X_upper : X_lower;
				auto offset = (k == 0) ? offset_upper : offset_lower;
				double max_cur_eb = max_eb_to_keep_position_and_type(U_row_pos[offset[0]], U_row_pos[offset[1]], U_row_pos[offset[2]],
					V_row_pos[offset[0]], V_row_pos[offset[1]], V_row_pos[offset[2]], X[0][0], X[1][0], X[2][0],
					X[0][1], X[1][1], X[2][1]);
				eb_row_pos[offset[0]] = MIN(eb_row_pos[offset[0]], max_cur_eb);
				eb_row_pos[offset[1]] = MIN(eb_row_pos[offset[1]], max_cur_eb);
				eb_row_pos[offset[2]] = MIN(eb_row_pos[offset[2]], max_cur_eb);
			}
			U_row_pos ++;
			V_row_pos ++;
			eb_row_pos ++;
		}
		U_pos += r2;
		V_pos += r2;
		eb_pos += r2;
	}
	printf("compute eb done\n");
	double * eb_u = (double *) malloc(num_elements * sizeof(double));
	double * eb_v = (double *) malloc(num_elements * sizeof(double));
	int * eb_quant_index = (int *) malloc(2*num_elements*sizeof(int));
	int * eb_quant_index_pos = eb_quant_index;
	const int base = 4;
	double log2_of_base = log2(base);
	const double threshold = std::numeric_limits<float>::epsilon();
	for(int i=0; i<num_elements; i++){
		eb_u[i] = fabs(U[i]) * eb[i];
		*(eb_quant_index_pos ++) = eb_exponential_quantize(eb_u[i], base, log2_of_base);
		if(eb_u[i] < threshold) eb_u[i] = 0;
	}
	for(int i=0; i<num_elements; i++){
		eb_v[i] = fabs(V[i]) * eb[i];
		*(eb_quant_index_pos ++) = eb_exponential_quantize(eb_v[i], base, log2_of_base);
		if(eb_v[i] < threshold) eb_v[i] = 0;
	}
	free(eb);
	printf("quantize eb done\n");
	unsigned char * compressed_eb = (unsigned char *) malloc(2*num_elements*sizeof(int));
	unsigned char * compressed_eb_pos = compressed_eb; 
	Huffman_encode_tree_and_data(2*256, eb_quant_index, 2*num_elements, compressed_eb_pos);
	size_t compressed_eb_size = compressed_eb_pos - compressed_eb;
	size_t compressed_u_size = 0;
	size_t compressed_v_size = 0;
	unsigned char * compressed_u = sz_compress_2d_with_eb(U, eb_u, r1, r2, compressed_u_size);
	unsigned char * compressed_v = sz_compress_2d_with_eb(V, eb_v, r1, r2, compressed_v_size);
	printf("eb_size = %ld, u_size = %ld, v_size = %ld\n", compressed_eb_size, compressed_u_size, compressed_v_size);
	free(eb_u);
	free(eb_v);
	compressed_size = sizeof(int) + sizeof(size_t) + compressed_eb_size + sizeof(size_t) + compressed_u_size + sizeof(size_t) + compressed_v_size;
	unsigned char * compressed = (unsigned char *) malloc(compressed_size);
	unsigned char * compressed_pos = compressed;
	write_variable_to_dst(compressed_pos, base);
	write_variable_to_dst(compressed_pos, compressed_eb_size);
	write_variable_to_dst(compressed_pos, compressed_u_size);
	write_variable_to_dst(compressed_pos, compressed_v_size);
	write_array_to_dst(compressed_pos, compressed_eb, compressed_eb_size);
	write_array_to_dst(compressed_pos, compressed_u, compressed_u_size);
	printf("compressed_pos = %ld\n", compressed_pos - compressed);
	write_array_to_dst(compressed_pos, compressed_v, compressed_v_size);
	free(compressed_eb);
	free(compressed_u);
	free(compressed_v);
	return compressed;
}

template
unsigned char *
sz_compress_cp_preserve_2d_offline(const float * U, const float * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

template
unsigned char *
sz_compress_cp_preserve_2d_offline(const double * U, const double * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

// compression with pre-computed error bounds in logarithmic domain
template<typename T>
unsigned char *
sz_compress_cp_preserve_2d_offline_log(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb){

	size_t num_elements = r1 * r2;
	double * eb = (double *) malloc(num_elements * sizeof(double));
	for(int i=0; i<num_elements; i++) eb[i] = 1;
	const T * U_pos = U;
	const T * V_pos = V;
	double * eb_pos = eb;
	// coordinates for triangle_coordinates
	const T X_upper[3][2] = {{0, 0}, {1, 0}, {1, 1}};
	const T X_lower[3][2] = {{0, 0}, {0, 1}, {1, 1}};
	const size_t offset_upper[3] = {0, r2, r2+1};
	const size_t offset_lower[3] = {0, 1, r2+1};
	printf("compute eb\n");
	for(int i=0; i<r1-1; i++){
		const T * U_row_pos = U_pos;
		const T * V_row_pos = V_pos;
		double * eb_row_pos = eb_pos;
		for(int j=0; j<r2-1; j++){
			for(int k=0; k<2; k++){
				auto X = (k == 0) ? X_upper : X_lower;
				auto offset = (k == 0) ? offset_upper : offset_lower;
				double max_cur_eb = max_eb_to_keep_position_and_type(U_row_pos[offset[0]], U_row_pos[offset[1]], U_row_pos[offset[2]],
					V_row_pos[offset[0]], V_row_pos[offset[1]], V_row_pos[offset[2]], X[0][0], X[1][0], X[2][0],
					X[0][1], X[1][1], X[2][1]);
				eb_row_pos[offset[0]] = MIN(eb_row_pos[offset[0]], max_cur_eb);
				eb_row_pos[offset[1]] = MIN(eb_row_pos[offset[1]], max_cur_eb);
				eb_row_pos[offset[2]] = MIN(eb_row_pos[offset[2]], max_cur_eb);
			}
			U_row_pos ++;
			V_row_pos ++;
			eb_row_pos ++;
		}
		U_pos += r2;
		V_pos += r2;
		eb_pos += r2;
	}
	printf("compute eb done\n");
	size_t sign_map_size = (num_elements - 1)/8 + 1;
	unsigned char * sign_map_compressed = (unsigned char *) malloc(2*sign_map_size);
	unsigned char * sign_map_compressed_pos = sign_map_compressed;
	unsigned char * sign_map = (unsigned char *) malloc(num_elements*sizeof(unsigned char));
	// Note the convert function has address auto increment
	T * log_U = log_transform(U, sign_map, num_elements);
	convertIntArray2ByteArray_fast_1b_to_result_sz(sign_map, num_elements, sign_map_compressed_pos);
	T * log_V = log_transform(V, sign_map, num_elements);
	convertIntArray2ByteArray_fast_1b_to_result_sz(sign_map, num_elements, sign_map_compressed_pos);
	free(sign_map);
	// transfrom eb to log(1 + eb) and the quantize
	int * eb_quant_index = (int *) malloc(num_elements*sizeof(int));
	int * eb_quant_index_pos = eb_quant_index;
	const int base = 2;
	double log2_of_base = log2(base);
	const double threshold = std::numeric_limits<float>::epsilon();
	for(int i=0; i<num_elements; i++){
		eb[i] = log2(1 + eb[i]);
		*(eb_quant_index_pos ++) = eb_exponential_quantize(eb[i], base, log2_of_base);
		if(eb[i] < threshold) eb[i] = 0;
	}
	printf("quantize eb done\n");
	unsigned char * compressed_eb = (unsigned char *) malloc(num_elements*sizeof(int));
	unsigned char * compressed_eb_pos = compressed_eb; 
	Huffman_encode_tree_and_data(2*256, eb_quant_index, num_elements, compressed_eb_pos);
	size_t compressed_eb_size = compressed_eb_pos - compressed_eb;
	size_t compressed_u_size = 0;
	size_t compressed_v_size = 0;
	unsigned char * compressed_u = sz_compress_2d_with_eb(log_U, eb, r1, r2, compressed_u_size);
	free(log_U);
	unsigned char * compressed_v = sz_compress_2d_with_eb(log_V, eb, r1, r2, compressed_v_size);
	free(log_V);
	printf("eb_size = %ld, log_u_size = %ld, log_v_size = %ld\n", compressed_eb_size, compressed_u_size, compressed_v_size);
	free(eb);
	compressed_size = sizeof(int) + 2*sign_map_size + sizeof(size_t) + compressed_eb_size + sizeof(size_t) + compressed_u_size + sizeof(size_t) + compressed_v_size;
	unsigned char * compressed = (unsigned char *) malloc(compressed_size);
	unsigned char * compressed_pos = compressed;
	write_variable_to_dst(compressed_pos, base);
	write_variable_to_dst(compressed_pos, compressed_eb_size);
	write_variable_to_dst(compressed_pos, compressed_u_size);
	write_variable_to_dst(compressed_pos, compressed_v_size);
	write_array_to_dst(compressed_pos, compressed_eb, compressed_eb_size);
	write_array_to_dst(compressed_pos, sign_map_compressed, 2*sign_map_size);
	write_array_to_dst(compressed_pos, compressed_u, compressed_u_size);
	write_array_to_dst(compressed_pos, compressed_v, compressed_v_size);
	free(sign_map_compressed);
	free(compressed_eb);
	free(compressed_u);
	free(compressed_v);
	return compressed;
}

template
unsigned char *
sz_compress_cp_preserve_2d_offline_log(const float * U, const float * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

template
unsigned char *
sz_compress_cp_preserve_2d_offline_log(const double * U, const double * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);
