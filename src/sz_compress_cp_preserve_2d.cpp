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
			log_data[i] = -100; //TODO???
		}
	}
	return log_data;
}

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
// where  B = - (c0 * (u0 - u2) + c1 * (u1 - u2) + c2 * (v0 - v2) + c3 * (v1 - v2))
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

// maximal error bound to keep the sign of W_i 
// e.g W0 = u1v2 - u2v1
// where u2, v2 is current data, u1, v1 is decompressed data
template<typename T>
inline double max_eb_to_keep_sign_online_W(const T u1v2, const T u2v1){
	if(u1v2 * u2v1 <= 0) return 1;
	double c = u1v2 / u2v1;
	return fabs(c - 1) / (c + 1);
}

// maximal error bound to keep the sign of W_i(W_(i+1)+W_(i+2))
// which preserves the FP/FN for positions
// W0 + W1 = u1v2 - u2v1 + u2v0 - u0v2
// W1 + W2 = u2v0 - u0v2 + u0v1 - u1v0
// W2 + W0 = u0v1 - u1v0 + u1v2 - u2v1
// where u2, v2 is current data, u0, v0, u1, v1 is decompressed data

// maximal error bound to keep the sign of A*(1 + e_1) + B*(1 + e_2) + C
template<typename T>
inline double max_eb_to_keep_sign_online(const T A, const T B, const T C){
	if(A * B > 0){
		// same sign
		T c = -C / (A + B);
		return (c < 1) ? 1 - c : c - 1;
	}
	else{
		T c = (A + B + C) / (A - B);
		return fabs(c);
	}
}

// W1 + W2
template<typename T>
inline double max_eb_to_keep_sign_online_W1_W2(const T u0v1, const T u1v0, const T u2v0, const T u0v2){
	if(u1v0 - u0v1 == 0){
		double positive = 0;
		double negative = 0;
		accumulate(u2v0, positive, negative);
		accumulate(-u0v2, positive, negative);
		return max_eb_to_keep_sign(positive, negative, 1);		
	}
	if(u2v0 - u0v2 == 0) return 1;
	return max_eb_to_keep_sign_online(u2v0, -u0v2, u0v1 - u1v0);
	// if(u2v0 * u0v2 < 0){ 
	// 	T c = (u1v0 - u0v1) / (u2v0 - u0v2);
	// 	if(c < 0) return 1;
	// 	return (u2v0 - u0v2 > 0) ? 1 - c : c - 1;
	// }
	// else{
	// 	T c = (u2v0 - u0v2 + u0v1 - u1v0) / (u2v0 + u0v2);
	// 	return fabs(c);
	// }
}

// W0 + W2
template<typename T>
inline double max_eb_to_keep_sign_online_W0_W2(const T u0v1, const T u1v0, const T u1v2, const T u2v1){
	if(u1v0 - u0v1 == 0){
		double positive = 0;
		double negative = 0;
		accumulate(u1v2, positive, negative);
		accumulate(-u2v1, positive, negative);
		return max_eb_to_keep_sign(positive, negative, 1);		
	}
	if(u1v2 - u2v1 == 0) return 1;
	return max_eb_to_keep_sign_online(u1v2, -u2v1, u0v1 - u1v0);
	// if(u1v2 - u2v1 < 0){ 
	// 	T c = (u1v0 - u0v1) / (u1v2 - u2v1);
	// 	if(c < 0) return 1;
	// 	return (u1v2 - u2v1 > 0) ? 1 - c : c - 1;
	// }
	// else{
	// 	T c = (u1v2 - u2v1 + u0v1 - u1v0) / (u1v2 + u2v1);
	// 	return fabs(c);
	// }
}

// W0 + W1
template<typename T>
inline double max_eb_to_keep_sign_online_W0_W1(const T u1v2, const T u2v1, const T u2v0, const T u0v2){
	double positive = 0;
	double negative = 0;
	accumulate(u1v2, positive, negative);
	accumulate(-u2v1, positive, negative);
	accumulate(u2v0, positive, negative);
	accumulate(-u0v2, positive, negative);
	if(positive * negative <= 0) return 1;
	double c = positive / negative;
	return MIN(1, fabs(c - 1) / (c + 1));
}

template<typename T>
inline double max_eb_to_keep_position_online(const T u0v1, const T u1v0, const T u1v2, const T u2v1, const T u2v0, const T u0v2){
	// for W_i
	double eb = MIN(max_eb_to_keep_sign_online_W(u1v2, u2v1), max_eb_to_keep_sign_online_W(u2v0, u0v2));
	// for W1 + W2, W2 + W0
	eb = MIN(eb, MIN(max_eb_to_keep_sign_online_W1_W2(u0v1, u1v0, u2v0, u0v2), max_eb_to_keep_sign_online_W0_W2(u0v1, u1v0, u1v2, u2v1)));
	// for W0 + W1
	eb = MIN(eb, max_eb_to_keep_sign_online_W0_W1(u1v2, u2v1, u2v0, u0v2));
	return eb;
}

// maximal error bound to keep the sign of B^2 - 4C
// where  B = - (c0 * (u0 - u2) + c1 * (u1 - u2) + c2 * (v0 - v2) + c3 * (v1 - v2))
//        C = det2x2 = u0v1 - u0v2 + u1v2 - u1v0 + u2v0 - u2v1
template<typename T>
inline double max_eb_to_keep_type_online(const T u0, const T u1, const T u2, const T v0, const T v1, const T v2, const T x0, const T x1, const T x2, const T y0, const T y1, const T y2){
	double eb = 1;
	T c[4] = {0};
	{
		get_adjugate_matrix_for_position(x0, x1, x2, y0, y1, y2, c);
		// keep sign for B
	    // coeff[0] = c[0]*u0;
	    // coeff[1] = c[1]*u1;
	    // coeff[2] = - (c[1] + c[0])*u2;
	    // coeff[3] = c[2]*v0;
	    // coeff[4] = c[3]*v1;
	    // coeff[5] = - (c[3] + c[2])*v2;
		eb = max_eb_to_keep_sign_online(-c[0]*u2 - c[1]*u2, -c[2]*v2 - c[3]*v2, c[0]*u0 + c[1]*u1 + c[2]*v0 + c[3]*v1);
		// keep sign for C
		eb = MIN(eb, max_eb_to_keep_sign_online(u2*v0 - u2*v1, u1*v2 - u0*v2, u0*v1 - u1*v0));
	}
	T m = c[1]*c[2] - c[0]*c[3];
	T C = (-m) * (u0*v1 - u0*v2 + u1*v2 - u1*v0 + u2*v0 - u2*v1);
	if(C <= 0) return eb;
	{
		// Note that meaning of B in the rhs changes here
		// keep sign for B^2 - 4*C
		// B = A*(1+e_1) + B*(1+e_2) + C
		// C = D*(1+e_1) + E*(1+e_2) + F
		double A = -c[0]*u2 - c[1]*u2, B = -c[2]*v2 - c[3]*v2, C = c[0]*u0 + c[1]*u1 + c[2]*v0 + c[3]*v1;
		double D = (-m)*(u2*v0 - u2*v1), E = (-m)*(u1*v2 - u0*v2), F = (-m)*(u0*v1 - u1*v0);
		// B = A*e_1 + B*e_2 + C'
		// C = D*e_1 + E*e_2 + F'
		C += A + B, F += D + E;
		// B^2 - 4C = (A*e_1 + B*e_2)^2 + (2AC' - 4D)e_1 + (2BC' - 4E)e_2 + C'^2 - 4F'
		double delta = C*C - 4*F;
		if(delta == 0) return 0;
		else if(delta > 0){
			// (|2AC' - 4D| + |2BC' - 4E|)* -e + delta > 0
			eb = MIN(eb, delta/(fabs(2*A*C - 4*D) + fabs(2*B*C - 4*E)));
		}
		else{
			// (|A| + |B|)*e^2 + (|2AC' - 4D| + |2BC' - 4E|)*e + delta < 0
			double a = fabs(A) + fabs(B);
			double b = fabs(2*A*C - 4*D) + fabs(2*B*C - 4*E);
			double c = delta;
			if(b*b - 4*a*c < 0){
				printf("impossible as a*c is always less than 0\n");
				exit(0);
			}
			eb = MIN(eb, (-b + sqrt(b*b - 4*a*c))/(2*a));
		}
	}
	return eb;
}

/*
triangle mesh x0, x1, x2, derive cp-preserving eb for x2 given x0, x1
*/
template<typename T>
double 
derive_cp_eb_for_positions_online(const T u0, const T u1, const T u2, const T v0, const T v1, const T v2,
								const T x0, const T x1, const T x2, const T y0, const T y1, const T y2){
	double u1v2 = u1*v2;
	double u2v1 = u2*v1;
	double u2v0 = u2*v0;
	double u0v2 = u0*v2;
	double u0v1 = u0*v1;
	double u1v0 = u1*v0;
	double eb = 0;
    double det = u0v1 - u1v0 + u1v2 - u2v1 + u2v0 - u0v2;
    if(det != 0){
		bool f1 = (det / (u2v0 - u0v2) >= 1);
		bool f2 = (det / (u1v2 - u2v1) >= 1); 
		bool f3 = (det / (u0v1 - u1v0) >= 1); 
		if(f1 && f2 && f3){
			// eb = max_eb_to_keep_position_online(u0v1, u1v0, u1v2, u2v1, u2v0, u0v2);
			eb = MIN(max_eb_to_keep_position_online(u0v1, u1v0, u1v2, u2v1, u2v0, u0v2), 
				max_eb_to_keep_type_online(u0, u1, u2, v0, v1, v2, x0, x1, x2, y0, y1, y2));
		}
		else{
			eb = 0;
			if(!f1){
				// W1(W0 + W2)
				double cur_eb = MIN(max_eb_to_keep_sign_online_W(u2v0, u0v2), max_eb_to_keep_sign_online_W0_W2(u0v1, u1v0, u1v2, u2v1));
				eb = MAX(eb, cur_eb);
			}
			if(!f2){
				// W0(W1 + W2)
				double cur_eb = MIN(max_eb_to_keep_sign_online_W(u1v2, u2v1), max_eb_to_keep_sign_online_W1_W2(u0v1, u1v0, u2v0, u0v2));
				eb = MAX(eb, cur_eb);				
			}
			if(!f3){
				// W2(W0 + W1)
				double cur_eb = MIN(max_eb_to_keep_sign_online_W(u0v1, u1v0), max_eb_to_keep_sign_online_W0_W1(u1v2, u2v1, u2v0, u0v2));
				eb = MAX(eb, cur_eb);				
			}
		}
	}
	return eb;
}

template<typename T>
inline bool 
inbound(T index, T lb, T ub){
	return (index >= lb) && (index < ub);
}

template<typename T>
unsigned char *
sz_compress_cp_preserve_2d_online(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb){
	size_t num_elements = r1 * r2;
	T * decompressed_U = (T *) malloc(num_elements*sizeof(T));
	memcpy(decompressed_U, U, num_elements*sizeof(T));
	T * decompressed_V = (T *) malloc(num_elements*sizeof(T));
	memcpy(decompressed_V, V, num_elements*sizeof(T));
	int * eb_quant_index = (int *) malloc(2*num_elements*sizeof(int));
	int * data_quant_index = (int *) malloc(2*num_elements*sizeof(int));
	int * eb_quant_index_pos = eb_quant_index;
	int * data_quant_index_pos = data_quant_index;
	// next, row by row
	const int base = 4;
	const double log_of_base = log2(base);
	const int capacity = 65536;
	const int intv_radius = (capacity >> 1);
	unpred_vec<T> unpred_data;
	T * U_pos = decompressed_U;
	T * V_pos = decompressed_V;
	// offsets to get six adjacent triangle indices
	// the 7-th rolls back to T0
	/*
			T3	T4
		T2	X 	T5
		T1	T0(T6)
	*/
	const int offsets[7] = {
		-(int)r2, -(int)r2 - 1, -1, (int)r2, (int)r2+1, 1, -(int)r2
	};
	const T x[6][3] = {
		{0, 0, 1},
		{0, 1, 1},
		{0, 1, 0},
		{1, 1, 0},
		{1, 0, 0},
		{1, 0, 1}
	};
	const T y[6][3] = {
		{1, 0, 1},
		{0, 0, 1},
		{0, 1, 1},
		{0, 1, 0},
		{1, 1, 0},
		{1, 0, 0}
	};
	for(int i=0; i<r1; i++){
		// printf("start %d row\n", i);
		T * cur_U_pos = U_pos;
		T * cur_V_pos = V_pos;
		for(int j=0; j<r2; j++){
			double required_eb = max_pwr_eb;
			// derive eb given six adjacent triangles
			for(int k=0; k<6; k++){
				if(inbound(i*(int)r2 + j + offsets[k], 0, (int)num_elements) && inbound(i*(int)r2 + j + offsets[k+1], 0, (int)num_elements)){
					required_eb = MIN(required_eb, derive_cp_eb_for_positions_online(cur_U_pos[offsets[k]], cur_U_pos[offsets[k+1]], cur_U_pos[0],
						cur_V_pos[offsets[k]], cur_V_pos[offsets[k+1]], cur_V_pos[0], x[k][0], x[k][1], x[k][2], y[k][0], y[k][1], y[k][2]));
				}
			}
			if(required_eb > 0){
				bool unpred_flag = false;
				T decompressed[2];
				// compress U and V
				for(int k=0; k<2; k++){
					T * cur_data_pos = (k == 0) ? cur_U_pos : cur_V_pos;
					T cur_data = *cur_data_pos;
					double abs_eb = fabs(cur_data) * required_eb;
					eb_quant_index_pos[k] = eb_exponential_quantize(abs_eb, base, log_of_base);
					// get adjacent data and perform Lorenzo
					/*
						d2 X
						d0 d1
					*/
					T d0 = (i && j) ? cur_data_pos[-1 - r2] : 0;
					T d1 = (i) ? cur_data_pos[-r2] : 0;
					T d2 = (j) ? cur_data_pos[-1] : 0;
					T pred = d1 + d2 - d0;
					double diff = cur_data - pred;
					double quant_diff = fabs(diff) / abs_eb + 1;
					if(quant_diff < capacity){
						quant_diff = (diff > 0) ? quant_diff : -quant_diff;
						int quant_index = (int)(quant_diff/2) + intv_radius;
						data_quant_index_pos[k] = quant_index;
						decompressed[k] = pred + 2 * (quant_index - intv_radius) * abs_eb; 
						// check original data
						if(fabs(decompressed[k] - cur_data) >= abs_eb){
							unpred_flag = true;
							break;
						}
					}
					else{
						unpred_flag = true;
						break;
					}
				}
				if(unpred_flag){
					// recover quant index
					*(eb_quant_index_pos ++) = 0;
					*(eb_quant_index_pos ++) = 0;
					*(data_quant_index_pos ++) = intv_radius;
					*(data_quant_index_pos ++) = intv_radius;
					unpred_data.push_back(*cur_U_pos);
					unpred_data.push_back(*cur_V_pos);
				}
				else{
					eb_quant_index_pos += 2;
					data_quant_index_pos += 2;
					// assign decompressed data
					*cur_U_pos = decompressed[0];
					*cur_V_pos = decompressed[1];
				}
			}
			else{
				// record as unpredictable data
				*(eb_quant_index_pos ++) = 0;
				*(eb_quant_index_pos ++) = 0;
				*(data_quant_index_pos ++) = intv_radius;
				*(data_quant_index_pos ++) = intv_radius;
				unpred_data.push_back(*cur_U_pos);
				unpred_data.push_back(*cur_V_pos);
			}
			cur_U_pos ++, cur_V_pos ++;
		}
		U_pos += r2;
		V_pos += r2;
	}
	free(decompressed_U);
	free(decompressed_V);
	printf("offsets eb_q, data_q, unpred: %ld %ld %ld\n", eb_quant_index_pos - eb_quant_index, data_quant_index_pos - data_quant_index, unpred_data.size());
	unsigned char * compressed = (unsigned char *) malloc(2*num_elements*sizeof(T));
	unsigned char * compressed_pos = compressed;
	write_variable_to_dst(compressed_pos, base);
	write_variable_to_dst(compressed_pos, intv_radius);
	size_t unpredictable_count = unpred_data.size();
	write_variable_to_dst(compressed_pos, unpredictable_count);
	write_array_to_dst(compressed_pos, (T *)&unpred_data[0], unpredictable_count);	
	Huffman_encode_tree_and_data(2*256, eb_quant_index, 2*num_elements, compressed_pos);
	free(eb_quant_index);
	Huffman_encode_tree_and_data(2*capacity, data_quant_index, 2*num_elements, compressed_pos);
	free(data_quant_index);
	compressed_size = compressed_pos - compressed;
	return compressed;	
}

template
unsigned char *
sz_compress_cp_preserve_2d_online(const float * U, const float * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

template
unsigned char *
sz_compress_cp_preserve_2d_online(const double * U, const double * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

template<typename T>
unsigned char *
sz_compress_cp_preserve_2d_online_log(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb){

	size_t num_elements = r1 * r2;
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

	T * decompressed_U = (T *) malloc(num_elements*sizeof(T));
	memcpy(decompressed_U, U, num_elements*sizeof(T));
	T * decompressed_V = (T *) malloc(num_elements*sizeof(T));
	memcpy(decompressed_V, V, num_elements*sizeof(T));

	int * eb_quant_index = (int *) malloc(num_elements*sizeof(int));
	int * data_quant_index = (int *) malloc(2*num_elements*sizeof(int));
	int * eb_quant_index_pos = eb_quant_index;
	int * data_quant_index_pos = data_quant_index;
	// next, row by row
	const int base = 2;
	const double log_of_base = log2(base);
	const int capacity = 65536;
	const int intv_radius = (capacity >> 1);
	unpred_vec<T> unpred_data;
	// offsets to get six adjacent triangle indices
	// the 7-th rolls back to T0
	/*
	|		T3	T4
	y	T2	X 	T5
	|	T1	T0(T6)
		-	x 	-
	*/
	const int offsets[7] = {
		-(int)r2, -(int)r2 - 1, -1, (int)r2, (int)r2+1, 1, -(int)r2
	};
	const T x[6][3] = {
		{0, 0, 1},
		{0, 1, 1},
		{0, 1, 0},
		{1, 1, 0},
		{1, 0, 0},
		{1, 0, 1}
	};
	const T y[6][3] = {
		{1, 0, 1},
		{0, 0, 1},
		{0, 1, 1},
		{0, 1, 0},
		{1, 1, 0},
		{1, 0, 0}
	};
	T * cur_log_U_pos = log_U;
	T * cur_log_V_pos = log_V;
	T * cur_U_pos = decompressed_U;
	T * cur_V_pos = decompressed_V;
	for(int i=0; i<r1; i++){
		// printf("start %d row\n", i);
		for(int j=0; j<r2; j++){
			double required_eb = max_pwr_eb;
			// derive eb given six adjacent triangles
			for(int k=0; k<6; k++){
				if(inbound(i*(int)r2 + j + offsets[k], 0, (int)num_elements) && inbound(i*(int)r2 + j + offsets[k+1], 0, (int)num_elements)){
					required_eb = MIN(required_eb, derive_cp_eb_for_positions_online(cur_U_pos[offsets[k]], cur_U_pos[offsets[k+1]], cur_U_pos[0],
						cur_V_pos[offsets[k]], cur_V_pos[offsets[k+1]], cur_V_pos[0], x[k][0], x[k][1], x[k][2], y[k][0], y[k][1], y[k][2]));
				}
			}
			if((required_eb > 0) && (*cur_U_pos != 0) && (*cur_V_pos != 0)){
				bool unpred_flag = false;
				T decompressed[2];
				double abs_eb = log(1 + required_eb);
				*eb_quant_index_pos = eb_exponential_quantize(abs_eb, base, log_of_base);
				// compress U and V
				for(int k=0; k<2; k++){
					T * cur_data_pos = (k == 0) ? cur_log_U_pos : cur_log_V_pos;
					T cur_data = *cur_data_pos;
					// get adjacent data and perform Lorenzo
					/*
						d2 X
						d0 d1
					*/
					T d0 = (i && j) ? cur_data_pos[-1 - r2] : 0;
					T d1 = (i) ? cur_data_pos[-r2] : 0;
					T d2 = (j) ? cur_data_pos[-1] : 0;
					T pred = d1 + d2 - d0;
					double diff = cur_data - pred;
					double quant_diff = fabs(diff) / abs_eb + 1;
					if(quant_diff < capacity){
						quant_diff = (diff > 0) ? quant_diff : -quant_diff;
						int quant_index = (int)(quant_diff/2) + intv_radius;
						data_quant_index_pos[k] = quant_index;
						decompressed[k] = pred + 2 * (quant_index - intv_radius) * abs_eb; 
						// check original data
						if(fabs(decompressed[k] - cur_data) >= abs_eb){
							unpred_flag = true;
							break;
						}
					}
					else{
						unpred_flag = true;
						break;
					}
				}
				if(unpred_flag){
					// recover quant index
					*(eb_quant_index_pos ++) = 0;
					*(data_quant_index_pos ++) = intv_radius;
					*(data_quant_index_pos ++) = intv_radius;
					unpred_data.push_back(*cur_U_pos);
					unpred_data.push_back(*cur_V_pos);
				}
				else{
					eb_quant_index_pos ++;
					data_quant_index_pos += 2;
					// assign decompressed data
					*cur_log_U_pos = decompressed[0];
					*cur_log_V_pos = decompressed[1];
					*cur_U_pos = (*cur_U_pos > 0) ? exp2(*cur_log_U_pos) : -exp2(*cur_log_U_pos);
					*cur_V_pos = (*cur_V_pos > 0) ? exp2(*cur_log_V_pos) : -exp2(*cur_log_V_pos);
				}
			}
			else{
				// record as unpredictable data
				*(eb_quant_index_pos ++) = 0;
				*(data_quant_index_pos ++) = intv_radius;
				*(data_quant_index_pos ++) = intv_radius;
				unpred_data.push_back(*cur_U_pos);
				unpred_data.push_back(*cur_V_pos);
			}
			cur_log_U_pos ++, cur_log_V_pos ++;
			cur_U_pos ++, cur_V_pos ++;
		}
	}
	free(log_U);
	free(log_V);
	free(decompressed_U);
	free(decompressed_V);
	printf("offsets eb_q, data_q, unpred: %ld %ld %ld\n", eb_quant_index_pos - eb_quant_index, data_quant_index_pos - data_quant_index, unpred_data.size());
	unsigned char * compressed = (unsigned char *) malloc(2*num_elements*sizeof(T));
	unsigned char * compressed_pos = compressed;
	write_variable_to_dst(compressed_pos, base);
	write_variable_to_dst(compressed_pos, intv_radius);
	write_array_to_dst(compressed_pos, sign_map_compressed, 2*sign_map_size);
	free(sign_map_compressed);
	size_t unpredictable_count = unpred_data.size();
	write_variable_to_dst(compressed_pos, unpredictable_count);
	write_array_to_dst(compressed_pos, (T *)&unpred_data[0], unpredictable_count);	
	Huffman_encode_tree_and_data(2*256, eb_quant_index, num_elements, compressed_pos);
	free(eb_quant_index);
	Huffman_encode_tree_and_data(2*capacity, data_quant_index, 2*num_elements, compressed_pos);
	free(data_quant_index);
	compressed_size = compressed_pos - compressed;
	return compressed;	
}

template
unsigned char *
sz_compress_cp_preserve_2d_online_log(const float * U, const float * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

template
unsigned char *
sz_compress_cp_preserve_2d_online_log(const double * U, const double * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);
