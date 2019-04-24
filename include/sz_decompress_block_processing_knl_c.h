#ifndef _sz_decompress_block_processing_knl_c_h
#define _sz_decompress_block_processing_knl_c_h
#ifdef __cplusplus
extern "C" {
#endif
#include <stddef.h>
#include <math.h>
#include <stdbool.h>

void
block_pred_and_decompress_regression_3d_knl(const float * reg_params_pos, float precision, int intv_radius, 
	int size_x, int size_y, int size_z, size_t dim0_offset, size_t dim1_offset, 
	const int * type_pos, int * unpred_count_buffer, const float * unpred_data_buffer, const int offset, float * dec_data_pos);

void
block_pred_and_decompress_regression_3d_with_buffer_knl_c(const float * reg_params_pos, float * buffer, float precision, int intv_radius, 
	int size_x, int size_y, int size_z, size_t buffer_dim0_offset, size_t buffer_dim1_offset,
	size_t dim0_offset, size_t dim1_offset, const int * type_pos, int * unpred_count_buffer, const float * unpred_data_buffer, const int offset, float * dec_data_pos);

void
block_pred_and_decompress_lorenzo_3d_knl_2d_pred_c(const bool use_mean, const float mean, float * buffer, float precision, int intv_radius, 
	int size_x, int size_y, int size_z, size_t buffer_dim0_offset, size_t buffer_dim1_offset, size_t dim0_offset, size_t dim1_offset, 
	const int * type_pos, int * unpred_count_buffer, const float * unpred_data_buffer, const int offset, float * dec_data_pos);

#ifdef __cplusplus
}
#endif
#endif
