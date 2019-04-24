#ifndef _sz_compress_block_processing_knl_c_h
#define _sz_compress_block_processing_knl_c_h
#ifdef __cplusplus
extern "C" {
#endif
#include <stddef.h>
#include <math.h>
#include <stdbool.h>

inline void
block_pred_and_quant_regression_3d_knl_c(const float * data_pos, const float * reg_params_pos, float precision, float recip_precision, int capacity, 
	int intv_radius, int size_x, int size_y, int size_z, size_t dim0_offset, size_t dim1_offset, 
	int * type_pos, int * unpred_count_buffer, float * unpred_buffer, size_t offset);

inline void
block_pred_and_quant_regression_3d_with_buffer_knl_c(const float * data_pos, const float * reg_params_pos, float * buffer, float precision, float recip_precision, int capacity, 
	int intv_radius, int size_x, int size_y, int size_z, size_t buffer_dim0_offset, size_t buffer_dim1_offset,
	size_t dim0_offset, size_t dim1_offset, int * type_pos, int * unpred_count_buffer, float * unpred_buffer, size_t offset);

inline void
block_pred_and_quant_lorenzo_3d_knl_2d_pred_c(const bool use_mean, const float mean, const float * data_pos, float * buffer, float precision, float recip_precision, int capacity, int intv_radius, 
	int size_x, int size_y, int size_z, size_t buffer_dim0_offset, size_t buffer_dim1_offset,
	size_t dim0_offset, size_t dim1_offset, int * type_pos, int * unpred_count_buffer, float * unpred_buffer, size_t offset);

#ifdef __cplusplus
}
#endif
#endif
