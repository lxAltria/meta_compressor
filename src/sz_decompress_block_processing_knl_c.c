#include "sz_decompress_block_processing_knl_c.h"

inline void
block_pred_and_decompress_regression_3d_knl(const float * reg_params_pos, float precision, int intv_radius, 
	int size_x, int size_y, int size_z, size_t dim0_offset, size_t dim1_offset, 
	const int * type_pos, int * unpred_count_buffer, const float * unpred_data_buffer, const int offset, float * dec_data_pos){
	for(int i=0; i<size_x; i++){
		float * cur_data_pos = dec_data_pos + i*dim0_offset;
		for(int j=0; j<size_y; j++){
			for(int k=0; k<size_z; k++){
				float pred = reg_params_pos[0] * (float)i + reg_params_pos[1] * (float)j + reg_params_pos[2] * (float)k + reg_params_pos[3];
				int index = j*size_z + k;
				int type_val = type_pos[index];
				if(type_val == 0){
					cur_data_pos[j*dim1_offset + k] = unpred_data_buffer[index*offset + unpred_count_buffer[index]];
					unpred_count_buffer[index] ++;
				}
				else{
					cur_data_pos[j*dim1_offset + k] = pred + (float)(2 * (type_val - intv_radius)) * precision;
				}
			}
		}
		type_pos += size_y * size_z;
	}
}

inline void
block_pred_and_decompress_regression_3d_with_buffer_knl(const float * reg_params_pos, float * buffer, float precision, int intv_radius, 
	int size_x, int size_y, int size_z, size_t buffer_dim0_offset, size_t buffer_dim1_offset,
	size_t dim0_offset, size_t dim1_offset, const int * type_pos, int * unpred_count_buffer, const float * unpred_data_buffer, const int offset, float * dec_data_pos){
	for(int i=0; i<size_x; i++){
		float * cur_data_pos = dec_data_pos + i*dim0_offset;
		float * buffer_pos = buffer + (i+1)*buffer_dim0_offset + buffer_dim1_offset + 1;
		for(int j=0; j<size_y; j++){
			for(int k=0; k<size_z; k++){
				float pred = reg_params_pos[0] * (float)i + reg_params_pos[1] * (float)j + reg_params_pos[2] * (float)k + reg_params_pos[3];
				int index = j*size_z + k;
				int type_val = type_pos[index];
				if(type_val == 0){
					cur_data_pos[j*dim1_offset + k] = buffer_pos[j*buffer_dim1_offset + k] = unpred_data_buffer[index*offset + unpred_count_buffer[index]];
					unpred_count_buffer[index] ++;
				}
				else{
					cur_data_pos[j*dim1_offset + k] = buffer_pos[j*buffer_dim1_offset + k] = pred + (float)(2 * (type_val - intv_radius)) * precision;
				}
			}
		}
		type_pos += size_y * size_z;
	}
}

// block-independant lorenzo pred & decompress
inline void
block_pred_and_decompress_lorenzo_3d_knl_2d_pred(const bool use_mean, const float mean, float * buffer, float precision, int intv_radius, 
	int size_x, int size_y, int size_z, size_t buffer_dim0_offset, size_t buffer_dim1_offset, size_t dim0_offset, size_t dim1_offset, 
	const int * type_pos, int * unpred_count_buffer, const float * unpred_data_buffer, const int offset, float * dec_data_pos){
	for(int i=0; i<size_x; i++){
		for(int j=0; j<size_y; j++){
			float * cur_data_pos = dec_data_pos + i*dim0_offset + j*dim1_offset;
			float * buffer_pos = buffer + (i+1)*buffer_dim0_offset + (j+1)*buffer_dim1_offset + 1;
			for(int k=0; k<size_z; k++){
				int index = j*size_z + k;
				int type_val = type_pos[index];
				if(use_mean && (type_val == 1)) cur_data_pos[k] = buffer_pos[k] = mean;
				else{
					float * cur_buffer_pos = buffer_pos + k;
					float pred = cur_buffer_pos[-buffer_dim0_offset] + cur_buffer_pos[-buffer_dim1_offset] - cur_buffer_pos[- buffer_dim0_offset - buffer_dim1_offset];
					if(type_val == 0){
						cur_data_pos[k] = *cur_buffer_pos = unpred_data_buffer[index*offset + unpred_count_buffer[index]];
						unpred_count_buffer[index] ++;
					}
					else{				
						cur_data_pos[k] = *cur_buffer_pos = pred + (float)(2 * (type_val - intv_radius)) * precision;
					}
				}
			}
		}
		type_pos += size_y * size_z;
	}
}