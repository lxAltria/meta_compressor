#include "sz_compress_block_processing_knl_c.h"
#include <stddef.h>

// optimizations on knl
void
block_pred_and_quant_regression_3d_knl_c(const float * data_pos, const float * reg_params_pos, float precision, float recip_precision, int capacity, 
	int intv_radius, int size_x, int size_y, int size_z, size_t dim0_offset, size_t dim1_offset, 
	int * type_pos, int * unpred_count_buffer, float * unpred_buffer, size_t offset){
	for(int i=0; i<size_x; i++){
		const float * cur_data_pos = data_pos + i*dim0_offset;
		for(int j=0; j<size_y; j++){
			for(int k=0; k<size_z; k++){
				float cur_data = cur_data_pos[j*dim1_offset + k];
				float pred = reg_params_pos[0] * (float)i + reg_params_pos[1] * (float)j + reg_params_pos[2] * (float)k + reg_params_pos[3];
				float diff = cur_data - pred;
				int quant_index = (int)(fabs(diff) * recip_precision) + 1;
				if(quant_index < capacity){
					quant_index >>= 1;
					int half_index = quant_index;
					quant_index <<= 1;
					if(diff < 0){
						quant_index = - quant_index;
						type_pos[j*size_z + k] = intv_radius - half_index;
					}
					else type_pos[j*size_z + k] = intv_radius + half_index;
					float decompressed_data = pred + (float)quant_index * precision;
					if(fabs(cur_data - decompressed_data) > precision){
						int index = j*size_z + k;
						type_pos[index] = 0;
						unpred_buffer[index*offset + unpred_count_buffer[index]] = cur_data;
						unpred_count_buffer[index] ++;
					}
			 	}
			 	else{
					int index = j*size_z + k;
					type_pos[index] = 0;
					unpred_buffer[index*offset + unpred_count_buffer[index]] = cur_data;
					unpred_count_buffer[index] ++;
			 	}
			}
		}
		type_pos += size_y * size_z;
	}
}

void
block_pred_and_quant_regression_3d_with_buffer_knl_c(const float * data_pos, const float * reg_params_pos, float * buffer, float precision, float recip_precision, int capacity, 
	int intv_radius, int size_x, int size_y, int size_z, size_t buffer_dim0_offset, size_t buffer_dim1_offset,
	size_t dim0_offset, size_t dim1_offset, int * type_pos, int * unpred_count_buffer, float * unpred_buffer, size_t offset){
	for(int i=0; i<size_x; i++){
		const float * cur_data_pos = data_pos + i*dim0_offset;
		float * buffer_pos = buffer + (i+1)*buffer_dim0_offset + buffer_dim1_offset + 1;
		for(int j=0; j<size_y; j++){
			for(int k=0; k<size_z; k++){
				float cur_data = cur_data_pos[j*dim1_offset + k];
				float pred = (float) (reg_params_pos[0] * (float)i + reg_params_pos[1] * (float)j + reg_params_pos[2] * (float)k + reg_params_pos[3]);
				float diff = cur_data - pred;
				int quant_index = (int) (fabs(diff) * recip_precision) + 1;
				if(quant_index < capacity){
					quant_index >>= 1;
					int half_index = quant_index;
					quant_index <<= 1;
					if(diff < 0){
						quant_index = - quant_index;
						type_pos[j*size_z + k] = intv_radius - half_index;
					}
					else type_pos[j*size_z + k] = intv_radius + half_index;
					float decompressed_data = pred + (float)quant_index * precision;
					if(fabs(cur_data - decompressed_data) > precision){
						int index = j*size_z + k;
						type_pos[index] = 0;
						unpred_buffer[index*offset + unpred_count_buffer[index]] = buffer_pos[j*buffer_dim1_offset + k] = cur_data;
						unpred_count_buffer[index] ++;
					}
					else buffer_pos[j*buffer_dim1_offset + k] = decompressed_data;
			 	}
			 	else{
					int index = j*size_z + k;
					type_pos[index] = 0;
					unpred_buffer[index*offset + unpred_count_buffer[index]] = buffer_pos[j*buffer_dim1_offset + k] = cur_data;
					unpred_count_buffer[index] ++;
			 	}
			}
		}
		type_pos += size_y * size_z;
	}
}

void
block_pred_and_quant_lorenzo_3d_knl_2d_pred_c(const bool use_mean, const float mean, const float * data_pos, float * buffer, float precision, float recip_precision, int capacity, int intv_radius, 
	int size_x, int size_y, int size_z, size_t buffer_dim0_offset, size_t buffer_dim1_offset,
	size_t dim0_offset, size_t dim1_offset, int * type_pos, int * unpred_count_buffer, float * unpred_buffer, size_t offset){
	const float * cur_data_pos = data_pos;
	for(int i=0; i<size_x; i++){
		for(int j=0; j<size_y; j++){
			const float * cur_data_pos = data_pos + i*dim0_offset + j*dim1_offset;
			float * buffer_pos = buffer + (i+1)*buffer_dim0_offset + (j+1)*buffer_dim1_offset + 1;
			for(int k=0; k<size_z; k++){
				if(use_mean && (fabs(cur_data_pos[k] - mean) < precision)){
					type_pos[k] = 1;
					buffer_pos[k] = mean;
				}
				else{
					// reduce to 2D prediction
					float * cur_buffer_pos = buffer_pos + k;
					float cur_data = cur_data_pos[k];
					float pred = cur_buffer_pos[-buffer_dim0_offset] + cur_buffer_pos[-buffer_dim1_offset] - cur_buffer_pos[- buffer_dim0_offset - buffer_dim1_offset];
					float diff = cur_data - pred;
					int quant_index = (int)(fabs(diff) * recip_precision) + 1;
					if(quant_index < capacity){
						quant_index >>= 1;
						int half_index = quant_index;
						quant_index <<= 1;
						if(diff < 0){
							quant_index = - quant_index;
							half_index = - half_index;
						}
						type_pos[k] = half_index + intv_radius;
						float decompressed_data = pred + (float)quant_index * precision;
						if(fabs(cur_data - decompressed_data) > precision){
							int index = j*size_z + k;
							unpred_buffer[index*offset + unpred_count_buffer[index]] = *cur_buffer_pos = cur_data;
							unpred_count_buffer[index] ++;
							type_pos[k] = 0;
						}
					 	else *cur_buffer_pos = decompressed_data;
				 	}
				 	else{
						int index = j*size_z + k;
						unpred_buffer[index*offset + unpred_count_buffer[index]] = *cur_buffer_pos = cur_data;
						unpred_count_buffer[index] ++;
						type_pos[k] = 0;
				 	}
				}
			}
			type_pos += size_z;
		}
	}
}

