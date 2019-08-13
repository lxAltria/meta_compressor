#ifndef _sz_compress_cp_preserve_3d_hpp
#define _sz_compress_cp_preserve_3d_hpp

#include <cstddef>

template<typename T>
unsigned char *
sz_compress_cp_preserve_3d_online_log(const T * U, const T * V, const T * W, size_t r1, size_t r2, size_t r3, size_t& compressed_size, bool transpose=false, double max_pwr_eb=0.1);

#endif