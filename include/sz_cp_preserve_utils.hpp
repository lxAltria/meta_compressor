#ifndef _sz_cp_preserve_utils_hpp
#define _sz_cp_preserve_utils_hpp

#include <cstddef>
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

template<typename T>
inline bool in_range(T pos, T n){
	return (pos >= 0) && (pos < n);
}

#endif