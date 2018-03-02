/**
 * simple_add is a simple test kernel file
 *
 * author: Tom Nordloh (nordloh@uni-bremen.de)
 */

__kernel void simple_add(__global const int* A,
						  __global const int* B,
						  __global int* C){
	C[get_global_id(0)] = A[get_global_id(0)] + B[get_global_id(0)];
}
