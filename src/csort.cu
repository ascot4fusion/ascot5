#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>

extern "C" {
        void sort_by_key_int_wrapper(int *data, int *value, int N)
        {
                thrust::device_ptr<int> dev_ptr_data(data);
                thrust::device_ptr<int> dev_ptr_value(value);
                thrust::sort_by_key(thrust::device,dev_ptr_data, dev_ptr_data + N, dev_ptr_value);
//                thrust::sort_by_key(thrust::host,data, data + N, value);
        }
}

extern "C" {
        void sort_int_wrapper(int *data, int N)
        {
                thrust::device_ptr<int> dev_ptr_data(data);
                thrust::sort(thrust::device,dev_ptr_data, dev_ptr_data + N);
        }
}