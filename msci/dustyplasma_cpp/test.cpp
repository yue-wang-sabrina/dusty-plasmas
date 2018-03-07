#include <iostream>
#include <vector>

void convertVecToArr (std::vector<int> vec, int *arr) {
    for (int i=0; i < vec.size(); i++) {
        arr[i*2] = vec[i];
    }
}


extern "C" void testing (int *output) {
        std::vector<int> v = {1,2,3};

        convertVecToArr(v, output);

}
