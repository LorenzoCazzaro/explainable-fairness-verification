#include "mytype.h"

#include <set>
#include <map>

#include <iostream>
#include <fstream>

#include <algorithm>
#include <numeric>

#include <chrono>
using namespace std::chrono;


#define SEED 17
using namespace std;


template <typename T>
T min_vect(std::vector<T>& v) {
    return *std::min_element(std::begin(v), std::end(v));
}

template <typename T>
T max_vect(std::vector<T>& v) {
    return *std::max_element(std::begin(v), std::end(v));
}

template <typename T>
T rand_val(T end, T start = 0) {
    return start + static_cast <T> (rand()) / ( static_cast <T> (RAND_MAX / end - start));
}

/* arg sort */
std::vector<size_t> arg_sort(std::vector<size_t>& v) {
    std::vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];} );
    return idx;
}

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

/* load JSON file */
std::string read_file(std::string file_name, std::string source) {
    std::ifstream ifile;
    ifile.open(file_name);
    if(!ifile) {
        std::cout << source << " path/file <" << file_name << "> doesn't exist" << std::endl;
        exit(0);
    }

    std::string s((std::istreambuf_iterator<char>(ifile)), std::istreambuf_iterator<char>());
    return s;    
}

/* max value in a map */
size_t max_map(std::map<float, size_t>& map, float* label=nullptr, bool flag_print=false) {
    std::map<float, size_t>::iterator itr = map.begin();
    if (label != nullptr)
        *label = itr->first;
    size_t max = itr->second;

    for (itr = map.begin(); itr != map.end(); ++itr) {
        if (flag_print)
            std::cout << "key: " << itr->first << " value: " << itr->second << std::endl;
        if (max < itr->second) {
            max = itr->second;
            if (label != nullptr)
                *label = itr->first;
        }
    }
    if (flag_print)
        std::cout << "max: " << max << std::endl;
    return max;
}

/* major voted value in with a map */
float majority_voting(std::map<float, size_t>& map) {
    std::map<float, size_t>::iterator itr = map.begin();
    size_t max = itr->second;
    float label = itr->first;

    for (itr = map.begin(); itr != map.end(); ++itr)
        if (max < itr->second) {
            max = itr->second;
            label = itr->first;
        }

    return label;
}

/* print vector of float */
void print_ist(V_f v) {
    for (size_t i = 0; i < v.size(); ++i)
        std::cout << v[i] << " ";

    std::cout << std::endl;
}