#ifndef LINSPACE
#define LINSPACE

#include <array>
#include <vector>
#include <cmath>

//template <size_t SIZE>
//void linspace(double interval_begin, double interval_end, std::array<double, SIZE>& arr);

//template <typename T>
//void linspace(T interval_begin, T interval_end, T arr[], size_t SIZE);

/* Fill the std::array with values from range [interval_begin, interval_end]. */
template <size_t SIZE>
void linspace(double interval_begin, double interval_end, std::array<double, SIZE>& arr) {
    double step = std::abs(interval_end - interval_begin) / static_cast<double>(SIZE - 1);
    for (size_t i = 0; i < SIZE; ++i) {
        arr[i] = interval_begin + i * step;
    }
}
/* Fill the std::vector with values from range [interval_begin, interval_end]. */
template <typename T>
void linspace(T interval_begin, T interval_end, std::vector<T>& arr) {
	size_t SIZE = arr.size();
    T step = std::abs(interval_end - interval_begin) / static_cast<T>(SIZE - 1);
    for (size_t i = 0; i < SIZE; ++i) {
        arr[i] = interval_begin + i * step;
    }
}
/* Fill the array of size SIZE with values from range [interval_begin, interval_end]. */
template <typename T>
void linspace(T interval_begin, T interval_end, T arr[], size_t SIZE) {
    T step = std::abs(interval_end - interval_begin) / static_cast<T>(SIZE - 1);
    for (size_t i = 0; i < SIZE; ++i) {
        arr[i] = interval_begin + i * step;
    }
}

#endif
