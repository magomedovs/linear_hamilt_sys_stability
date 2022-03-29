#include "isapprox.h"

bool isapprox(double u, double v, double tolerance = 1.0e-5) {
    if (u == 0 || v == 0) {
        return (std::abs(u - v) <= tolerance);
    } else {
        return ( (std::abs(u - v)/std::abs(u) <= tolerance) && (std::abs(u - v)/std::abs(v) <= tolerance) );
    }
}

