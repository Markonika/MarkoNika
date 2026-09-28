#pragma once

#include <cassert>
#include <vector>

namespace mfs {

// Thomas algorithm for a tridiagonal system A*x = d, with
// A(i, i-1) = a[i] (a[0] unused), A(i,i) = b[i], A(i,i+1) = c[i] (c[N-1] unused).
// Solves in place; returns x.
inline std::vector<double> solveTridiagonal(std::vector<double> a,
                                             std::vector<double> b,
                                             std::vector<double> c,
                                             std::vector<double> d) {
    const std::size_t n = b.size();
    assert(a.size() == n && c.size() == n && d.size() == n);
    if (n == 0) return {};

    for (std::size_t i = 1; i < n; ++i) {
        const double w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }

    std::vector<double> x(n);
    x[n - 1] = d[n - 1] / b[n - 1];
    for (std::size_t ii = n - 1; ii-- > 0;) {
        x[ii] = (d[ii] - c[ii] * x[ii + 1]) / b[ii];
    }
    return x;
}

} // namespace mfs
