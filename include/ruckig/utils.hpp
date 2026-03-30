#pragma once
// TwinCAT RTOS Stripped Version — All stdlib I/O removed

#include <array>
#include <tuple>
#include <type_traits>

namespace ruckig {

//! Constant for indicating a dynamic (run-time settable) number of DoFs
constexpr static size_t DynamicDOFs {0};

//! Vector data type: fixed-size uses std::array (no heap allocation)
template<class T, size_t DOFs> using StandardVector = std::array<T, DOFs>;
template<class T, size_t DOFs, size_t SIZE> using StandardSizeVector = std::array<T, SIZE>;

//! Integrate with constant jerk for duration t.
inline std::tuple<double, double, double> integrate(double t, double p0, double v0, double a0, double j) {
    return std::make_tuple(
        p0 + t * (v0 + t * (a0 / 2 + t * j / 6)),
        v0 + t * (a0 + t * j / 2),
        a0 + t * j
    );
}

} // namespace ruckig
