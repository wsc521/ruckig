#pragma once

// Ruckig Cloud Calculator - Stubbed out for TwinCAT RTOS
// This feature is not supported in the stripped version

namespace ruckig {

template<size_t DOFs, template<class, size_t> class CustomVector = StandardVector>
class WaypointsCalculator {
public:
    explicit WaypointsCalculator() { }
    explicit WaypointsCalculator(size_t) { }

    template<bool throw_error>
    Result calculate(const InputParameter<DOFs, CustomVector>&, Trajectory<DOFs, CustomVector>&, double, bool&) {
        return Result::Error;
    }

    template<bool throw_error>
    Result continue_calculation(const InputParameter<DOFs, CustomVector>&, Trajectory<DOFs, CustomVector>&, double, bool&) {
        return Result::Error;
    }
};

} // namespace ruckig
