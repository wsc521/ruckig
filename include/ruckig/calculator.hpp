#pragma once

#include <ruckig/calculator_target.hpp>
#include <ruckig/input_parameter.hpp>
#include <ruckig/trajectory.hpp>


namespace ruckig {

//! Internal interface for the main calculator - (RTOS Stripped version)
template<size_t DOFs, template<class, size_t> class CustomVector = StandardVector>
class Calculator {
public:
    //! Calculator for state-to-state trajectories
    TargetCalculator<DOFs, CustomVector> target_calculator;

    template<size_t D = DOFs, typename std::enable_if<(D >= 1), int>::type = 0>
    explicit Calculator() { }

    template<size_t D = DOFs, typename std::enable_if<(D == 0), int>::type = 0>
    explicit Calculator(size_t dofs): target_calculator(TargetCalculator<DOFs, CustomVector>(dofs)) { }

    //! Calculate the time-optimal state-to-state trajectory
    template<bool throw_error>
    Result calculate(const InputParameter<DOFs, CustomVector>& input, Trajectory<DOFs, CustomVector>& trajectory, double delta_time, bool& was_interrupted) {
        return target_calculator.template calculate<throw_error>(input, trajectory, delta_time, was_interrupted);
    }

    //! Continue the trajectory calculation
    template<bool throw_error>
    Result continue_calculation(const InputParameter<DOFs, CustomVector>& input, Trajectory<DOFs, CustomVector>& trajectory, double delta_time, bool& was_interrupted) {
        return target_calculator.template continue_calculation<throw_error>(input, trajectory, delta_time, was_interrupted);
    }
};

} // namespace ruckig
