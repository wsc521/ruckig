#pragma once

#include <algorithm>
#include <array>
#include <limits>
#include <math.h>
#include <optional>
#include <tuple>

#include <ruckig/calculator.hpp>
#include <ruckig/error.hpp>
#include <ruckig/input_parameter.hpp>
#include <ruckig/output_parameter.hpp>
#include <ruckig/trajectory.hpp>


namespace ruckig {

//! Main interface for the Ruckig algorithm - (RTOS Stripped version)
template<size_t DOFs = 0, template<class, size_t> class CustomVector = StandardVector, bool throw_error = false>
class Ruckig {
    //! Current input, only for comparison for recalculation
    InputParameter<DOFs, CustomVector> current_input;

    //! Flag that indicates if the current_input was properly initialized
    bool current_input_initialized {false};

public:
    //! Calculator for new trajectories
    Calculator<DOFs, CustomVector> calculator;

    //! Degrees of freedom
    const size_t degrees_of_freedom;

    //! Time step between updates (cycle time) in [s]
    double delta_time {0.0};

    template<size_t D = DOFs, typename std::enable_if<(D >= 1), int>::type = 0>
    explicit Ruckig():
        degrees_of_freedom(DOFs),
        delta_time(-1.0)
    {
    }

    template<size_t D = DOFs, typename std::enable_if<(D >= 1), int>::type = 0>
    explicit Ruckig(double delta_time):
        degrees_of_freedom(DOFs),
        delta_time(delta_time)
    {
    }

    template<size_t D = DOFs, typename std::enable_if<(D == 0), int>::type = 0>
    explicit Ruckig(size_t dofs):
        current_input(InputParameter<DOFs, CustomVector>(dofs)),
        calculator(Calculator<DOFs, CustomVector>(dofs)),
        degrees_of_freedom(dofs),
        delta_time(-1.0)
    {
    }

    template<size_t D = DOFs, typename std::enable_if<(D == 0), int>::type = 0>
    explicit Ruckig(size_t dofs, double delta_time):
        current_input(InputParameter<DOFs, CustomVector>(dofs)),
        calculator(Calculator<DOFs, CustomVector>(dofs)),
        degrees_of_freedom(dofs),
        delta_time(delta_time)
    {
    }

    //! Reset the instance (e.g. to force a new calculation in the next update)
    void reset() {
        current_input_initialized = false;
    }

    //! Validate the input as well as the Ruckig instance for trajectory calculation
    template<bool throw_validation_error = true>
    bool validate_input(const InputParameter<DOFs, CustomVector>& input, bool check_current_state_within_limits = false, bool check_target_state_within_limits = true) const {
        if (!input.template validate<throw_validation_error>(check_current_state_within_limits, check_target_state_within_limits)) {
            return false;
        }

        if (delta_time <= 0.0 && input.duration_discretization != DurationDiscretization::Continuous) {
            return false;
        }

        return true;
    }

    //! Calculate a new trajectory for the given input
    Result calculate(const InputParameter<DOFs, CustomVector>& input, Trajectory<DOFs, CustomVector>& trajectory) {
        bool was_interrupted {false};
        return calculate(input, trajectory, was_interrupted);
    }

    //! Calculate a new trajectory for the given input and check for interruption
    Result calculate(const InputParameter<DOFs, CustomVector>& input, Trajectory<DOFs, CustomVector>& trajectory, bool& was_interrupted) {
        if (!validate_input<throw_error>(input, false, true)) {
            return Result::ErrorInvalidInput;
        }

        return calculator.template calculate<throw_error>(input, trajectory, delta_time, was_interrupted);
    }

    //! Get the next output state (with step delta_time) along the calculated trajectory for the given input
    Result update(const InputParameter<DOFs, CustomVector>& input, OutputParameter<DOFs, CustomVector>& output) {
        output.new_calculation = false;

        Result result {Result::Working};
        if (!current_input_initialized || input != current_input) {
            result = calculate(input, output.trajectory, output.was_calculation_interrupted);
            if (result != Result::Working && result != Result::ErrorPositionalLimits) {
                return result;
            }

            current_input = input;
            current_input_initialized = true;
            output.time = 0.0;
            output.new_calculation = true;
        }

        const size_t old_section = output.new_section;
        output.time += delta_time;
        output.trajectory.at_time(output.time, output.new_position, output.new_velocity, output.new_acceleration, output.new_jerk, output.new_section);
        output.did_section_change = (output.new_section > old_section); 

        output.calculation_duration = 0.0;
        output.pass_to_input(current_input);

        if (output.time > output.trajectory.get_duration()) {
            return Result::Finished;
        }

        return result;
    }
};

} // namespace ruckig
