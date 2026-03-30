#pragma once

#include <array>
#include <type_traits>

#include <ruckig/trajectory.hpp>
#include <ruckig/utils.hpp>


namespace ruckig {

//! Output of the Ruckig algorithm - (RTOS Stripped version)
template<size_t DOFs, template<class, size_t> class CustomVector = StandardVector>
class OutputParameter {
    template<class T> using Vector = CustomVector<T, DOFs>;

public:
    size_t degrees_of_freedom;

    //! Current trajectory
    Trajectory<DOFs, CustomVector> trajectory;

    // Current kinematic state
    Vector<double> new_position, new_velocity, new_acceleration, new_jerk;

    //! Current time on trajectory
    double time {0.0};

    //! Index of the current section (Stays 0 in single-segment mode)
    size_t new_section {0};

    //! Was a new section reached?
    bool did_section_change {false};

    //! Was a new trajectory calculation performed?
    bool new_calculation {false};

    //! Was the trajectory calculation interrupted? (Legacy)
    bool was_calculation_interrupted {false};

    //! Computational duration
    double calculation_duration {0.0};

    template<size_t D = DOFs, typename std::enable_if<(D >= 1), int>::type = 0>
    OutputParameter(): degrees_of_freedom(DOFs) { }

    template<size_t D = DOFs, typename std::enable_if<(D == 0), int>::type = 0>
    OutputParameter(size_t dofs):
        degrees_of_freedom(dofs),
        trajectory(Trajectory<0, CustomVector>(dofs))
    {
    }

    void pass_to_input(InputParameter<DOFs, CustomVector>& input) const {
        input.current_position = new_position;
        input.current_velocity = new_velocity;
        input.current_acceleration = new_acceleration;
    }
};

} // namespace ruckig
