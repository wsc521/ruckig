#pragma once

#include <array>
#include <limits>
#include <optional>
#include <type_traits>

#include <ruckig/error.hpp>
#include <ruckig/result.hpp>
#include <ruckig/utils.hpp>

namespace ruckig {

enum class ControlInterface {
    Position, ///< Position-control: Full control over the entire kinematic state (Default)
    Velocity, ///< Velocity-control: Ignores the current position, target position, and velocity limits
};

enum class Synchronization {
    Time, ///< Always synchronize the DoFs to reach the target at the same time (Default)
    TimeIfNecessary, ///< Synchronize only when necessary (e.g. for non-zero target velocity or acceleration)
    Phase, ///< Phase synchronize the DoFs when this is possible, else fallback to "Time" strategy. Phase synchronization will result a straight-line trajectory
    None, ///< Calculate every DoF independently
};

enum class DurationDiscretization {
    Continuous, ///< Every trajectory synchronization duration is allowed (Default)
    Discrete, ///< The trajectory synchronization duration must be a multiple of the control cycle
};


//! Input of the Ruckig algorithm - (RTOS Stripped version)
template<size_t DOFs, template<class, size_t> class CustomVector = StandardVector>
class InputParameter {
    template<class T> using Vector = CustomVector<T, DOFs>;

    inline static double v_at_a_zero(double v0, double a0, double j) {
        return v0 + (a0 * a0) / (2 * j);
    }

    void initialize() {
        for (size_t dof = 0; dof < degrees_of_freedom; ++dof) {
            current_velocity[dof] = 0.0;
            current_acceleration[dof] = 0.0;
            target_velocity[dof] = 0.0;
            target_acceleration[dof] = 0.0;
            max_acceleration[dof] = std::numeric_limits<double>::infinity();
            max_jerk[dof] = std::numeric_limits<double>::infinity();
            enabled[dof] = true;
        }
    }

    void resize(size_t dofs) {
        // resize is usually for dynamic sizes, in TwinCAT we use fixed DOFs
        // But we keep it for template compatibility
    }

public:
    size_t degrees_of_freedom;

    ControlInterface control_interface {ControlInterface::Position};
    Synchronization synchronization {Synchronization::Time};
    DurationDiscretization duration_discretization {DurationDiscretization::Continuous};

    //! Current state
    Vector<double> current_position, current_velocity, current_acceleration;

    //! Target state
    Vector<double> target_position, target_velocity, target_acceleration;

    //! Kinematic constraints
    Vector<double> max_velocity, max_acceleration, max_jerk;
    std::optional<Vector<double>> min_velocity, min_acceleration;

    //! Waypoints/Intermediate positions disabled for TwinCAT RTOS
    // std::vector<Vector<double>> intermediate_positions; 

    //! Is the DoF considered for calculation?
    Vector<bool> enabled;

    //! Per-DoF control_interface (overwrites global control_interface)
    std::optional<Vector<ControlInterface>> per_dof_control_interface;

    //! Per-DoF synchronization (overwrites global synchronization)
    std::optional<Vector<Synchronization>> per_dof_synchronization;

    //! Optional minimum trajectory duration
    std::optional<double> minimum_duration;

    template<size_t D = DOFs, typename std::enable_if<(D >= 1), int>::type = 0>
    InputParameter(): degrees_of_freedom(DOFs) {
        initialize();
    }

    template<size_t D = DOFs, typename std::enable_if<(D == 0), int>::type = 0>
    InputParameter(size_t dofs): degrees_of_freedom(dofs) {
        // resize(dofs);
        initialize();
    }

    //! Validate the input for trajectory calculation
    template<bool throw_validation_error = true>
    bool validate(bool check_current_state_within_limits = false, bool check_target_state_within_limits = true) const {
        for (size_t dof = 0; dof < degrees_of_freedom; ++dof) {
            const double jMax = max_jerk[dof];
            if (std::isnan(jMax) || jMax < 0.0) {
                return false;
            }

            const double aMax = max_acceleration[dof];
            if (std::isnan(aMax) || aMax < 0.0) {
                return false;
            }

            const double aMin = min_acceleration ? min_acceleration.value()[dof] : -max_acceleration[dof];
            if (std::isnan(aMin) || aMin > 0.0) {
                return false;
            }

            const double a0 = current_acceleration[dof];
            if (std::isnan(a0)) { return false; }
            const double af = target_acceleration[dof];
            if (std::isnan(af)) { return false; }

            if (check_current_state_within_limits) {
                if (a0 > aMax || a0 < aMin) return false;
            }
            if (check_target_state_within_limits) {
                if (af > aMax || af < aMin) return false;
            }

            const double v0 = current_velocity[dof];
            if (std::isnan(v0)) { return false; }
            const double vf = target_velocity[dof];
            if (std::isnan(vf)) { return false; }

            auto control_interface_ = per_dof_control_interface ? per_dof_control_interface.value()[dof] : control_interface;
            if (control_interface_ == ControlInterface::Position) {
                const double p0 = current_position[dof];
                if (std::isnan(p0)) { return false; }
                const double pf = target_position[dof];
                if (std::isnan(pf)) { return false; }

                const double vMax = max_velocity[dof];
                if (std::isnan(vMax) || vMax < 0.0) return false;

                const double vMin = min_velocity ? min_velocity.value()[dof] : -max_velocity[dof];
                if (std::isnan(vMin) || vMin > 0.0) return false;

                if (check_current_state_within_limits) {
                    if (v0 > vMax || v0 < vMin) return false;
                }
                if (check_target_state_within_limits) {
                    if (vf > vMax || vf < vMin) return false;
                }
            }
        }
        return true;
    }

    bool operator!=(const InputParameter<DOFs, CustomVector>& rhs) const {
        return !(
            current_position == rhs.current_position
            && current_velocity == rhs.current_velocity
            && current_acceleration == rhs.current_acceleration
            && target_position == rhs.target_position
            && target_velocity == rhs.target_velocity
            && target_acceleration == rhs.target_acceleration
            && max_velocity == rhs.max_velocity
            && max_acceleration == rhs.max_acceleration
            && max_jerk == rhs.max_jerk
            && enabled == rhs.enabled
            && minimum_duration == rhs.minimum_duration
            && min_velocity == rhs.min_velocity
            && min_acceleration == rhs.min_acceleration
            && control_interface == rhs.control_interface
            && synchronization == rhs.synchronization
            && duration_discretization == rhs.duration_discretization
            && per_dof_control_interface == rhs.per_dof_control_interface
            && per_dof_synchronization == rhs.per_dof_synchronization
        );
    }
};

} // namespace ruckig
