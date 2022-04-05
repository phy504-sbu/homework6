#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

struct OrbitState {
    double t;
    double x;
    double y;
    double vx;
    double vy;
};

OrbitState rhs(const OrbitState& state);
void write_history(const std::vector<OrbitState>& history);
std::vector<OrbitState> integrate(const double a, const double tmax, const double dt);

const double GM = 4.0 * M_PI * M_PI;   // G * Mass in AU, year, solar mass units

OrbitState rhs(const OrbitState& state) {

    OrbitState dodt{};
    
    // dx/dt = vx; dy/dt = vy

    dodt.x = state.vx;
    dodt.y = state.vy;
    
    // d(vx)/dt = - GMx/r**3; d(vy)/dt = - GMy/r**3

    double r = std::sqrt(state.x * state.x + state.y * state.y);

    dodt.vx = - GM * state.x / std::pow(r, 3);
    dodt.vy = - GM * state.y / std::pow(r, 3);

    return dodt;

}

void write_history(const std::vector<OrbitState>& history, const std::string& prefix) {

    std::ofstream of(prefix + "orbit.dat");

    for (auto o : history) {
        of << std::setw(12) << o.t << " "
           << std::setw(12) << o.x << " "
           << std::setw(12) << o.y << " "
           << std::setw(12) << o.vx << " "
           << std::setw(12) << o.vy << std::endl;
    }

}

std::vector<OrbitState> integrate(const double a, const double tmax, const double dt) {

    // how the history of the orbit

    std::vector<OrbitState> orbit_history{};

    // set initial conditions
    OrbitState state{};

    // assume circular orbit on the x-axis, counter-clockwise orbit

    state.t = 0.0;
    state.x = a;
    state.y = 0.0;
    state.vx = 0.0;
    state.vy = std::sqrt(GM / a);

    orbit_history.push_back(state);

    // integration loop
    while (state.t < tmax) {

        // get the derivatives
        auto state_derivs = rhs(state);

        // update the state
        state.t += dt;
        state.x += state_derivs.x * dt;
        state.y += state_derivs.y * dt;
        state.vx += state_derivs.vx * dt;
        state.vy += state_derivs.vy * dt;

        orbit_history.push_back(state);
    }

    return orbit_history;

}

std::vector<OrbitState> integrate_rk2(const double a, const double tmax, const double dt) {

    // how the history of the orbit

    std::vector<OrbitState> orbit_history{};

    // set initial conditions
    OrbitState state{};

    // assume circular orbit on the x-axis, counter-clockwise orbit

    state.t = 0.0;
    state.x = a;
    state.y = 0.0;
    state.vx = 0.0;
    state.vy = std::sqrt(GM / a);

    orbit_history.push_back(state);

    // integration loop
    while (state.t < tmax) {

        // get the derivatives
        auto state_derivs = rhs(state);

        OrbitState state_mid{};

        // predict the state at the midpoint in time

        state_mid.t = state.t + 0.5 * dt;
        state_mid.x = state.x + state_derivs.x * 0.5 * dt;
        state_mid.y = state.y + state_derivs.y * 0.5 * dt;
        state_mid.vx = state.vx + state_derivs.vx * 0.5 * dt;
        state_mid.vy = state.vy + state_derivs.vy * 0.5 * dt;

        // compute the acceleration at the midpoint in time

        state_derivs = rhs(state_mid);

        // do the final update

        state.t += dt;
        state.x += state_derivs.x * dt;
        state.y += state_derivs.y * dt;
        state.vx += state_derivs.vx * dt;
        state.vy += state_derivs.vy * dt;

        orbit_history.push_back(state);
    }

    return orbit_history;

}

int main() {

    double tmax = 1.0;
    double dt = 0.001;
    double a = 1.0;      // 1 AU

    auto orbit_history_euler = integrate(a, tmax, dt);
    write_history(orbit_history_euler, "euler_");

    auto orbit_history_rk2 = integrate_rk2(a, tmax, dt);
    write_history(orbit_history_rk2, "rk2_");    

}