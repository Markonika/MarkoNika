#pragma once

#include <vector>

namespace mfs {

// Staggered-grid (Harlow & Welch, 1965) state for the four-field model.
//
// N control volumes span the pipe, indexed 0..N-1 (cell centres). There are
// N+1 faces, indexed 0..N (face i sits at the boundary between cell i-1 and
// cell i for 1<=i<=N-1; face 0 is the inlet, face N is the outlet).
//
// Cell-centred quantities: volume fractions, pressure, densities, geometry.
// Face-centred quantities: the four field velocities.
struct FlowState {
    int N = 0;
    double L = 0.0;
    double dz = 0.0;

    // Cell-centred, size N.
    std::vector<double> el;   // continuous liquid volume fraction
    std::vector<double> ed;   // dispersed liquid (droplet) volume fraction
    std::vector<double> eg;   // continuous gas volume fraction
    std::vector<double> eb;   // dispersed gas (bubble) volume fraction
    std::vector<double> P;    // pressure [Pa]
    std::vector<double> theta; // pipe inclination from horizontal [rad], can vary with z

    // Face-centred, size N+1.
    std::vector<double> u1;   // layer-1 (liquid continuous + bubbles) centre-of-mass velocity
    std::vector<double> u2;   // layer-2 (gas continuous + droplets) centre-of-mass velocity
    std::vector<double> ud;   // droplet velocity
    std::vector<double> ub;   // bubble velocity
    // Derived (recomputed each step, kept here for diagnostics / output).
    std::vector<double> ul;   // continuous liquid velocity, back-substituted from u1, ub
    std::vector<double> ug;   // continuous gas velocity, back-substituted from u2, ud

    void resize(int nCells) {
        N = nCells;
        el.assign(N, 0.0);
        ed.assign(N, 0.0);
        eg.assign(N, 0.0);
        eb.assign(N, 0.0);
        P.assign(N, 0.0);
        theta.assign(N, 0.0);

        u1.assign(N + 1, 0.0);
        u2.assign(N + 1, 0.0);
        ud.assign(N + 1, 0.0);
        ub.assign(N + 1, 0.0);
        ul.assign(N + 1, 0.0);
        ug.assign(N + 1, 0.0);
    }

    double eL(int i) const { return el[i] + ed[i]; }
    double eG(int i) const { return eg[i] + eb[i]; }
    double e1(int i) const { return el[i] + eb[i]; } // layer-1 area fraction
    double e2(int i) const { return eg[i] + ed[i]; } // layer-2 area fraction

    double cellCenter(int i) const { return (i + 0.5) * dz; }
    double faceCenter(int i) const { return i * dz; }
};

// Cumulative mass-conservation tracker, Eq. (24): sums (Mdot_out - Mdot_in)*dt
// over the run and compares against the field's storage change.
struct MassBalanceTracker {
    double integratedNetOutflow = 0.0; // sum over steps of (out - in) * dt
    double initialMass = 0.0;
    bool initialised = false;

    void accumulate(double massInRate, double massOutRate, double dt) {
        integratedNetOutflow += (massOutRate - massInRate) * dt;
    }

    double error(double currentMass) const {
        if (!initialised) return 0.0;
        return (currentMass - initialMass) - integratedNetOutflow;
    }
};

} // namespace mfs
