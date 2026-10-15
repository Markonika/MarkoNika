// Copyright (c) Marko Nika. All rights reserved. Proprietary and confidential.
// Frequency-dependent added mass A(omega) and radiation damping B(omega) in the time domain (Cummins equation):
//   (M + A_inf) xi'' + integral_0^t K(t - tau) xi'(tau) dtau + C xi = F,
//   K(t) = (2/pi) integral_0^inf B(omega) cos(omega t) domega  (retardation function),
//   A(omega) = A_inf - (1/omega) integral_0^inf K(t) sin(omega t) dt,  B(omega) = integral_0^inf K(t) cos(omega t) dt   (Ogilvie).
#pragma once
#include <deque>
#include <memory>
#include <utility>
#include <vector>
#include "mooring/rigid_body.hpp"

namespace mooring {

struct HydroTable {
    std::vector<double> omega;       // ascending [rad/s]
    std::vector<Mat6> A, B;          // added mass and radiation damping at each omega
    bool hasAinf{false};
    Mat6 Ainf;                       // infinite-frequency added mass, if given
};

class RadiationModel {
public:
    // Retardation functions on a uniform grid of step hK up to length TK from the piecewise-linear B(omega) of the table (B(0) = 0 is assumed; above
    // the last frequency B follows B_N (omega_N/omega)^2 if 'tail' is set, else it is zero). A_inf is taken from the table or evaluated from the Ogilvie
    // relation at the highest table frequency.
    RadiationModel(const HydroTable& t, double TK, double hK, bool tail = true);
    const Mat6& Ainf() const { return Ainf_; }
    double length() const { return TK_; }
    double step() const { return hK_; }
    // Kernel entry K_ij(s) (linear interpolation, zero for s >= TK) and its value at zero lag.
    double K(int i, int j, double s) const;
    double K0(int i, int j) const { return K(i, j, 0.0); }
    // Tabulated frequency-domain data, linearly interpolated (clamped outside the table).
    Mat6 A(double omega) const;
    Mat6 B(double omega) const;
    // Frequency-domain check: the A(omega), B(omega) implied by the retardation function and A_inf (Ogilvie), entry (i, j).
    double AFromKernel(int i, int j, double omega) const;
    double BFromKernel(int i, int j, double omega) const;
    const std::vector<std::pair<int, int>>& pairs() const { return pairs_; }

private:
    HydroTable t_;
    double TK_, hK_;
    Mat6 Ainf_;
    std::vector<std::pair<int, int>> pairs_;
    std::vector<std::vector<double>> kern_;       // per pair, n+1 samples
    int idx_[36];
};

// Convolution state of the body: the velocity history and the memory force  I(t) = integral K(t - tau) v(tau) dtau  (trapezoid rule on the stored steps).
class RadiationMemory {
public:
    explicit RadiationMemory(std::shared_ptr<const RadiationModel> m) : m_(std::move(m)) {}
    void reset(double t0, const Vec6& v0);
    void push(double t, const Vec6& v);                           // append the accepted state
    // I(tNew) without the contribution of the unknown v(tNew): sum over the stored history plus 0.5 K(dt) v_n dt in the last interval.
    Vec6 known(double tNew) const;
    // Zero-lag coefficient of the unknown velocity: 0.5 dt K(0), as a matrix.
    Mat6 zeroLag(double dt) const;
private:
    std::shared_ptr<const RadiationModel> m_;
    std::deque<double> t_;
    std::deque<Vec6> v_;
};

// Natural frequencies of the body with the frequency-dependent added mass: for every mode omega_k^2 (Mrb + A(omega_k)) x = K x is solved by fixed-point
// iteration on omega_k starting from the A_inf modes (A(omega) interpolated from the table; Mrb = rigid-body mass without the radiation added mass).
BodyModes bodyModesRadiation(const Mat6& Mrb, const RadiationModel& rad, const Mat6& K);

}  // namespace mooring
