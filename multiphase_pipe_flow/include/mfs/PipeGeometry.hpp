#pragma once

namespace mfs {

// Circular pipe cross-section geometry for the stratified "two layer" picture
// used by the four-field model (Bonizzi, Andreussi & Banerjee, 2009).
//
// Layer 1 (continuous liquid + dispersed gas bubbles) is assumed, purely for
// the purpose of computing wetted perimeters / interfacial width / hydraulic
// diameters, to occupy the bottom of the pipe as a circular segment of area
// A1, exactly as in the classical Taitel & Dukler (1976) stratified-flow
// geometry. Layer 2 (continuous gas + dispersed liquid droplets) occupies the
// remainder of the cross-section, A2 = A - A1.
//
// This geometric picture is unambiguous and well validated for horizontal
// and near-horizontal pipes. As the pipe inclination approaches vertical the
// notion of a gravity-segregated "bottom layer" loses physical meaning
// (there is no preferred low side of the cross-section), but the governing
// momentum equations remain well defined because the hydrostatic
// "level-gradient" term is multiplied by cos(theta) and therefore vanishes
// smoothly as theta -> 90 deg (see FourFieldSolver.hpp). The perimeters
// computed here are then best interpreted as an effective, area-based
// wetted/interfacial-perimeter closure rather than a literal geometric
// picture -- this is the same pragmatic extension used by most
// industrial-grade transient multiphase pipe-flow codes when applied outside
// the strictly horizontal regime, and is flagged explicitly in the README.
class PipeGeometry {
public:
    explicit PipeGeometry(double diameter);

    double diameter() const { return diameter_; }
    double radius() const { return 0.5 * diameter_; }
    double area() const { return area_; }

    struct StratifiedGeometry {
        double A1;      // cross-sectional area occupied by layer 1 [m^2]
        double A2;      // cross-sectional area occupied by layer 2 [m^2]
        double Swp1;    // wall perimeter wetted by layer 1 [m]
        double Swp2;    // wall perimeter wetted by layer 2 [m]
        double Si;      // interfacial width (chord length) [m]
        double h1;       // "liquid" layer height measured from the bottom [m]
        double D1;      // hydraulic diameter of layer 1, 4*A1/Swp1 [m]
        double D2;      // hydraulic diameter of layer 2, 4*A2/(Swp2+Si) [m]
    };

    // Build the full stratified geometry from the layer-1 area fraction
    // e1 = A1 / A, e1 in (0, 1).
    StratifiedGeometry fromAreaFraction(double e1) const;

    // dA1/dh1 at a given height, i.e. the local chord width. Identically
    // equal to Si(h1); provided separately since it is needed to convert a
    // gradient of e1 along the pipe axis into a gradient of liquid height
    // (the "hydraulic head" term in the layer-1/layer-2 momentum equations).
    double chordWidth(double h1) const;

    // Invert area -> height (monotinic, safeguarded Newton/bisection).
    double heightFromArea(double A1) const;

private:
    double diameter_;
    double area_;
};

} // namespace mfs
