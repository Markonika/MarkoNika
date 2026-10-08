#pragma once

#include "mfs/FluidProperties.hpp"
#include "mfs/PipeGeometry.hpp"

namespace mfs {

// -----------------------------------------------------------------------
// Closure relationships, Tables 1-3 of Bonizzi, Andreussi & Banerjee (2009),
// "Flow regime independent, high resolution multi-field modelling of
// near-horizontal gas-liquid flows in pipelines", Int. J. Multiphase Flow
// 35 (2009) 34-46.
//
// All closures here are explicitly INDEPENDENT of flow regime, by design --
// that is the central point of the paper. None of them may branch on a
// "flow pattern" label; they are functions only of local field values
// (velocities, volume fractions, geometry) and fluid/non-dimensional
// numbers (Reynolds, Eotvos, ...).
//
// Two entries (droplet entrainment rate, and the interfacial "wave
// celerity" feeding the bubble entrainment closure) could not be recovered
// unambiguously from the OCR text extraction of the PDF (minus signs and
// square-root radicals are frequently lost by text extractors on this
// particular typeset). Clearly-flagged, literature-consistent substitutes
// are used in their place (see comments at each function) and can be
// swapped out trivially since every closure is a free function taking
// only physical arguments.
// -----------------------------------------------------------------------

// Non-dimensional numbers, Eq. (17), (20).
double reynoldsNumber(double density, double velocity, double hydraulicDiameter, double viscosity);
double eotvosNumber(double densityLiquid, double densityGas, double particleDiameter, double surfaceTension);

// ---- Table 1 / Table 2: wall and interfacial friction factors ----------

enum class FrictionCorrelation {
    TaitelDukler1976,      // Table 1 (reference / default)
    SpeddingHand1997,      // Table 2, liquid wall only
    MoodyRoughWall,        // Table 2, liquid & gas wall
    AndreussiPersen1987,   // Table 2, interfacial
    AndritsosHanratty1987, // Table 2, interfacial
    CohenHanratty1968,     // Table 2, interfacial, constant
    Wallis1969             // Table 2, interfacial
};

struct WallFrictionInputs {
    double density;
    double velocity;
    double hydraulicDiameter;
    double viscosity;
    double roughness = 0.0; // absolute wall roughness [m], for Moody rough-wall closure
};

double wallFrictionFactor(const WallFrictionInputs& in, FrictionCorrelation correlation = FrictionCorrelation::TaitelDukler1976);

struct InterfacialFrictionInputs {
    double densityGas;
    double densityLiquid;
    double viscosityGas;
    double velocityGas;          // superficial-equivalent centre-of-mass velocity of layer 2, or ug
    double velocityLiquid;
    double hydraulicDiameter1;   // D1, layer-1 hydraulic diameter
    double hydraulicDiameter2;   // D2, layer-2 hydraulic diameter
    double liquidHeight;         // h1
    double pipeDiameter;
    double area1, area2;         // A1, A2
    double areaFraction2;        // e2 = A2 / A
    double dA1dh1;               // dA1/dh1 = local chord width (Si)
    double inclination;          // theta, [rad], from horizontal
    double atmPressure = 101325.0;
    double localPressure = 101325.0;
    double gasWallFriction;      // f_gw, needed as the smooth-interface base value

    // Additive shift to the AndreussiPersen1987 branch's own F0=0.36
    // Kelvin-Helmholtz onset threshold (F > F0 + f0Shift triggers the
    // enhancement instead of F > F0). Zero by default, a true no-op
    // (adding 0.0 is exact in IEEE754) -- see
    // mfs::SolverOptions::closureF0Shift for why this exists and how its
    // value, when non-zero, was fit.
    double f0Shift = 0.0;
};

// Returns f_i. Andreussi & Persen (1987) is the correlation used for the
// F > F0 branch by default; Andritsos & Hanratty (1987) is offered as an
// alternative critical-velocity-based closure. Both fall back to f_i = f_gw
// below their respective thresholds, matching Table 2.
double interfacialFrictionFactor(const InterfacialFrictionInputs& in,
                                  FrictionCorrelation correlation = FrictionCorrelation::AndreussiPersen1987);

// Andreussi & Persen (1987) inviscid Kelvin-Helmholtz parameter F, Eq. (19).
double kelvinHelmholtzParameterF(const InterfacialFrictionInputs& in);

// Wall / interfacial shear stress, Eq. (8), (9). tau = 0.5 * f * rho * |u| * u
// (interfacial stress uses the gas density and the gas-liquid relative
// velocity, per Eq. 9).
double wallShearStress(double frictionFactor, double density, double velocity);
double interfacialShearStress(double frictionFactor, double densityGas, double velocityGas, double velocityLiquid);

// ---- Table 3: entrainment / disengagement / deposition -----------------

struct BubbleEntrainmentInputs {
    double densityGas;
    double interfacialWidth;   // Si
    double waveVelocity;       // u_wave -- see waveCelerity() below
    double velocityLiquid;     // ul
    double area;               // pipe cross-section A
};
// phi_e = rho_g * A * [0.076 * (Si/D) * (u_wave - ul) - 0.15], clamped >= 0.
// Nydal & Andreussi (1991).
double bubbleEntrainmentRate(const BubbleEntrainmentInputs& in, double pipeDiameter);

struct BubbleDisengagementInputs {
    double densityGas;
    double densityLiquid;
    double surfaceTension;
    double interfacialWidth;   // Si
    double continuousLiquidFraction; // el
    double K = 0.28;
};
// phi_de = rho_g * K * 1.18 * [sigma*g*(rho_l-rho_g)/rho_l^2]^0.25 * Si * (1-el)
// (Harmathy-type terminal rise velocity scale). Andreussi et al. (1993a,b).
double bubbleDisengagementRate(const BubbleDisengagementInputs& in);

// Interfacial "large-wave"/slug-front celerity used to drive bubble
// entrainment beneath a growing interfacial disturbance. The paper defines
// u_wave in words ("wave celerity") in Table 3 but the OCR text does not
// carry a closed-form expression for it. This implementation uses the same
// Bendiksen (1984) translational front-velocity correlation the paper
// itself uses later (Eq. 25-28) to validate elongated-bubble/slug-front
// propagation -- a physically consistent choice, since the entity that
// entrains gas beneath it is exactly that advancing large-wave/slug front.
// Swap this out for a direct Kelvin-Helmholtz long-wave celerity if a more
// precise closure is required.
struct WaveCelerityInputs {
    double mixtureVelocity; // u_mix = usL + usG (or e1 u1 + e2 u2)
    double pipeDiameter;
    double inclination;     // theta [rad] from horizontal
};
double waveCelerity(const WaveCelerityInputs& in);

struct DropletEntrainmentInputs {
    double pipeDiameter;
    double densityGas;
    double densityLiquid;
    double velocityGas;
    double viscosityLiquid;
    double ke = 7.7e-8;
};
// Pan & Hanratty (2002) type entrainment law, cast (per the excess
// kinetic-energy-above-threshold form widely used for this closure) as
//   Ue = (pi*D/4) * ke * sqrt(rho_g) * max( sqrt(rho_g*rho_l)*ug^2*D - 100*mu_l, 0 )
// The multiplicative structure and constants (ke, the factor 100) are taken
// directly from the legible parts of Table 3; the exact bracketed
// combination could not be fully recovered from the OCR text and is a
// literature-consistent reconstruction -- see Closures.hpp header comment.
double dropletEntrainmentRate(const DropletEntrainmentInputs& in);

struct DropletDepositionInputs {
    double pipeDiameter;
    double densityLiquid;
    double dispersedLiquidFraction; // ed
    double continuousGasFraction;   // eg
    double depositionVelocity = 0.1; // kd [m/s]
};
// Ud = (D/4) * kd * (ed/eg) * rho_l
double dropletDepositionRate(const DropletDepositionInputs& in);

// ---- Drag coefficients & particle sizes ---------------------------------

// Alipchenkov et al. (2004): CD = 18.5/Red^0.6 for 2<Red<500, else 0.44.
double dropletDragCoefficient(double dropletReynolds);

// Tomiyama et al. (1995): CD = max( 24/Reb*(1+0.15*Reb^0.687), 8/3 * Eo/(Eo+4) )
double bubbleDragCoefficient(double bubbleReynolds, double eotvos);

struct DropletDiameterInputs {
    double frictionVelocityGas; // u*g (shear velocity based on gas-wall friction)
    double densityGas;
    double densityLiquid;
    double surfaceTension;
    double pipeDiameter;
};
// Sarkhi & Hanratty (2002): dd = (1/U*) * [4.848*sigma/(rho_g*U*^2*D) + 0.0038] * D
// (non-dimensional droplet-size correlation in terms of the gas friction
// velocity U* = sqrt(tau_wg/rho_g)).
double dropletDiameter(const DropletDiameterInputs& in);

struct BubbleDiameterInputs {
    double wallFrictionFactorLiquid; // f_lw
    double densityLiquid;
    double velocityLiquid1;          // u1 (layer-1 mixture velocity)
    double continuousLiquidFraction; // el
    double weberCritical = 1.05;     // We'_crit
};
// Andreussi et al. (1999): db = sigma/(0.5*f_lw*rho_l*u1^2) * We'_crit * (1 + 51.7*el^1.5)
double bubbleDiameter(const BubbleDiameterInputs& in, double surfaceTension);

// ---- Interfacial drag force on the dispersed fields, Eq. (10) ----------
// F_drag = (3/4) * (CD/d_p) * rho_continuous * e_dispersed * |u_c - u_p| * (u_c - u_p)
// This is not tabulated explicitly in the paper (it is folded into
// "F_drag" in Eq. 10) but is the standard closure consistent with the
// Reynolds-number and CD definitions the paper does give (Eq. 20, Table 3).
double dragForcePerVolume(double dragCoefficient, double particleDiameter,
                           double continuousDensity, double dispersedFraction,
                           double relativeVelocity);

} // namespace mfs
