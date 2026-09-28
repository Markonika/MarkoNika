#include "mfs/Closures.hpp"
#include "mfs/Constants.hpp"

#include <algorithm>
#include <cmath>

namespace mfs {

namespace {
    inline double sgn(double x) { return (x > 0.0) - (x < 0.0); }
}

double reynoldsNumber(double density, double velocity, double hydraulicDiameter, double viscosity) {
    return density * std::fabs(velocity) * hydraulicDiameter / std::max(viscosity, constants::tiny);
}

double eotvosNumber(double densityLiquid, double densityGas, double particleDiameter, double surfaceTension) {
    return constants::gravity * (densityLiquid - densityGas) * particleDiameter * particleDiameter /
           std::max(surfaceTension, constants::tiny);
}

// ---- Table 1 / Table 2 --------------------------------------------------

namespace {
    // Taitel & Dukler (1976) single-phase-analogy friction factor, common
    // form re-used (with different Reynolds numbers) for liquid-wall,
    // gas-wall and gas-liquid interfacial friction in Table 1.
    double taitelDuklerFactor(double Re) {
        Re = std::max(Re, constants::tiny);
        return (Re < 2100.0) ? 16.0 / Re : 0.046 * std::pow(Re, -0.2);
    }

    double speddingHandFactor(double Re, double e1ReSl) {
        Re = std::max(Re, constants::tiny);
        if (Re < 2100.0) return 24.0 / Re;
        return 0.0262 * std::pow(std::max(e1ReSl, constants::tiny), -0.139);
    }

    // Moody (1947) rough-wall friction factor, as tabulated (Table 2, Hall
    // 1957): f = 0.001375 * [1 + (2e4*(k/D) + 1e6/Re)^(1/3)].
    double moodyRoughWallFactor(double Re, double roughnessOverD) {
        Re = std::max(Re, constants::tiny);
        if (Re < 2100.0) return 16.0 / Re;
        const double bracket = 2.0e4 * roughnessOverD + 1.0e6 / Re;
        return 0.001375 * (1.0 + std::cbrt(bracket));
    }
}

double wallFrictionFactor(const WallFrictionInputs& in, FrictionCorrelation correlation) {
    const double Re = reynoldsNumber(in.density, in.velocity, in.hydraulicDiameter, in.viscosity);
    switch (correlation) {
        case FrictionCorrelation::SpeddingHand1997: {
            // e1*Re_sl in the turbulent branch is approximated here using
            // the supplied Reynolds number directly (the holdup-weighting
            // e1 must be folded into `in.velocity`/`in.hydraulicDiameter`
            // by the caller if the exact Spedding & Hand form is required).
            return speddingHandFactor(Re, Re);
        }
        case FrictionCorrelation::MoodyRoughWall:
            return moodyRoughWallFactor(Re, in.roughness / std::max(in.hydraulicDiameter, constants::tiny));
        case FrictionCorrelation::TaitelDukler1976:
        default:
            return taitelDuklerFactor(Re);
    }
}

double kelvinHelmholtzParameterF(const InterfacialFrictionInputs& in) {
    // F = (ug - ul) * sqrt( rho_g / ((rho_l - rho_g) * A2 * g * cos(theta)) * dA1/dh1 )
    // Eq. (19).
    //
    // This F parameter models the inviscid Kelvin-Helmholtz instability of
    // a gravity-stratified interface: cos(theta) is exactly the mechanism
    // by which gravity resists the interface being displaced across the
    // pipe cross-section. As theta -> 90 deg that stabilising mechanism
    // vanishes (there is no cross-sectional "low side" left for gravity to
    // hold the interface against), so F formally diverges -- a genuine
    // physical breakdown of this near-horizontal closure, not just a
    // numerical artefact. A floor on |cos(theta)| keeps F (and hence the
    // interfacial friction factor built from it) bounded as that limit is
    // approached, consistent with this correlation only being applied, per
    // this codebase's own documented scope, as an extrapolation at steep
    // or vertical inclination -- see PipeGeometry.hpp and the README.
    const double dRho = std::max(in.densityLiquid - in.densityGas, constants::tiny);
    const double cosThetaFloored = std::max(std::fabs(std::cos(in.inclination)), 0.05);
    const double denom = dRho * std::max(in.area2, constants::tiny) * constants::gravity * cosThetaFloored;
    const double inside = (in.densityGas / std::max(denom, constants::tiny)) * in.dA1dh1;
    const double relVel = in.velocityGas - in.velocityLiquid;
    return relVel * std::sqrt(std::max(inside, 0.0));
}

double interfacialFrictionFactor(const InterfacialFrictionInputs& in, FrictionCorrelation correlation) {
    switch (correlation) {
        case FrictionCorrelation::AndreussiPersen1987: {
            const double F0 = 0.36;
            const double F = kelvinHelmholtzParameterF(in);
            if (F <= F0) return in.gasWallFriction;
            const double hOverD1 = in.liquidHeight / std::max(in.hydraulicDiameter1, constants::tiny);
            return in.gasWallFriction * (1.0 + 29.7 * std::pow(F - F0, 0.67) * std::pow(hOverD1, 0.2));
        }
        case FrictionCorrelation::AndritsosHanratty1987: {
            // U_G,crit = 5 * sqrt(rho_g(P)/rho_g(P_atm)) [Andritsos & Hanratty, 1987];
            // for an ideal gas at fixed temperature rho_g(P)/rho_g(P_atm) = P/P_atm.
            const double UGcrit = 5.0 * std::sqrt(std::max(in.localPressure, constants::tiny) /
                                                   std::max(in.atmPressure, constants::tiny));
            const double e2u2 = in.areaFraction2 * in.velocityGas;
            if (e2u2 < UGcrit) return in.gasWallFriction;
            const double hOverD1 = in.liquidHeight / std::max(in.hydraulicDiameter1, constants::tiny);
            return in.gasWallFriction * (1.0 + 15.0 * std::sqrt(hOverD1) * (e2u2 / UGcrit - 1.0));
        }
        case FrictionCorrelation::CohenHanratty1968:
            return 0.014;
        case FrictionCorrelation::Wallis1969: {
            const double e1 = in.area1 / std::max(in.area1 + in.area2, constants::tiny);
            return 0.005 * (1.0 + 75.0 * e1);
        }
        case FrictionCorrelation::TaitelDukler1976:
        default: {
            const double Re = reynoldsNumber(in.densityGas, in.velocityGas - in.velocityLiquid,
                                              in.hydraulicDiameter2, in.viscosityGas);
            return taitelDuklerFactor(Re);
        }
    }
}

double wallShearStress(double frictionFactor, double density, double velocity) {
    return 0.5 * frictionFactor * density * std::fabs(velocity) * velocity;
}

double interfacialShearStress(double frictionFactor, double densityGas, double velocityGas, double velocityLiquid) {
    const double relVel = velocityGas - velocityLiquid;
    return 0.5 * frictionFactor * densityGas * std::fabs(relVel) * relVel;
}

// ---- Table 3 --------------------------------------------------------------

double bubbleEntrainmentRate(const BubbleEntrainmentInputs& in, double pipeDiameter) {
    const double bracket = 0.076 * (in.interfacialWidth / std::max(pipeDiameter, constants::tiny)) *
                                (in.waveVelocity - in.velocityLiquid) - 0.15;
    return in.densityGas * in.area * std::max(bracket, 0.0);
}

double bubbleDisengagementRate(const BubbleDisengagementInputs& in) {
    const double dRho = std::max(in.densityLiquid - in.densityGas, 0.0);
    const double riseVelScale = 1.18 * std::pow(constants::gravity * in.surfaceTension * dRho /
                                                 (in.densityLiquid * in.densityLiquid), 0.25);
    const double el = std::clamp(in.continuousLiquidFraction, 0.0, 1.0);
    return in.densityGas * in.K * riseVelScale * in.interfacialWidth * (1.0 - el);
}

double waveCelerity(const WaveCelerityInputs& in) {
    const double Fr = in.mixtureVelocity / std::sqrt(std::max(constants::gravity * in.pipeDiameter, constants::tiny));
    double C0, u0;
    if (Fr <= 3.5) {
        C0 = 1.05 + 0.15 * std::sin(in.inclination) * std::sin(in.inclination);
        u0 = (0.35 * std::sin(in.inclination) + 0.54 * std::cos(in.inclination)) *
             std::sqrt(constants::gravity * in.pipeDiameter);
    } else {
        C0 = 1.2;
        u0 = 0.35 * std::sqrt(constants::gravity * in.pipeDiameter) * std::sin(in.inclination);
    }
    return C0 * in.mixtureVelocity + u0;
}

double dropletEntrainmentRate(const DropletEntrainmentInputs& in) {
    const double kinetic = std::sqrt(in.densityGas * in.densityLiquid) * in.velocityGas * in.velocityGas *
                            in.pipeDiameter;
    const double excess = std::max(kinetic - 100.0 * in.viscosityLiquid, 0.0);
    return (constants::pi * in.pipeDiameter / 4.0) * in.ke * std::sqrt(in.densityGas) * excess;
}

double dropletDepositionRate(const DropletDepositionInputs& in) {
    const double eg = std::max(in.continuousGasFraction, constants::small_e);
    return (in.pipeDiameter / 4.0) * in.depositionVelocity * (in.dispersedLiquidFraction / eg) * in.densityLiquid;
}

double dropletDragCoefficient(double dropletReynolds) {
    const double Red = std::max(dropletReynolds, constants::tiny);
    if (Red > 500.0) return 0.44;
    return 18.5 * std::pow(Red, -0.6);
}

double bubbleDragCoefficient(double bubbleReynolds, double eotvos) {
    const double Reb = std::max(bubbleReynolds, constants::tiny);
    const double a = (24.0 / Reb) * (1.0 + 0.15 * std::pow(Reb, 0.687));
    const double b = (8.0 / 3.0) * eotvos / (eotvos + 4.0);
    return std::max(a, b);
}

double dropletDiameter(const DropletDiameterInputs& in) {
    const double Ustar = std::max(in.frictionVelocityGas, constants::tiny);
    const double bracket = 4.848 * in.surfaceTension / (in.densityGas * Ustar * Ustar * in.pipeDiameter) + 0.0038;
    return (1.0 / Ustar) * bracket * in.pipeDiameter;
}

double bubbleDiameter(const BubbleDiameterInputs& in, double surfaceTension) {
    const double denom = 0.5 * std::max(in.wallFrictionFactorLiquid, constants::tiny) * in.densityLiquid *
                          in.velocityLiquid1 * in.velocityLiquid1;
    const double el = std::clamp(in.continuousLiquidFraction, 0.0, 1.0);
    return (surfaceTension / std::max(denom, constants::tiny)) * in.weberCritical *
           (1.0 + 51.7 * std::pow(el, 1.5));
}

double dragForcePerVolume(double dragCoefficient, double particleDiameter,
                           double continuousDensity, double dispersedFraction,
                           double relativeVelocity) {
    const double d = std::max(particleDiameter, constants::tiny);
    return 0.75 * (dragCoefficient / d) * continuousDensity * dispersedFraction *
           std::fabs(relativeVelocity) * relativeVelocity;
}

} // namespace mfs
