#pragma once

namespace mfs {

// Fluid property set. Liquid is treated as incompressible (constant density);
// gas obeys the isothermal ideal-gas law, consistent with the pressure
// equation (Eq. 22 in the paper) which requires d(rho_g)/dP.
struct FluidProperties {
    double rhoLiquid   = 998.0;   // [kg/m^3]
    double muLiquid    = 1.0e-3;  // [Pa.s]
    double sigma       = 0.072;   // gas-liquid surface tension [N/m]

    double gasConstant = 287.0;   // specific gas constant [J/(kg.K)]
    double temperature = 293.0;   // [K]
    double muGas       = 1.8e-5;  // [Pa.s]

    double rhoGas(double pressure) const {
        return pressure / (gasConstant * temperature);
    }

    // d(rho_g)/dP for the isothermal ideal-gas EOS, used to linearize the
    // pressure equation.
    double drhoGasdP() const {
        return 1.0 / (gasConstant * temperature);
    }
};

} // namespace mfs
