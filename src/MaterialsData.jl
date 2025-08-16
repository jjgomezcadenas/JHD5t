using Unitful
using UnitfulEquivalences
using PhysicalConstants.CODATA2018
using Interpolations

import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
    μW, mW, W,
    A, N, mol, mmol, V, L, M,
    Bq, mBq, μBq 


MU = Dict("H2O" => 0.045 * cm^2/g,
          "Fe316Ti" => 0.039 * cm^2/g,
          "Ti" => 0.038 * cm^2/g,
          "Cu" => 0.039 * cm^2/g,
          "Pb" => 0.044 * cm^2/g,
          "LXe" => 3.8e-2 * cm^2/g,
          "Kevlar" => 0.045 * cm^2/g,
          "Poly" => 1e-6 * cm^2/g,
          "PTFE" => 1e-6 * cm^2/g,
)

RHO = Dict("H2O" => 1.00 * g/cm^3,
           "Fe316Ti" => 7.87 * g/cm^3,
           "Ti" => 4.54 * g/cm^3,
           "Cu" => 8.96 * g/cm^3,
           "Pb" => 11.33 * g/cm^3,
           "LXe" => 2.98 * g/cm^3,
           "Kevlar" => 1440 * kg/m^3,
           "Poly" => 1.0 * g/cm^3,
           "PTFE" => 2.0 * g/cm^3,
)


BI214 = Dict("Fe316Ti" => 1.9e-3 * Bq/kg,
             "Ti" => 0.93e-3 * Bq/kg,
             "Cu" => 12.0e-6 * Bq/kg,
             "Cu-X" => 3.0e-6 * Bq/kg,
             "H2O" => 1.e-5 * Bq/kg,
             "LXe" => 1.e-12 * Bq/kg,
             "Pb" => 370.0e-6 * Bq/kg,
             "Kevlar" => 0.17 * Bq/kg,
             "Poly" => 62.0e-6 * Bq/kg,
             "PTFE" => 25e-6 * Bq/kg
)

TL208 = Dict("Fe316Ti" => 0.4e-3 * Bq/kg,
             "Ti" => 0.22e-3 * Bq/kg,
             "Cu" => 1.4e-6 * Bq/kg,
             "Cu-X" => 1.4e-6 * Bq/kg,
             "H2O" => 1.e-5/3 * Bq/kg,
             "LXe" => 1.e-12/3 * Bq/kg,
             "Pb" => 73e-6 * Bq/kg,
             "Kevlar" => 0.17 * Bq/kg,
             "Poly" => 8.e-6 * Bq/kg,
             "PTFE" => 10.0e-6 * Bq/kg
)




rhoxe = Dict("rho_2020" => 124.3 * kg/m^3, 
               "rho_3020" => 203.35 * kg/m^3,
               "rho_1520" =>  89.9 * kg/m^3,
               "rho_1020" =>   58.0 * kg/m^3,
               "rho_0520" =>   30.0 * kg/m^3,
               "rho_0720" =>   40.0 * kg/m^3,
               "rho_0920" =>   50.0 * kg/m^3)

function rhoxegr(rhoxe_dict)
    return Dict(k => uconvert(g/cm^3, v) for (k, v) in rhoxe_dict)
end

rhoxe_gcm3 = rhoxegr(rhoxe)