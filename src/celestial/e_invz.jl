# Date created: 24/02/2019
# Julia Conversion: Gerry Gralton
# Original author: Aaron Robotham

function e_invz(z::Float64,
                omegaM::Float64,
                omegaL::Float64,
                omegaR::Float64,
                omegaK::Float64,
                w0::Float64,
                wprime::Float64)

    e_invz = 1/(sqrt(omegaR*(1+z)^4 +
            omegaM * (1 + z)^3 +
            omegaK * (1 + z)^2 +
            omegaL*cosgrowRhoDE(z, w0=w0, wprime=wprime, rhoDE=1.)) * (1 + z))
    return e_invz
end
