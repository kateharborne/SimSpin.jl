# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

using Interpolations

"""
    assign_flux(particle, filter, redshift)

    Constructs a flux profile using a Julia implementation of `ProSpect`
    for a `Galaxy_ssp` particle when observed with specified filter and redshift.

    Filter to be created using `get_filter()`.
"""
function assign_flux(particle::Galaxy_ssp,
                        filter::Interpolations.FilledExtrapolation{},
                        z::Float64,
                        redshiftCoef::Float64,
                        lum_dist::Float64,
                        mass2light::Float64)

    spectra = part_spectra(particle)
    flux = ProSpect.photom_lum(spectra, filter, z, lum_dist)

    return flux
end

function assign_flux(particle::Galaxy_ssp,
                    filter::Nothing,
                    z::Float64,
                    redshiftCoef::Float64,
                    lumDist::Float64,
                    mass2light::Float64)
    error("IFU filter must be specified if SSP particles are used")
end

"""
    assign_flux(particle, filter, redshift)

    Constructs a flux profile using for a `Galaxy_lum` particle using a mass to flux conversion.
"""
function assign_flux(particle::Galaxy_lum,
                        filter::Union{Nothing, Interpolations.FilledExtrapolation{}},
                        z::Float64,
                        redshiftCoef::Float64,
                        lumDist::Float64,
                        mass2light::Float64)

    flux = mass_to_flux(particle, redshiftCoef, mass2light)
    return flux
end

"""
    assign_flux(particle, filter, redshift)

    Assigns a flux profile of zero for a `Galaxy_dark` particle (ie dark matter or gas).
"""
function assign_flux(particle::Galaxy_dark,
                    filter::Union{Nothing, Interpolations.FilledExtrapolation{}},
                    z::Float64,
                    redshiftCoef::Float64,
                    lumDist::Float64,
                    mass2light::Float64)

    flux = 0
    return flux
end

"""
    part_spectra(particle)

Returns spectra for input particle of type Galaxy_ssp.
Julia implementation of R package, ProSpect.
"""
function part_spectra(particle::Galaxy_ssp)

    metallicity = particle.metallicity
    age = particle.age * 1e9
    mass = particle.mass * 1e10

    Z = ProSpect.interp_param([metallicity], ProSpect.BC03lr["Z"], log = true)
    A = ProSpect.interp_param([age], ProSpect.BC03lr["Age"], log = true)

    weights = Dict(     "hihi" => Z["weight_hi"][1] * A["weight_hi"][1],
                        "hilo" => Z["weight_hi"][1] * A["weight_lo"][1],
                        "lohi" => Z["weight_lo"][1] * A["weight_hi"][1],
                        "lolo" => Z["weight_lo"][1] * A["weight_lo"][1])

    part_spec = zeros(length(ProSpect.BC03lr["Wave"]))
    part_spec = (   (ProSpect.BC03lr["ZSpec"][Z["ID_hi"][1]][A["ID_hi"][1],:] * weights["hihi"]) +
                    (ProSpect.BC03lr["ZSpec"][Z["ID_hi"][1]][A["ID_lo"][1],:] * weights["hilo"]) +
                    (ProSpect.BC03lr["ZSpec"][Z["ID_lo"][1]][A["ID_hi"][1],:] * weights["lohi"]) +
                    (ProSpect.BC03lr["ZSpec"][Z["ID_lo"][1]][A["ID_lo"][1],:] * weights["lolo"])) * mass

    return part_spec
end

"""
    mass_to_flux(particle, redshiftCoef, m2l)

Converted particle mass to a flux value in jansky using the provided mass to light ratio.
"""
function mass_to_flux(particle::Galaxy_lum, redshiftCoef::Float64, m2l::Float64)

    mass = particle.mass
    lum = mass * 1e10 / m2l      #Convert particle masses to cell luminosity

    flux = lum * redshiftCoef * ProSpect.cgs_to_jansky

    return flux
end
