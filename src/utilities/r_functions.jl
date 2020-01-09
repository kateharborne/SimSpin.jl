# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

using RCall
using StaticArrays
using CodecXz
using RData

"""
function sftToAge(sft)

    Returns particle age from stellar formation time input.
    Uses R package, celestial.
"""
function sftToAge(sft::Array{Float64,1})
    R"library(celestial)"
    r_ages = R"cosdistTravelTime($sft)"
    ages = rcopy(r_ages)

    return ages
end

"""
function angleSize(z)

    Returns angular galaxy size in kpc from redshift, z, input.
    Uses R package, celestial.
"""
function angleSize(z::Float64)

    ang_size = cosdistAngScale(z, ref="Planck")

    return ang_size
end

"""
function cosdistLumDist(z)

    Returns luminosity distance in Mpc from redshift, z, input.
    Uses julia implementation of R package, celestial.
"""
function lumDist(z::Float64=0.1,
                        H0::Float64=67.8,
                        omegaM::Float64=0.308,
                        omegaL::Float64=1-omegaM,
                        ref::String="Planck")

    lum_dist = cosdistLumDist(z=z,
                                H0=H0,
                                omegaM=omegaM,
                                omegaL=omegaL,
                                ref=ref)

    return lum_dist
end

"""
function photom_lum(spectra, filter, z)


"""
function photom_lum(spectra, filter, z::Float64)

    particle_flux = photom_lum(wave=BC03lr["Wave"],
                            lum=spectra,
                            filter=filter,
                            outtype="Janksy",
                            z=redshift,
                            ref="Planck")

    return particle_flux
end

"""
function part_spectra(particle)

    Returns spectra for input particle of type Galaxy_ssp.
    Julia implementation of R package, ProSpect.
"""
function part_spectra(particle::Galaxy_ssp)

    metallicity = particle.metallicity
    age = particle.age
    mass = particle.mass

    Z = interp_param(metallicity, BC03lr["Z"], log = true)
    A = interp_param(age, BC03lr["Age"], log = true)

    weights = Dict(     "hihi" => Z["weight_hi"] * A["weight_hi"],
                        "hilo" => Z["weight_hi"] * A["weight_lo"],
                        "lohi" => Z["weight_lo"] * A["weight_hi"],
                        "lolo" => Z["weight_lo"] * A["weight_lo"])

    part_spec = zeros(MVector{length(BC03lr["Wave"])})
    part_spec = (   (BC03lr["ZSpec"][Z["ID_hi"]][A["ID_hi"],:] * weights["hihi"]) +
                    (BC03lr["ZSpec"][Z["ID_hi"]][A["ID_lo"],:] * weights["hilo"]) +
                    (BC03lr["ZSpec"][Z["ID_lo"]][A["ID_hi"],:] * weights["lohi"]) +
                    (BC03lr["ZSpec"][Z["ID_lo"]][A["ID_lo"],:] * weights["lolo"])) * mass

    return part_spec
end

#TODO: m2l ratio needs user input
function mass_to_flux(particle::Galaxy_lum, redshiftCoef::Float64; m2l::Float64=1.)

    mass = particle.mass
    lum = mass * 1e10 / m2l      #Convert particle masses to cell luminosity

    flux = lum * redshiftCoef * cgs_to_jansky

    return flux
end
