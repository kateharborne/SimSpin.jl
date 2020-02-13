# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
function sftToAge(sft)

    Returns particle age from stellar formation time input.
    Uses julia implementation of R package, celestial.
"""
function sftToAge(sft::Float64)
    ages = cosdistTravelTime((1 /sft) - 1)

    return ages
end

"""
function angleSize(z)

    Returns angular galaxy size in kpc from redshift, z, input.
    Uses julia implementation of R package, celestial.
"""
function angleSize(z::Float64)

    ang_size = cosdistAngScale(z, ref="Planck")

    return ang_size
end

"""
function lumDist(z)

    Returns luminosity distance in Mpc from redshift, z, input.
    Uses julia implementation of R package, celestial.
"""
function lumDist(;z::Float64=0.1,
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

    Calls julia implementation of ProSpect's photom_lum function
"""
function photom_lum(spectra::Array{Float64,1},
                    filter::Interpolations.FilledExtrapolation{},
                    z::Float64,
                    lum_dist::Float64)

    particle_flux = photom_lum(BC03lr["Wave"],
                                spectra,
                                filter,
                                outtype="Janksy",
                                z=z,
                                ref="Planck",
                                lum_dist=lum_dist)

    return particle_flux
end

"""
function part_spectra(particle)

    Returns spectra for input particle of type Galaxy_ssp.
    Julia implementation of R package, ProSpect.
"""
function part_spectra(particle::Galaxy_ssp)

    metallicity = particle.metallicity
    age = particle.age * 1e9
    mass = particle.mass * 1e10

    Z = interp_param([metallicity], BC03lr["Z"], log = true)
    A = interp_param([age], BC03lr["Age"], log = true)

    weights = Dict(     "hihi" => Z["weight_hi"][1] * A["weight_hi"][1],
                        "hilo" => Z["weight_hi"][1] * A["weight_lo"][1],
                        "lohi" => Z["weight_lo"][1] * A["weight_hi"][1],
                        "lolo" => Z["weight_lo"][1] * A["weight_lo"][1])

    part_spec = zeros(length(BC03lr["Wave"]))
    part_spec = (   (BC03lr["ZSpec"][Z["ID_hi"][1]][A["ID_hi"][1],:] * weights["hihi"]) +
                    (BC03lr["ZSpec"][Z["ID_hi"][1]][A["ID_lo"][1],:] * weights["hilo"]) +
                    (BC03lr["ZSpec"][Z["ID_lo"][1]][A["ID_hi"][1],:] * weights["lohi"]) +
                    (BC03lr["ZSpec"][Z["ID_lo"][1]][A["ID_lo"][1],:] * weights["lolo"])) * mass

    return part_spec
end

#TODO: m2l ratio needs user input
function mass_to_flux(particle::Galaxy_lum, redshiftCoef::Float64; m2l::Float64=1.)

    mass = particle.mass
    lum = mass * 1e10 / m2l      #Convert particle masses to cell luminosity

    flux = lum * redshiftCoef * cgs_to_jansky

    return flux
end
