# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

using DelimitedFiles

function assign_flux(particle::Galaxy_ssp, filter::String, redshift::Float64)

    if(filter == "g")   tempfilt = readdlm(joinpath(@__DIR__, "../ProSpect/data/filt_g_SDSS.tab"), header=true)
    elseif(filter == "r")   tempfilt = readdlm(joinpath(@__DIR__, "../ProSpect/data/filt_r_SDSS.tab"), header=true)
    else error("An unsupported filter type was specified.")
    end

    flux = 0

    spectra = part_spectra(particle)
    flux += photom_lum(spectra, tempfilt, redshift)

    return flux
end

function assign_flux(particle::Galaxy_lum, filter, redshiftCoef::Float64)

    flux = mass_to_flux(particle, redshiftCoef)
    return flux
end

function assign_flux(particle::Galaxy_dark, filter, redshift::Float64)

    flux = 0
    return flux
end
