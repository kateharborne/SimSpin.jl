# Date created: 23/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

using FITSIO

"""
    sim_FITS(data_cube,
                observe,
                out_file)

Outputs FITS file containing 3D, velocity binned datacube.

Parameters:\n
    data_cube   3D array of a simulated IFU datacube with spatial and velocity binned fluxes
    observe     A struct of type `Observation` summarising the observational properties used.
    out_file    String denoting path and name of FITS file to be output.
"""
function sim_FITS(out_data::Array{Float64,3},
                    observe::Observation,
                    out_file::String)

    if out_file[end-4:end] != ".fits"       #Append .fits to file name
        out_file = string(out_file, ".fits")
    end

    if isnothing(observe.blur); blurfwhm = 0
    else blurfwhm = observe.blur.fwhm
    end

    headers = ["REDSHIFT" observe.z "redshift, z";
               "SPIXSIZE" observe.sbinsize "spatial size, kpc/pixel";
               "PSFFWHM" blurfwhm "FWHM, arcsec";
               "SIMR200" observe.r200 "virial radius from sim, kpc"]
    keywords::Array{String,1} = headers[:, 1]
    comments::Array{String,1} = headers[:, 3]

    file = FITS(out_file, "r+")

    head_sim = FITSHeader(keywords, headers[:, 2], comments)

    write(file, out_data, header=head_sim, name=obs_name)

    close(file)
end
