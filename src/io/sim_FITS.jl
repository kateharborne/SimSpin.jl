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
                    out_file::String;
                    obs_name::String="SimSpin datacube")

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

    file = FITS(out_file, "w")

    head_sim = FITSHeader(keywords, headers[:, 2], comments)

    write(file, out_data, header=head_sim, name=obs_name)

    close(file)
end

"""
    sim_FITS(out_data; folder)

Can be used to export multiple datacubes in different FITS files.

Parameters:\n
    out_data    Array of (datacube, observation) tuples. As output by [`build_datacube`](@ref) when taking multiple observations.
    folder      Optional. The file path to where the output should be. Defaults to "~".
"""
function sim_FITS(out_data::Array{Tuple, 1}; folder::String="~")

    len = length(out_data)

    for tuple in out_data
        cube = tuple[1]
        obs = tuple[2]

        filename = string("SimSpin_z", obs.z, "_INC", obs.inc_deg, "_R200", obs.r200, "_M2L", obs.mass2light, "_BLURFWHM", obs.blur.fwhm)
        file = joinpath(folder, filename)

        sim_FITS(cube, obs, file)
    end
end
