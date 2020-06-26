# Date created: 23/01/2020
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
    values::Array{Any,1} = headers[:, 2]
    comments::Array{String,1} = headers[:, 3]

    file = FITS(out_file, "w")

    head_sim = FITSHeader(keywords, values, comments)

    write(file, out_data, header=head_sim, name=obs_name)

    close(file)

    println("Completed output of ", out_file, ".")
end

"""
    sim_FITS(data_cube,
                out_file)

Outputs FITS file containing 3D, velocity binned datacube.

Parameters:\n
    data_cube   Tuple with simulated IFU datacube and Observation as output from `build_datacube()`.
    out_file    String denoting path and name of FITS file to be output.
"""
function sim_FITS(out_data::Tuple{Array{Float64,3},Observation},
                    out_file::String;
                    obs_name::String="SimSpin datacube")
    return sim_FITS(out_data[1], out_data[2], out_file, obs_name=obs_name)
end

"""
    sim_FITS(out_data;
                folder = "~",
                file_prefix = "SimSpin")

Can be used to export multiple datacubes in different FITS files.

Parameters:\n
    out_data    Array of (datacube, observation) tuples. As output by `build_datacube()` when taking multiple observations.
    folder      Optional. The file path to where the output should be. Defaults to "~".
    file_prefix Optional. Prefix for the start of the .fits file's name.
"""
function sim_FITS(out_data::Array{Tuple{Array{Float64,3},Observation}, 1}; folder::String="~", file_prefix::String="SimSpin")

    for tuple in out_data
        cube = tuple[1]
        obs = tuple[2]

        if isnothing(obs.mass2light) && isnothing(obs.blur)
            filename = string(file_prefix, "_z", obs.z, "_INC", obs.inc_deg, "_R200", obs.r200)
        elseif !isnothing(obs.mass2light) && isnothing(obs.blur)
            filename = string(file_prefix, "_z", obs.z, "_INC", obs.inc_deg, "_R200", obs.r200, "_M2L", obs.mass2light)
        elseif isnothing(obs.mass2light) && !isnothing(obs.blur)
            filename = string(file_prefix, "_z", obs.z, "_INC", obs.inc_deg, "_R200", obs.r200, "_BLURFWHM", obs.blur.fwhm)
        else
            filename = string(file_prefix, "_z", obs.z, "_INC", obs.inc_deg, "_R200", obs.r200, "_M2L", obs.mass2light, "_BLURFWHM", obs.blur.fwhm)
        end

        file = joinpath(folder, filename)

        sim_FITS(cube, obs, file)
    end
end
