# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    build_datacube(galaxy_data, ifu, envir)

Returns a simulated ifu datacube for input, `galaxy_data`, an array of `Particle` and a summary of observation properties used.
The instrument to be used for observation is given in the `Telescope` class parameter
and the environmental variables (redshift, inclination, seeing, etc) are given by the `Environment` class parameter.

Parameters:\n
    galaxy_data         Array of `Particle` describing galaxy
    ifu                 Struct of type `Telescope`
    envir               Struct of type `Environment`

Returns:\n
    data_cube           3D array of a simulated IFU datacube with spatial and velocity binned fluxes
    observe             A struct of type `Observation` summarising the observational properties used.
"""
function build_datacube(galaxy_data::Array{Galaxy_particle, 1},
                        ifu::Telescope,
                        envir::Environment)

    galaxy_data, parts_in_cell, observe = obs_data_prep(galaxy_data, ifu, envir)

    fluxes = flux_grid(parts_in_cell, observe, ifu.filter)
    data_cube = ifu_cube(fluxes, parts_in_cell, observe)

    if isnothing(observe.blur)    #No spatial blurring
        return data_cube, observe
    else                        #Blur image
        blur_imgs = blur_cube(data_cube, observe)
        return blur_imgs, observe
    end
end

function build_datacube(sim_data::Array{Sim_particle, 1},
                        ifu::Telescope,
                        envir::Environment)

    galaxy_data = sim_to_galaxy(sim_data)
    return build_datacube(galaxy_data, ifu, envir)
end
