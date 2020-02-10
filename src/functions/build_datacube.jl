# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    build_datacube(galaxy_data, observation, ifu, blur)

Returns a simulated ifu datacube for input, `galaxy_data`, an array of `Particle`.
Observational variables used are as specified in `Observation` and `IFU` types. Blurring can be added if desired using `Blur` type.

Parameters:\n
    galaxy_data         Array of `Particle` describing galaxy
    observation         Struct of type `Observation`
    ifu                 Struct of type `Telescope`
    blur                Optional. Struct of type `Blur`. If omitted no blurring occurs.
"""
function build_datacube(galaxy_data::Array{Galaxy_particle, 1},
                        obs::Observation,
                        ifu::Telescope,
                        blur::Union{Blur, Nothing})

    galaxy_data,
    parts_in_cell,
    ap_region, sbin,
    vbin, vseq, lsf_size,
    ang_size, sbin_size = obs_data_prep(galaxy_data, obs, ifu)

    fluxes = flux_grid(parts_in_cell, ap_region, sbin, vbin, obs.z, ifu.filter)
    data_cube = ifu_cube(fluxes, parts_in_cell, sbin, vbin, vseq, lsf_size)

    if isnothing(blur)  #No spatial blurring
        return data_cube
    else                #Blur image
        blur_imgs = blur_cube(data_cube, blur, ap_region, sbin, vbin, ang_size, sbin_size)
        return blur_imgs
    end
end


function build_datacube(sim_data::Array{Sim_particle, 1},
                        obs::Observation,
                        ifu::Telescope,
                        blur::Union{Blur, Nothing})

    galaxy_data = sim_to_galaxy(sim_data)
    return build_datacube(galaxy_data, obs, ifu, blur)
end

function build_datacube(galaxy_data::Array{<:Particle, 1},
                        obs::Observation,
                        ifu::Telescope)
    blur = nothing
    return build_datacube(galaxy_data, obs, ifu, blur)
end
