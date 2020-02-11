# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    build_datacube(galaxy_data, observation, ifu, blur)

Returns a simulated ifu datacube for input, `galaxy_data`, an array of `Particle`.
Observational variables used are as specified in `Observation` and `IFU` types. Blurring can be added if desired using `Blur` type.

Parameters:\n
    galaxy_data         Array of `Particle` describing galaxy
    ifu                 Struct of type `Telescope`
    z                   The projected redshift at which the observation is made.
    inc_deg             The inclination at which to observe the galaxy in degrees.
    r200                The virial radius specified in the simulation, kpc.
    blur                Optional. Struct of type `Blur`. If omitted no blurring occurs.
"""
function build_datacube(galaxy_data::Array{Galaxy_particle, 1},
                        ifu::Telescope,
                        z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        blur::Union{Blur, Nothing})

    galaxy_data, parts_in_cell, observe = obs_data_prep(galaxy_data, ifu, z, inc_deg, r200)

    fluxes = flux_grid(parts_in_cell, observe, ifu.filter)
    data_cube = ifu_cube(fluxes, parts_in_cell, observe)

    if isnothing(blur)  #No spatial blurring
        return data_cube
    else                #Blur image
        blur_imgs = blur_cube(data_cube, blur, observe)
        return blur_imgs
    end
end

function build_datacube(sim_data::Array{Sim_particle, 1},
                        ifu::Telescope,
                        z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        blur::Union{Blur, Nothing})

    galaxy_data = sim_to_galaxy(sim_data)
    return build_datacube(galaxy_data, ifu, z, inc_deg, r200, blur)
end

function build_datacube(galaxy_data::Array{<:Particle, 1},
                        ifu::Telescope,
                        z::Float64,
                        inc_deg::Real,
                        r200::Real)
    blur = nothing
    return build_datacube(galaxy_data, ifu, z, inc_deg, r200, blur)
end
