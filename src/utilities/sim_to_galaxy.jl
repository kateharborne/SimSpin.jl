# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
function sim_to_galaxy(sim_data)

    Returns particle data wrt to galaxy reference frame rather
    than simulation reference frame
"""
function sim_to_galaxy(sim_data::Array{Sim_particle, 1})

    galaxy_data = Array{Galaxy_particle, 1}(undef, length(sim_data))
    centre = galaxy_centre(sim_data)
    rot_mat = galaxy_orient(sim_data, centre)

    Threads.@threads for index in eachindex(sim_data)
        new = galaxy_particle(sim_data[index], centre, rot_mat)
        galaxy_data[index] = new
    end
    sim_data = nothing

    return galaxy_data
end
