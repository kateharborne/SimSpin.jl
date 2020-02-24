# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    sim_to_galaxy(sim_data)

Returns particle data wrt centre of galaxy in both spatial and velocity space instead of wrt arbitrary point in simulation.
"""
function sim_to_galaxy(sim_data::Array{Sim_particle, 1})

    galaxy_data = Array{Galaxy_particle, 1}(undef, length(sim_data))

    lum_data = sim_data[findall(x -> typeof(x) != Sim_dark, sim_data)]   #get all luminous particles
    centre = galaxy_centre(lum_data)                    #use luminous particles only to find centre
    rot_mat = galaxy_orient(lum_data, centre)           #use luminous particles only to find rotation matrix

    Threads.@threads for index in eachindex(sim_data)
        galaxy_data[index] = galaxy_particle(sim_data[index], centre, rot_mat)  #convert each particle from Sim_particle to Galaxy_particle
    end

    return galaxy_data
end

function sim_to_galaxy(galaxy_data::Array{Galaxy_particle, 1})

    @warn("This galaxy was already in the galaxy reference frame.")
    return galaxy_data
end
