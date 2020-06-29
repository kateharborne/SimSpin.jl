# Date created: 10/01/2020
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    sim_to_galaxy(sim_data; centre)

Returns particle data wrt centre of galaxy in both spatial and velocity space instead of wrt arbitrary point in simulation.
Centre in position or centre in both position and velocity can be optionally specified. Else, the median of each distribution will be used as centre.
"""
function sim_to_galaxy(sim_data::Array{Sim_particle, 1}; centre::Array{Float64, 1} = Float64[])

    galaxy_data = Array{Galaxy_particle, 1}(undef, length(sim_data))

    lum_data = sim_data[findall(x -> typeof(x) != Sim_dark, sim_data)]   #get all luminous particles
    if isempty(centre)
        centre = galaxy_centre(lum_data)                    #use luminous particles only to find centre
    elseif length(centre) == 3
        vel_centre = galaxy_centre(lum_data)[4:6]
        centre = [centre; vel_centre]
    elseif length(centre) != 6
        error("Unable to process centre. Expected format [xcen, ycen, zcen] OR [xcen, ycen, zcen, vxcen, vycen, vzcen].")
    end
    rot_mat = galaxy_orient(lum_data, centre)           #use luminous particles only to find rotation matrix

    Threads.@threads for index in eachindex(sim_data)
        galaxy_data[index] = galaxy_particle(sim_data[index], centre, rot_mat)  #convert each particle from Sim_particle to Galaxy_particle
    end

    return galaxy_data
end

function sim_to_galaxy(galaxy_data::Array{Galaxy_particle, 1}; centre::Array{Float64, 1} = Float64[])

    @warn("This galaxy was already in the galaxy reference frame.")
    return galaxy_data
end
