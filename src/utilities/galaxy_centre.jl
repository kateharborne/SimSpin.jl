# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

using Statistics

"""
function galaxy_centre(sim_data)

    Returns median of galaxy in velocity and spatial domains.
"""
function galaxy_centre(sim_data::Array{Sim_particle, 1})

    x = getfield.(sim_data, :x)
    y = getfield.(sim_data, :y)
    z = getfield.(sim_data, :z)

    vx = getfield.(sim_data, :vx)
    vy = getfield.(sim_data, :vy)
    vz = getfield.(sim_data, :vz)

    centre = [median(x), median(y), median(z), median(vx), median(vy), median(vz)]

    return centre
end
