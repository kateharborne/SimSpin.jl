# Date created: 10/01/2020
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

using LinearAlgebra

"""
    circular_ap(sbin)

Returns Array of size sbin * sbin
Array is populated with ones denoting aperture pixels
in circular arrangement and zeros denoting no pixel.
"""
function circular_ap(sbin::Int64)

    max_rad = (sbin / 2) - 0.5

    y = repeat((-max_rad:max_rad).^2, outer = [1, sbin])
    x = y'

    ap_region = sqrt.(x + y)

    ap_region[ap_region .<= max_rad + 0.5] .= 1
    ap_region[ap_region .> max_rad + 0.5] .= 0

    return ap_region
end

"""
    circular_ap_cut(galaxy_data, ap_size)

Returns all particles that are within circular aperature of radius ap_size.
"""
function circular_ap_cut(galaxy_data::Array{Galaxy_particle, 1}, ap_size::Float64)

    rad = ap_size/2
    filter!(part -> part.obs.r_obs < ap_size/2, galaxy_data)
    return galaxy_data
end

"""
    square_ap(sbin)

Returns Array of size sbin * sbin
Array is populated with ones denoting aperture pixels
in square arrangement.
"""
function square_ap(sbin::Int64)

    ap_region = ones(sbin, sbin)
    return ap_region
end

"""
    square_ap_cut(galaxy_data, sbin, sbinsize)

Returns all particles that are within square aperature
"""
function square_ap_cut(galaxy_data::Array{Galaxy_particle, 1}, sbin::Int64, sbinsize::Float64)
    threshold = sbin * sbinsize / 2

    filter!(part -> abs(part.x) < threshold && abs(part.obs.z_obs) < threshold, galaxy_data)
    return galaxy_data
end


"""
    hexagonal_ap(sbin)

Returns Array of size sbin * sbin
Array is populated with ones denoting aperture pixels
in hexagonal arrangement and zeros denoting no pixel.
"""
function hexagonal_ap(sbin::Int64)

    max_rad = (sbin / 2) - 0.5

    quart = sbin / 4
    qsqrt3 = quart * sqrt(3)

    x = Vector{Float64}(abs.(-max_rad:max_rad))     #row
    y = x'                                          #column

    len = length(x)

    ap_region = zeros(len, len)
    rr = @. (2 * quart * qsqrt3) - (quart * y) - (qsqrt3 * x)

    ap_region[rr .>= 0] .= 1

    return ap_region
end

"""
    hexagonal_ap_cut(galaxy_data, sbin, sbinsize)

Returns all particles that are within hexagonal aperature
"""
function hexagonal_ap_cut(galaxy_data::Array{Galaxy_particle, 1}, sbin::Int64, sbinsize::Float64)

    quart = sbin / 4
    qsqrt3 = quart * sqrt(3)
    threshold = sbin * sbinsize
    half_threshold = threshold / 2
    vert_threshold = threshold * sqrt(3) / 4

    abs_z_obs = abs.(getfield.(getfield.(galaxy_data, :obs), :z_obs))
    abs_x = abs.(getfield.(galaxy_data, :x))

    indexs = findall(ind -> abs_x[ind] < half_threshold && abs_z_obs[ind] .< vert_threshold, 1:length(galaxy_data))

    trimmed = galaxy_data[indexs]
    trimmed_abs_x = abs_x[indexs]
    trimmed_abs_z_obs = abs_z_obs[indexs]

    dotprod = @. (2 * quart * sbinsize * qsqrt3 * sbinsize) - (quart * sbinsize * trimmed_abs_z_obs) - (qsqrt3 * sbinsize * trimmed_abs_x)

    trimmed = trimmed[dotprod .>= 0]
    return trimmed
    #TODO: probably uncessary. just multiply parts_in_cell by ap_region
end
