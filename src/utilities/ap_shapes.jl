# Date created: 10/01/2019
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

    trimmed = galaxy_data[findall(part -> part.obs.r_obs < ap_size/2, galaxy_data)]
    return trimmed
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

    trimmed = galaxy_data[findall(part -> abs(part.x) < threshold)]
    trimmed = trimmed[findall(part -> (abs(part.obs.z_obs)) < threshold)]
    return trimmed
    #TODO: check if can be one line
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
    threshold = sbin * sbinsize

    trimmed = galaxy_data[findall(part -> abs(part.x) < threshold / 2, galaxy_data)]
    trimmed = trimmed[findall(part -> abs(part.obs.z_obs) < threshold * sqrt(3) / 4, trimmed)]

    z_obs = getfield.(getfield.(trimmed, :obs), :z_obs)
    x = getfield.(trimmed, :x)

    dotprod = @. (2 * (sbin / 4) * sbinsize * (sbin * sqrt(3) / 4) * sbinsize) - ((sbin / 4) * sbinsize) * abs(z_obs) - ((sbin * sqrt(3) / 4) * sbinsize) * abs(x)

    trimmed = trimmed[dotprod .>= 0]
end
