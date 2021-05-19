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
