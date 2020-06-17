# Date created: 16/06/2020
# Julia Conversion: Gerry Gralton
# Original author: Aaron Robotham

using OffsetArrays
using HCubature

"""
    cuba_moffat(β, α, mag, psf_dim)

Returns a moffat kernel of size psf_dim * psf_dim.
Integrates a moffat function over each pixel.
"""
function moffat_cuba(β::Real, α::Real, mag::Real, psf_dim::Int64)

    ind::Int64 = floor(psf_dim/2)
    kernel = zeros(Float64, psf_dim, psf_dim)
    kernel = OffsetArray(kernel, -ind:ind, -ind:ind)

    func = moffat_construct(β, α)
    for coords in CartesianIndices(kernel)
        xmin = coords[1] - 0.5
        xmax = coords[1] + 0.5

        ymin = coords[2] - 0.5
        ymax = coords[2] + 0.5

        kernel[coords] = hcubature(func, [xmin, ymin], [xmax, ymax]; rtol = 1e-3, atol=1e-10)[1]
    end

    return kernel .* moffat_scale(β, α, mag)
end

function moffat_construct(β::Real, α::Real)
    α2 = α^2
    func =  function(coords)
                rad2 = sum(coords .^ 2)
                return ((1 + rad2/α2)^-β)
            end
    return func
end

function moffat_scale(β::Real, α::Real, mag::Real)
    lumtot = pi * α^2/(β - 1)
    magtot = -2.5 * log10(lumtot)
    factor = 1/(10^(0.4 * (mag - magtot)))
    return factor
end
