# Date created: 16/06/2020
# Julia Conversion: Gerry Gralton
# Original author: Aaron Robotham

using OffsetArrays
using HCubature
using Distributions
using SpecialFunctions

"""
    sersic_cuba(nser, re, mag, psf_dim)

Returns a sersic kernel of size psf_dim * psf_dim.
Integrates a sersic function over each pixel.
"""
function sersic_cuba(nser::Real, re::Real, mag::Real, psf_dim::Int64)

    gamma_dist = Distributions.Gamma(nser * 2)
    bn  = quantile(gamma_dist, 0.5)

    ind::Int64 = floor(psf_dim/2)
    kernel = zeros(Float64, psf_dim, psf_dim)
    kernel = OffsetArray(kernel, -ind:ind, -ind:ind)

    func = sersic_construct(nser, re, bn)
    for coords in CartesianIndices(kernel)
        xmin = coords[1] - 0.5
        xmax = coords[1] + 0.5

        ymin = coords[2] - 0.5
        ymax = coords[2] + 0.5

        kernel[coords] = hcubature(func, [xmin, ymin], [xmax, ymax]; rtol = 1e-3, atol=1e-10)[1]
    end

    return kernel .* sersic_scale(nser, re, mag, bn)
end

function sersic_construct(nser::Real, re::Real, bn::Real)
    inv_nser = 1/nser

    func =  function(coords)
                rad = sqrt(sum(coords .^ 2))
                return exp(-bn * ( (rad/re)^inv_nser - 1))
            end
    return func
end

function sersic_scale(nser::Real, re::Real, mag::Real, bn::Real)
    #TODO: remove SpecialFunctions dependency
    
    lumtot = re^2 * 2 * π * nser *(exp(bn)/(bn^(2*nser)))*SpecialFunctions.gamma(2*nser)
    magtot = -2.5 * log10(lumtot)
    factor = 1/(10^(0.4 * (mag - magtot)))
    return factor
end
