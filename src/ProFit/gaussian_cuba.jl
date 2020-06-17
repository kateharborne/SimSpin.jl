# Date created: 16/06/2020
# Julia Conversion: Gerry Gralton
# Original author: Aaron Robotham

using OffsetArrays

"""
    gaussian_cuba(fwhm, psf_dim)

Returns a gaussian kernel of size psf_dim * psf_dim.
Integrates a sersic function with nser=0.5 over each pixel.
"""
function gaussian_cuba(fwhm::Real, psf_dim::Int64)
    nser = 0.5
    re = fwhm/2
    mag = 0.
    kernel = sersic_cuba(nser, re, mag, psf_dim)

    return gaussian_scale(kernel, psf_dim)
end

function gaussian_scale(kernel::OffsetArrays.OffsetArray{Float64,2,Array{Float64,2}}, psf_dim::Int64)
    size = psf_dim^2

    trim = 1 - π / 4
    cut = floor(trim * size / 8) * 8 / size

    kernel_vec = reshape(kernel, size)
    kernel[kernel .< quantile(kernel_vec, cut)] .= 0.

    kernel = kernel / sum(kernel)
    return kernel
end
