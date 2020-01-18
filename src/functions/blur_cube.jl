# Date created: 13/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

using ImageFiltering
using StaticArrays

function blur_cube(ifu_cube::Array{Float64, 3},
                    blur::Blur,
                    ap_region::MArray,
                    sbin::Int64,
                    vbin::Int64,
                    ang_size::Float64,
                    sbinsize::Float64)

    sd_scaled = blur.sigma * ang_size / sbinsize # sd scaled to image pixel dimensions

    if (sbin < 25 && (sbin % 2) != 0)
        psf_dim = sbin
    elseif (sbin < 25 && (sbin % 2) == 0)
        psf_dim = sbin-1
    else psf_dim = 25
    end

    if lowercase(blur.psf) == "gaussian"
        psf_k = Kernel.gaussian([sd_scaled], [psf_dim])
    #elseif lowercase(blur.psf) == "moffat"
        #psf_k = profitCubaMoffat(fwhm = fwhm_scaled, mag = 1, con = 5, dim = [psf_dim,psf_dim])
    else
        error("Blur PSF type:", blur.psf, "is not supported.")
    end

    kernel = kernelfactors((psf_k, psf_k))
    blur_cube = zeros(Float64, (sbin, sbin, vbin))

    for bin in 1:vbin
        blur_cube[:,:,bin] = imfilter(ifu_cube[:,:,bin], kernel)
        blur_cube[:,:,bin] = blur_cube[:,:,bin] .* ap_region #Set all cells outside aperture to zero
    end

    return blur_cube
end
