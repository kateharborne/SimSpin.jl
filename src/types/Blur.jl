# Date created: 13/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    Blur(psf; sigma, fwhm)

    Create a struct containing the blurring to be used.
    Must specify blurring point spread function, `psf`,
    currently only "Gaussian" is supported.

    Keyword arguments (one must be specified, sigma is prioritised):\n
        sigma       The standard deviation of the point spread function
        fwhm        The full width half max of the point spread function
"""
struct Blur
    psf::String
    sigma::Float64

    function Blur(psf::String;
                    sigma::Float64=-1.,
                    fwhm::Float64=-1.)

        if (lowercase(psf) != "moffat" && lowercase(psf) != "gaussian")
            error("Blur PSF type:", blur.psf, "is not supported.")
        elseif sigma < 0 && fwhm < 0
            error("Please specify either standard deviation, sigma, or full width half max, fwhm, as a positive float.")
        elseif !isnothing(sigma)
            new(psf, sigma)
        else
            sigma = fwhm/2.355
            new(psf, sigma)
        end
    end
end
