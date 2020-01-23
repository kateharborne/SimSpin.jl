# Date created: 13/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

abstract type Blur end

"""
    Gaussian_blur(;sigma, fwhm)

Create a struct containing seeing information.

Keyword arguments (at least one must be specified, sigma is prioritised):\n
    sigma       The standard deviation of the point spread function
    fwhm        The full width half max of the point spread function
"""
struct Gaussian_blur <: Blur
    sigma::Float64

    function Gaussian_blur(;sigma::Float64=-1.,
                            fwhm::Float64=-1.)

        if sigma < 0 && fwhm < 0
            error("Please specify standard deviation, sigma, or full width half max, fwhm, as a positive float.")
        elseif sigma > 0
            new(sigma)
        else
            sigma = fwhm/2.355
            new(sigma)
        end
    end
end

"""
    Moffat_blur(β;
                α,
                fwhm)

Create a struct containing seeing information. β and either α or fwhm must be specified.
If both α and fwhm are specified, α is prioritised.

Arguments:\n
    β           The power component in the Moffat distribution
    α           The core width of the Moffat distribution (optional)
    fwhm        The full width half max of the Moffat distribution (optional)
"""
struct Moffat_blur <: Blur
    β::Float64
    α::Float64
    function Moffat_blur(β::Float64;
                            α::Union{Float64, Nothing} = nothing,
                            fwhm::Union{Float64, Nothing} = nothing)

        if isnothing(α) && !isnothing(fwhm)
            α = fwhm / (2 * sqrt(2^(1/β) - 1))
            new(β, α)
        elseif !isnothing(α)
            new(β, α)
        else
            error("Please specify either core width, α, or full width half max, fwhm.")
        end
    end
end
