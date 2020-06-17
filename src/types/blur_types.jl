# Date created: 13/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

abstract type Blur end

mutable struct Scaled_lsf
    scaled_lsf::Float64
end

"""
    Gaussian_blur(;sigma, fwhm)

Create a struct containing seeing information.

Keyword arguments (at least one must be specified, sigma is prioritised):\n
    sigma       The standard deviation of the point spread function
    fwhm        The full width half max of the point spread function
"""
struct Gaussian_blur <: Blur
    σ::Float64
    fwhm::Float64
    scaled_fwhm::Scaled_lsf

    function Gaussian_blur(;σ::Float64=-1.,
                            fwhm::Float64=-1.)

        if σ < 0 && fwhm < 0
            error("Please specify standard deviation, σ, or full width half max, fwhm, as a positive float.")
        elseif σ > 0
            fwhm = σ * 2 * sqrt(2 * log(2))
            new(σ, fwhm, Scaled_lsf(0.))
        else
            σ = fwhm / 2 * sqrt(2 * log(2))
            new(σ, fwhm, Scaled_lsf(0.))
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
    fwhm::Float64
    scaled_α::Scaled_lsf

    function Moffat_blur(β::Float64;
                            α::Union{Float64, Nothing} = nothing,
                            fwhm::Union{Float64, Nothing} = nothing)

        if isnothing(α) && !isnothing(fwhm)
            α = fwhm / (2 * sqrt(2^(1/β) - 1))
            new(β, α, fwhm, Scaled_lsf(0.))
        elseif !isnothing(α)
            fwhm = α * 2 * sqrt(2^(1/β) - 1)
            new(β, α, fwhm, Scaled_lsf(0.))
        else
            error("Please specify either core width, α, or full width half max, fwhm.")
        end
    end
end

function scale_lsf(blur::Gaussian_blur, ang_size::Float64, sbinsize::Float64)
    blur.scaled_fwhm.scaled_lsf = blur.fwhm * ang_size / sbinsize
end

function scale_lsf(blur::Moffat_blur, ang_size::Float64, sbinsize::Float64)
    blur.scaled_α.scaled_lsf = blur.α * ang_size / sbinsize
end
