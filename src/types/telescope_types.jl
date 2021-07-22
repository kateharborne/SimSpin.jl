# Date created: 10/02/2020
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

using Interpolations

abstract type Telescope end

"""
    IFU(fov, ap_shape, central_wvl, lsf_fwhm, pixel_sscale, pixel_vscale,  filter)

Creates a customisable, mock IFU telescope which can be used to "observe" simulations.

Parameters:\n
    fov             The field of view of the IFU, diameter in arcseconds.
    ap_shape        The shape of the field of view, with options "circular", "square" or "hexagonal".
    central_wvl     The central filter wavelength used for the observation, given in angstroms.
    lsf_fwhm        The line spread function full-width half-max, given in angstroms.
    pixel_sscale    The corresponding spatial pixel scale associated with a given telescope output in arcseconds.
    pixel_vscale    The corresponding velocity pixel scale associated with a given telescope filter output in angstroms.
    filter          Optional. If particles type is ssp, the filter within which the SED is generated. Options include "r" and "g"  for SDSS-r and SDSS-g bands respectively.
"""
struct IFU <: Telescope
    fov::Real
    ap_shape::String
    central_wvl::Real
    lsf_fwhm::Float64
    lsf_size::Float64
    pixel_sscale::Float64
    pixel_vscale::Float64
    vbinsize::Float64
    sbin::Int64
    ap_region::Array{Float64, 2}
    filter::Union{String, Nothing}
    filter_value::Union{Interpolations.FilledExtrapolation, Nothing}

    function IFU(fov::Real,
                    ap_shape::String,
                    central_wvl::Real,
                    lsf_fwhm::Float64,
                    pixel_sscale::Float64,
                    pixel_vscale::Float64,
                    filter::Union{String, Nothing})


        if isnothing(filter)
            filter_value = nothing
        elseif filter != "r" && filter != "g"
            error("Currently supported filters are only either \"r\", \"g\" or, for use without SSP, simply left blank.")
        else
            filter_value = ProSpect.get_filter(filter)
        end

        vbinsize    = (pixel_vscale / central_wvl) * (3e8 / 1e3)                        # km/s per velocity bin
        lsf_size    = ((lsf_fwhm / central_wvl) * (3e8 / 1e3)) / (2 * sqrt(2*log(2)))   # velocity uncertainty (sd)
        sbin::Int64 = floor(fov / pixel_sscale)                                         # number of spatial bins

        if (ap_shape == "circular")       # circular apperture mask
          ap_region = circular_ap(sbin)
        elseif (ap_shape == "square")     # square apperture mask
          ap_region = square_ap(sbin)
        elseif (ap_shape == "hexagonal")  # hexagonal apperture mask
          ap_region = hexagonal_ap(sbin)
        else
           error("Unsupported aperture shape specified.")
        end

        new(fov, ap_shape, central_wvl, lsf_fwhm, lsf_size, pixel_sscale, pixel_vscale, vbinsize, sbin, ap_region, filter, filter_value)
    end

    function IFU(fov::Real,
                    ap_shape::String,
                    central_wvl::Real,
                    lsf_fwhm::Float64,
                    pixel_sscale::Float64,
                    pixel_vscale::Float64)
        IFU(fov, ap_shape, central_wvl, lsf_fwhm, pixel_sscale, pixel_vscale, nothing)
    end
end


"""
    SAMI(;filter)

Creates an IFU using parameters of the SAMI survey.
Optional filters "r" or "g" can also be specified for use with SSP particles.
"""
function SAMI(;filter::Union{String, Nothing} = nothing)
    return IFU(15,
                "circular",
                4800,
                2.65,
                0.5,
                1.04,
                filter)
end

"""
    MaNGA(;filter)

Creates an IFU using parameters of the MaNGA survey.
Optional filters "r" or "g" can also be specified for use with SSP particles.
"""
function MaNGA(;filter::Union{String, Nothing} = nothing)
    return IFU(22,
                "hexagonal",
                4950,
                2.8,
                0.25,
                1.2,
                filter)
end

"""
    MUSE(;filter)

Creates an IFU using parameters of the MUSE survey.
Optional filters "r" or "g" can also be specified for use with SSP particles.
"""
function MUSE(;filter::Union{String, Nothing} = nothing)
    return IFU(15,
                "square",
                6975,
                2.63,
                0.2,
                1.25,
                filter)
end

"""
    CALIFA(;filter)

Creates an IFU using parameters of the CALIFA survey.
Optional filters "r" or "g" can also be specified for use with SSP particles.
"""
function CALIFA(;filter::Union{String, Nothing} = nothing)
    return IFU(30,
                "hexagonal",
                4200,
                5.65,
                1.,
                2.,
                filter)
end

"""
    Hector(;filter)

Creates an IFU using parameters of the CALIFA survey.
Optional filters "r" or "g" can also be specified for use with SSP particles.
"""
function Hector(;filter::Union{String, Nothing} = nothing)
    return IFU(30,
                "hexagonal",
                4800,
                1.3,
                1.6,
                0.05,
                filter)
end
