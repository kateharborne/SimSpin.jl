# Date created: 10/02/2020
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

abstract type Telescope end

"""
    IFU(fov, ap_shape, pixel_sscale, pixel_vscale, central_wvl, lsf_fwhm)

Creates a mock IFU telescope which can be used to \"observe\" simulations.

Parameters:\n
    fov             The field of view of the IFU, diameter in arcseconds.
    ap_shape        The shape of the field of view, with options "circular", "square" or "hexagonal".
    central_wvl     The central filter wavelength used for the observation, given in angstroms.
    lsf_fwhm        The line spread function full-width half-max, given in angstroms.
    pixel_sscale    The corresponding spatial pixel scale associated with a given telescope output in arcseconds.
    pixel_vscale    The corresponding velocity pixel scale associated with a given telescope filter output in angstroms.
"""
struct IFU <: Telescope
    fov::Float64
    ap_shape::String
    central_wvl::Float64
    lsf_fwhm::Float64
    pixel_sscale::Float64
    pixel_vscale::Float64

    function IFU(fov::Float64,
                        ap_shape::String,
                        central_wvl::Int64,
                        lsf_fwhm::Float64,
                        pixel_sscale::Float64,
                        pixel_vscale::Float64,)

        if (ap_shape != "circular" && ap_shape != "square" & ap_shape != "hexagonal")
            error("Please specify ap_shape as either 'circular', 'square' or 'hexagonal'.")
        end

        new(fov, ap_shape, central_wvl, lsf_fwhm, pixel_sscale, pixel_vscale)
    end
end



#function SAMI()
#    return IFU()
#end
