# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    build_datacube(galaxy_data;
                    r200 = 200,
                    z = 0.05,
                    fov = 15.,
                    ap_shape = "circular",
                    central_wvl = 4800,
                    lsf_fwhm = 2.65,
                    pixel_sscale = 0.5,
                    pixel_vscale = 1.04,
                    inc_deg = 70,
                    filter = "r",
                    blur = nothing)

Returns a simulated ifu datacube for input, `galaxy_data`, an array of `Particle`.
Observational parameters specified as keyword arguments.

Keyword arguments (optional):\n
    r200            The virial radius specified in the simulation, kpc.
    z               The projected redshift at which the observation is made.
    fov             The field of view of the IFU, diameter in arcseconds.
    ap_shape        The shape of the field of view, with options "circular", "square" or "hexagonal".
    central_wvl     The central filter wavelength used for the observation, given in angstroms.
    lsf_fwhm        The line spread function full-width half-max, given in angstroms.
    pixel_sscale    The corresponding spatial pixel scale associated with a given telescope output in arcseconds.
    pixel_vscale    The corresponding velocity pixel scale associated with a given telescope filter output in angstroms.
    inc_deg         The inclination at which to observe the galaxy in degrees.
    filter          If particles type is ssp, the filter within which the SED is generated. Options include "r" and "g"  for SDSS-r and SDSS-g bands respectively.
    blur            Specified to apply observational seeing effects to the cube. A list of the form "psf" = "Moffat", "fwhm" = 0.5. "psf" specifies the shape of the PSF chosen and may be either "Moffat" or "Gaussian". "fwhm" is a numeric specifying the full-width half-maximum of the PSF given in units of arcseconds.
"""
function build_datacube(galaxy_data::Array{Galaxy_particle, 1};
                        r200::Int64 = 200,
                        z::Float64=0.05,
                        fov::Float64=15.,
                        ap_shape::String="circular",
                        central_wvl::Int64=4800,
                        lsf_fwhm::Float64=2.65,
                        pixel_sscale::Float64=0.5,
                        pixel_vscale::Float64=1.04,
                        inc_deg::Int64=70,
                        filter::String="r",
                        blur=nothing)

    galaxy_data,
    parts_in_cell,
    ap_region, sbin,
    vbin, vseq, lsf_size = obs_data_prep(galaxy_data, r200=r200, z=z, fov=fov,
                                    ap_shape=ap_shape, central_wvl=central_wvl,
                                    lsf_fwhm=lsf_fwhm, pixel_sscale=pixel_sscale,
                                    pixel_vscale=pixel_vscale, inc_deg=inc_deg)

    fluxes = flux_grid(parts_in_cell, ap_region, sbin, vbin, z, filter)
    ifu_imgs = ifu_cube(fluxes, parts_in_cell, sbin, vbin, vseq, lsf_size)

    if(isnothing(blur)) #No spatial blurring
        return ifu_imgs
    else()
        error("Image blur is not currently supported.")
        blur_imgs = blur_cube(observe_data, ifu_imgs, blur.psf, blur.fwhm)

        return blur_imgs
    end
end


function build_datacube(sim_data::Array{Sim_particle, 1};
                        r200::Int64 = 200,
                        z::Float64=0.1,
                        fov::Float64=15.,
                        ap_shape::String="circular",
                        central_wvl::Int64=4800,
                        lsf_fwhm::Float64=2.65,
                        pixel_sscale::Float64=0.5,
                        pixel_vscale::Float64=1.04,
                        inc_deg::Int64=0,
                        filter::String="r",
                        blur=nothing)

galaxy_data = sim_to_galaxy(sim_data)

return build_datacube(galaxy_data, r200=r200, z=z, fov=fov,
                        ap_shape=ap_shape, central_wvl=central_wvl,
                        lsf_fwhm=lsf_fwhm,
                        pixel_sscale=pixel_sscale,
                        pixel_vscale=pixel_vscale,
                        inc_deg=inc_deg,
                        filter=filter,
                        blur=blur)
end
