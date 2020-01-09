# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

function build_datacube(galaxy_data::Array{Galaxy_particle, 1};
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

    if(isnothing(blur))     #Without spatial blurring
        galaxy_data, parts_in_cell, ap_region, sbin, vbin = obs_data_prep(
                                                            galaxy_data, r200=r200, z=z, fov=fov,
                                                            ap_shape=ap_shape, central_wvl=central_wvl,
                                                            lsf_fwhm=lsf_fwhm, pixel_sscale=pixel_sscale,
                                                            pixel_vscale=pixel_vscale, inc_deg=inc_deg)

        fluxes = flux_grid(parts_in_cell, ap_region, sbin, vbin, z, filter)
        return fluxes
        #ifu_imgs = ifu_cube(galaxy_data, fluxes)

        #return ifu_imgs
    else()                  #With spatial blurring
        observe_data = obs_data_prep(galaxy_data, r200, z, fov, ap_shape, central_wvl,
                                        lsf_fwhm, pixel_sscale, pixel_vscale, inc_deg)

        fluxes = flux_grid(observe_data, filter)

        ifu_imgs = ifu_cube(observe_data, fluxes)
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
