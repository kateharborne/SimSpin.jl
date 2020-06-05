# Date created: 10/02/2020
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    Observation(z, inc_deg, r200, blur, ap_region, sbin, vbin, vseq, lsf_size, ang_size, sbinsize, disc_mass2light, bulge_mass2light, lum_dist, redshift_coef)

Creates a `struct` containing all the parameters required for a mock observation of a simulated galaxy.

Parameters:\n
    z                   The projected redshift at which the observation is made.
    inc_deg             The inclination at which to observe the galaxy in degrees.
    r200                The virial radius specified in the simulation, kpc.
    blur                Struct of type `Blur`. Contains seeing information.
    ap_region           The aperture region mask used to remove flux outside of the specified aperture.
    sbin                The number of spatial bins in the aperture.
    vbin                The number of velocity bins in the flux grid.
    vseq                The bounds of each velocity bin in the flux grid.
    lsf_size            The Gaussian standard deviation of the line spread function in km/s.
    ang_size            The angular size given redshift z in kpc.
    sbinsize            The spatial bin size in kpc per bin.
    disc_mass2light     The mass to light ratio for disk particles.
    bulge_mass2light    The mass to light ratio for bulge particles.
    lum_dist            The luminosity distance in Mpc.
    redshift_coef       The redshift coefficient.
"""
struct Observation
    z::Float64
    inc_deg::Real
    r200::Real
    blur::Union{Blur, Nothing}
    ap_region::Array{Float64}
    sbin::Int64
    sseq::StepRangeLen
    vbin::Int64
    vseq::StepRangeLen
    lsf_size::Float64
    ang_size::Float64
    sbinsize::Float64
    disc_mass2light::Real
    bulge_mass2light::Real
    lum_dist::Real
    redshift_coef::Real

    function Observation(z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        blur::Union{Blur, Nothing},
                        ap_region::Array{Float64, 2},
                        sbin::Int64,
                        sseq::StepRangeLen,
                        vbin::Int64,
                        vseq::StepRangeLen,
                        lsf_size::Float64,
                        ang_size::Float64,
                        sbinsize::Float64,
                        disc_mass2light::Real,
                        bulge_mass2light,
                        lum_dist::Real,
                        redshift_coef::Real)

    new(z,
        inc_deg,
        r200,
        blur,
        ap_region,
        sbin,
        sseq,
        vbin,
        vseq,
        lsf_size,
        ang_size,
        sbinsize,
        disc_mass2light,
        bulge_mass2light,
        lum_dist,
        redshift_coef)
    end
end
