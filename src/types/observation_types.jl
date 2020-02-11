# Date created: 10/02/2020
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    Observation(z, inc_deg, r200, ap_region, sbin, vbin, vseq, lsf_size, ang_size, sbinsize)

Creates a `struct` containing all the parameters required for a mock observation of a simulated galaxy.

Parameters:\n
    z               The projected redshift at which the observation is made.
    inc_deg         The inclination at which to observe the galaxy in degrees.
    r200            The virial radius specified in the simulation, kpc.
    ap_region       The aperture region mask used to remove flux outside of the specified aperture.
    sbin            The number of spatial bins in the aperture.
    vbin            The number of velocity bins in the flux grid.
    vseq            The bounds of each velocity bin in the flux grid.
    lsf_size        The Gaussian standard deviation of the line spread function in km/s.
    ang_size        The angular size given redshift z in kpc
    sbinsize        The spatial bin size in kpc per bin
"""
struct Observation
    z::Float64
    inc_deg::Real
    r200::Real
    ap_region::Array{Float64}
    sbin::Int64
    vbin::Int64
    vseq::StepRangeLen
    lsf_size::Float64
    ang_size::Float64
    sbinsize::Float64

    function Observation(z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        ap_region::Array{Float64, 2},
                        sbin::Int64,
                        vbin::Int64,
                        vseq::StepRangeLen,
                        lsf_size::Float64,
                        ang_size::Float64,
                        sbinsize::Float64)

    new(z, inc_deg, r200, ap_region, sbin, vbin, vseq, lsf_size, ang_size, sbinsize)
    end
end
