# Date created: 10/02/2020
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne


"""
    Observation(z, inc_deg, r200)

Creates a 'struct' specifying the setup of a galaxy for observation. To be used in conjunction with a 'Telescope'.

Parameters:\n
    z               The projected redshift at which the observation is made.
    inc_deg         The inclination at which to observe the galaxy in degrees.
    r200            The virial radius specified in the simulation, kpc.
"""
struct Observation
    z::Float64
    inc_deg::Float64
    r200::Int64

    function Observation(z::Float64,
                        inc_deg::Float64,
                        r200::Int64)
        new(z, inc_deg, r200)
    end
end
