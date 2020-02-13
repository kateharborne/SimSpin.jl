# Date created: 12/02/2020
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    Environment(z, inc_deg, r200, blur)

Creates a `struct` containing environmental parameters required for a mock observation of a simulated galaxy.

Parameters:\n
    z               The projected redshift at which the observation is made.
    inc_deg         The inclination at which to observe the galaxy in degrees.
    r200            The virial radius specified in the simulation, kpc.
    blur            Optional. Struct of type `Blur` containing seeing information. If ommitted no blurring is used.
"""
struct Environment

    z::Float64
    inc_deg::Real
    r200::Real
    blur::Union{Blur, Nothing}

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        blur::Blur)

        new(z, inc_deg, r200, blur)
    end

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real)

        blur = nothing
        new(z, inc_deg, r200, blur)
    end
end
