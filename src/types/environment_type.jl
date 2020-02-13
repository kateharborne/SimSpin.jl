# Date created: 12/02/2020
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    Environment(z, inc_deg, r200, mass2light, blur)

Creates a `struct` containing environmental parameters required for a mock observation of a simulated galaxy.

Parameters:\n
    z               The projected redshift at which the observation is made.
    inc_deg         The inclination at which to observe the galaxy in degrees.
    r200            The virial radius specified in the simulation, kpc.
    mass2light      Optional. The mass to light ratio for non-ssp, luminous particles. Defaults to 1.
    blur            Optional. Struct of type `Blur` containing seeing information. If ommitted no blurring is used.
"""
struct Environment

    z::Float64
    inc_deg::Real
    mass2light::Float64
    r200::Real
    blur::Union{Blur, Nothing}

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        mass2light::Float64,
                        blur::Blur)

        new(z, inc_deg, r200, mass2light, blur)
    end

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        blur::Blur)

        mass2light=1.
        new(z, inc_deg, r200, mass2light, blur)
    end

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        mass2light::Float64)

        blur = nothing
        new(z, inc_deg, r200, mass2light, blur)
    end

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real)
        mass2light=1.
        blur = nothing
        new(z, inc_deg, r200, blur)
    end
end
