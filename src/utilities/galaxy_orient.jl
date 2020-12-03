# Date created: 10/01/2020
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

using LinearAlgebra

"""
    galaxy_orient(sim_data, centre)

Returns rotation matrix to convert from simulation reference frame to
the galaxy reference frame where galaxy's angular momentum vector is aligned with x.
"""
function galaxy_orient(sim_data::Array{Sim_particle, 1}, centre::Array{Float64, 1})

    J = galaxy_angmom(sim_data, centre)
    J_norm = J ./ sqrt(sum(J.^2))

    if J_norm[3] == -1
        rot_mat =  [1 0 0;
                    0 -1 0;
                    0 0 -1]
        return rot_mat
    elseif J_norm[3] == 1
        return I(3)
    end

    v = [J_norm[2], -J_norm[1], 0.] # unit vector normal to J and z-axis, about which we want to rotate
    c = J_norm[3]                   # cos(angle)
    s = sqrt(sum(v.^2))             # sin(angle)

    v_x =  [0. -v[3] v[2];          # skew-symmetric cross product
            v[3] 0. -v[1];
            -v[2] v[1] 0.]

    rot_mat = I(3) + v_x + (1/(1+c)) * (v_x * v_x) # rotation matrix via Rodrigues Rotation Formula:
    println(rot_mat)                                            # wikipedia.org/wiki/Rodrigues'_rotation_formula
    return rot_mat
end

function galaxy_angmom(sim_data::Array{Sim_particle, 1}, centre::Array{Float64, 1})

    x = getfield.(sim_data, :x) .- centre[1]
    y = getfield.(sim_data, :y) .- centre[2]
    z = getfield.(sim_data, :z) .- centre[3]

    radii = @. sqrt(x^2 + y^2 + z^2)
    close_part = radii .< 0.33 * maximum(radii)

    x = x[close_part]
    y = y[close_part]
    z = z[close_part]

    vx = getfield.(sim_data[close_part], :vx)
    vy = getfield.(sim_data[close_part], :vy)
    vz = getfield.(sim_data[close_part], :vz)
    mass = getfield.(sim_data[close_part], :mass)

    Jx = sum(@. mass * ((y * vz) - (z * vy)))
    Jy = sum(@. mass * ((x * vz) - (z * vx)))
    Jz = sum(@. mass * ((x * vy) - (y * vx)))

    return [Jx, Jy, Jz]
end
