# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

using LinearAlgebra

"""
    galaxy_orient(sim_data, centre)

Returns rotation matrix to convert from simluation reference frame to
the galaxy reference frame where x is aligned with the semi-major axis.
"""
function galaxy_orient(sim_data::Array{Sim_particle, 1}, centre::Array{Float64, 1})

    x = getfield.(sim_data, :x) .- centre[1]
    y = getfield.(sim_data, :y) .- centre[2]
    z = getfield.(sim_data, :z) .- centre[3]

    mass = getfield.(sim_data, :mass)

    # step 1: define the moment of inertia tensor
    inertiaTensor = randn(3,3)
    inertiaTensor[1,1] = sum(mass .* (y.^2 .+ z.^2))
    inertiaTensor[2,2] = sum(mass .* (x.^2 .+ z.^2))
    inertiaTensor[3,3] = sum(mass .* (x.^2 .+ y.^2))
    inertiaTensor[1,2] = -sum(mass .* x .* y)
    inertiaTensor[1,3] = -sum(mass .* x .* z)
    inertiaTensor[2,3] = -sum(mass .* y .* z)
    inertiaTensor[2,1] = inertiaTensor[1,2]
    inertiaTensor[3,1] = inertiaTensor[1,3]
    inertiaTensor[3,2] = inertiaTensor[2,3]

    # step 2: find eigen vectors and reorder such that x is the major axis
    eigen_vec = eigvecs(inertiaTensor)
    rot_mat = eigen_vec' .* [-1 -1 -1;1 1 1;-1 -1 -1]   # tranpose to get rotation matrix

    return rot_mat
end
