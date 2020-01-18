abstract type Shell end

struct Spherical_3d <: Shell
    rad::Float64
    mass::Float64
    logp::Float64
    vc::Float64
    Jx::Float64
    Jy::Float64
    Jz::Float64
    J::Float64
    vr::Float64
    sigma_vr::Float64
    simga_vt::Float64
    sigma2_vx::Float64
    sigma2_vz::Float64
    b::Float64
    vrot::Float64
    lambda::Float64

    function Spherical_3d(rad::Float64
                            mass::Float64
                            logp::Float64
                            vc::Float64
                            Jx::Float64
                            Jy::Float64
                            Jz::Float64
                            J::Float64
                            vr::Float64
                            sigma_vr::Float64
                            simga_vt::Float64
                            sigma2_vx::Float64
                            sigma2_vz::Float64
                            b::Float64
                            vrot::Float64
                            lambda::Float64)
        new(rad, mass, logp, vc,
            Jx, Jy, Jz, J,
            vr, sigma_vr, sigma_vt,
            sigma2_vx, sigma2_vz,
            b, vrot, lamda)
    end
end

struct Circ_radial_2d <: Shell

    cr::Float64
    mass::Float64
    logp::Float64
    Jx::Float64
    Jy::Float64
    Jz::Float64
    J::Float64
    vcr::Float64
    sigma_vcr::Float64

    function Circ_radial_2d(cr::Float64
                            mass::Float64
                            logp::Float64
                            Jx::Float64
                            Jy::Float64
                            Jz::Float64
                            J::Float64
                            vcr::Float64
                            sigma_vcr::Float64)
        new(cr, mass, logp,
            Jx, Jy, Jz, J,
            vcr, sigma_vcr)
    end
end

struct Circ_planar_2d <: Shell

    z::Float64
    mass::Float64
    logp::Float64
    Jx::Float64
    Jy::Float64
    Jz::Float64
    J::Float64
    vz::Float64
    sigma_vz::Float64

    function Circ_planar_2d(z::Float64
                            mass::Float64
                            logp::Float64
                            Jx::Float64
                            Jy::Float64
                            Jz::Float64
                            J::Float64
                            vz::Float64
                            sigma_vz::Float64)
        new(z, mass, logp,
            Jx, Jy, Jz, J,
            vz, sigma_vz)
    end
end
