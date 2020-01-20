# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

function get_r(x, y, z)

    r = sqrt(x^2 + y^2 + z^2)
    return r
end

function get_cr(x, y)

    cr = sqrt(x^2 + y^2)
    return cr
end

function get_spherical(x, y, z, r)

    theta = acos(z/r)
    phi = atan(y/x)
    return theta, phi
end

function get_vr(x, y, z, vx, vy, vz)

    vr = (x*vx) + (y * vy) + (z * vz)
    return vr
end

function get_vt(x, y, z, vx, vy, vz, r)

    vt = sqrt( (y * vz - z * vy)^2 +
                (z * vx - x * vz)^2 +
                (x * vy - y * vx)^2) / r
    return vt
end

function get_vcr(x, y, vx, vy, cr)

    vcr = (x * vx + y *vy) / cr
    return vcr
end

function get_vtheta(z, vz, r, vr)

    vtheta = -(r / sqrt(1 - (z^2 / r^2))) *
                ((vz / r) - ((vr * z) / r^2))
    return vtheta
end

function get_vphi(x, y, vx, vy, r, theta)

    vphi = (r / ((y^2 / x^2) + 1)) * sin(theta) *
            ((vy/x) - vx * y / x^2)
    return vphi
end

function get_j(x, y, z, vx, vy, vz, mass)

    jx = ((y * vz) - (z * vy)) * mass * pc_to_m
    jy = ((z * vx) - (x * vz)) * mass * pc_to_m
    jz = ((x * vy) - (y * vx)) * mass * pc_to_m

    return jx, jy, jz
end
