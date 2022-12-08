module Field
    using LinearAlgebra
    include("functions.jl")
    using .Functions

    """
    dipole(r, theta)

    Dipole magnetic field

    # Arguments

    - r: in stellar radius units
    - theta: in radians

    """
    function dipole(r, theta)
        b_r = 2. * cos(theta) / (r ^ 3)
        b_theta = sin(theta) / (r ^ 3)
        return (b_r, b_theta, 0.)
    end


    """
    dipolar component of magnetic field at the surface based on x, y components
    """
    function bd(x, y, psr)
        d = sqrt(x ^ 2 + y ^ 2)
        theta = asin(d / psr.r)
        bd_sph = dipole(1, theta)
        bd = spherical2cartesian(bd_sph)
        bd /= norm(bd)
        return bd
    end

end
