module Field

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

end
