module Functions
    #using CoordinateTransformations # different convention
    using Unitful
    using PhysicalConstants.CODATA2018
    export spherical2cartesian, spherical2cartesian2, rdp, rlc, theta_max

    """
    spherical2cartesian(spherical)

Converts spherical cordinates to cartesian ones...
    """
    function spherical2cartesian(spherical)
        spherical
        x = spherical[1] * sin(spherical[2]) * cos(spherical[3])
        y = spherical[1] * sin(spherical[2]) * sin(spherical[3])
        z = spherical[1] * cos(spherical[2])
        return [x, y, z]
    end


    """
    spherical2cartesian2(spherical)

https://github.com/JuliaGeometry/CoordinateTransformations.jl/issues/25

Not used  [3, 0, 0] --> [3.0, 0.0, 0.0], but expecting [0.0, 0.0, 3.0] (along z axis - wiki ISO definition, not math wolfram...)
    """
    function spherical2cartesian2julianotation(spherical)
        return CartesianFromSpherical()(Spherical(spherical[1], spherical[2], spherical[3]))
    end


    """
    r_dp(p, r)

Polar cap radius assuming purely dipolar configuration of the magnetic field at the surface

# Arguments

- p: pulsar period
- r: neutron star radius

returns radius of the polar in meters
    """
    function rdp(p, r)
        #https://juliaphysics.github.io/PhysicalConstants.jl/stable/reference/
        c = SpeedOfLightInVacuum.val # no units hereafter
        #println((2. * pi * r ^ 3. / (SpeedOfLightInVacuum * p)) ^ 0.5)
        return (2. * pi * r ^ 3. / (c * p)) ^ 0.5
    end


    """
    theta_max(z, psr)

Eq. 3.23 in the handbook

Note this is not an opening angle! (Note 1.5 differecence E.g 3.29 in the handbook)

# Arguments

- z: distance from the star's center [in stellar radius]
- psr: pulsar class (struct) - to get period and star r

returns the theta component of the last open magnetic field line for a given distance from the star's
    center and pulsar period [in radians]
    """
    function theta_max(z, psr)
        return asin(sqrt(z * psr.r / rlc(psr.p)))
    end


    """
    rlc(p)

Calculates radius of the light cylinder [in meters]

# Arguments

- p: pulsar period
    """
    function rlc(p)
        c = SpeedOfLightInVacuum.val # no units hereafter
        return c * p / (2 * pi)
    end


end
