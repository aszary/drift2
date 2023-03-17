module Field
    using LinearAlgebra
    using PhysicalConstants.CODATA2018

    include("functions.jl")
    using .Functions


    mutable struct Vacuum
        size # number of points to calculate (size*size*size) or size*size for lines
        rmax # radius in meters # TODO change?
        locations # point locations
        magnetic # global star's magnetic field in vacuum
        electric # global electric field in vacuum
        beq # Magnetic field strength at the stellar equator
        magnetic_lines # magnetic field lines
        electric_lines # electric field lines
        locs2 # point location (ns interior)
        eint # electric field in the interior (only radial component?)
        function Vacuum(; size=20, rmax=50e3)
            return new(size, rmax, [], [], [], nothing, [], [], [], [])
        end
    end


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
    |B_d| = 2 at the pole
    """
    function bd(x, y, psr)
        d = sqrt(x ^ 2 + y ^ 2)
        theta = asin(d / psr.r)
        bd_sph = dipole(1, theta)
        #println(dipole(1, 0))
        bd = spherical2cartesian(bd_sph)
        return bd
    end


    """
    Magnetic field strength at the stellar equator (Handbook p. 267)
    """
    function beq(p, pdot)
        return 3.2e19 * sqrt(p * pdot)
    end


    """
    Spherical components of magnetic field for an aligned rotator (Cerutti, 2016 p.3) for Q=0
    """
    function bvac(pos_sph, rstar, beq)
        #println(beq)
        r = pos_sph[1]
        theta = pos_sph[2]
        phi = pos_sph[3]
        br = beq * (rstar / r) ^ 3 * 2 * cos(theta)
        btheta = beq * (rstar / r) ^ 3 * sin(theta)
        bphi = 0
        return(br, btheta, bphi)
    end


    """
    Spherical components of electric field for an aligned rotator (Cerutti, 2016 p.3) for Q=0
    """
    function evac(pos_sph, rstar, beq, omega)
        r = pos_sph[1]
        theta = pos_sph[2]
        phi = pos_sph[3]
        c = SpeedOfLightInVacuum.val # no units hereafter
        pre = omega * rstar / c * beq * (rstar / r) ^ 4
        er = pre * (1 - 3 * cos(theta)^2)
        etheta = pre * (- sin(2 * theta))
        ephi = 0
        return(er, etheta, ephi)
    end


    """
    Electric field caused by the interior charge Equation in text (between 3/4) in Cerutti (2017) ONLY Er?
    """
    function e_int(theta, rstar, beq, omega)
        c = SpeedOfLightInVacuum.val # no units hereafter
        mu = beq * rstar ^ 3 
        er_int = omega * mu * sin(theta)^2 / (c * rstar ^2) 
        return [er_int, 0, 0]
    end


    """
    Calculates magnetic and electric fields in vacuum (using Vacuum class).
    """
    function calculate_vac!(psr)
        fv = psr.field_vacuum
        fv.beq = beq(psr.p, psr.pdot)
        #println(fv)

        rs = LinRange(psr.r, fv.rmax, fv.size)
        thetas = LinRange(0, pi, fv.size)
        phis = LinRange(0, 2pi, fv.size)

        for i in 1:fv.size
            for j in 1:fv.size
                for k in 1:fv.size
                    #println(rs[i], " ", thetas[j], " ", phis[k])
                    pos_sph = [rs[i], thetas[j], phis[k]]
                    b_sph = bvac(pos_sph, psr.r, fv.beq)
                    e_sph = evac(pos_sph, psr.r, fv.beq, psr.omega)
                    push!(fv.locations, Functions.spherical2cartesian(pos_sph))
                    push!(fv.magnetic, Functions.vec_spherical2cartesian(pos_sph, b_sph))
                    push!(fv.electric, Functions.vec_spherical2cartesian(pos_sph, e_sph))
                end
            end
        end
    end


    """
    Calculates electric field in the interior of the neutron star
    """
    function calculate_eint!(psr)
        fv = psr.field_vacuum
        fv.beq = beq(psr.p, psr.pdot)

        rs = LinRange(0, psr.r, fv.size)
        thetas = LinRange(0, pi, fv.size)
        phis = LinRange(0, 2pi, fv.size)

        for i in 1:fv.size
            for j in 1:fv.size
                for k in 1:fv.size
                    #println(rs[i], " ", thetas[j], " ", phis[k])
                    pos_sph = [rs[i], thetas[j], phis[k]]
                    er = evac(pos_sph, psr.r, fv.beq, psr.omega)
                    push!(fv.locs2, Functions.spherical2cartesian(pos_sph))
                    push!(fv.eint, Functions.vec_spherical2cartesian(pos_sph, er))
                end
            end
        end

    end


end
