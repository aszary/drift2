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
        locs # point location (ns interior)
        eint # electric field in the interior (only radial component?)
        locs2 # point location 2 
        eint2 # internal electric field
        locs3 # point location 3 
        gj3 # Goldreich-Julian density
        xgj # x coordinate for the heatmap plot
        zgj # z coordinate for the heatmap plot
        gj # Goldreich-Julian density fir the heatmap plot
        function Vacuum(; size=16, rmax=50e3)
            return new(size, rmax, [], [], [], nothing, [], [], [], [], [], [], [], [], [], [], [])
        end
    end


    mutable struct ForceFree 
        size # number of points to calculate (size*size*size) or size*size for lines
        rmax # radius in meters
        locations # point locations
        magnetic # global star's magnetic field in vacuum
        electric # global electric field in vacuum
        beq # Magnetic field strength at the stellar equator
        magnetic_lines # magnetic field lines
        electric_lines # electric field lines
        locs3 # point location 3 
        gj3 # Goldreich-Julian density
        xgj # x coordinate for the heatmap plot
        zgj # z coordinate for the heatmap plot
        gj # Goldreich-Julian density fir the heatmap plot
        function ForceFree(; size=16, rmax=50e3)
            return new(size, rmax, [], [], [], nothing, [], [], [], [], [], [], [])
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
        #println(rstar, r, beq)
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
    function e_int(pos_sph, rstar, beq, omega)
        theta = pos_sph[1]
        c = SpeedOfLightInVacuum.val # no units hereafter
        mu = beq * rstar ^ 3 
        er_int = omega * mu * sin(theta)^2 / (c * rstar ^2) 
        return [er_int, 0, 0]
    end

    """
    Internal electric field Equation 1 in Cerutti (2017) in cartesian coordinates !
    """
    function e_int2(pos_sph, psr)
        fv = psr.field_vacuum
        omega_vec = psr.omega_vec
        r = Functions.spherical2cartesian(pos_sph)
        v = cross(omega_vec, r)
        b = Functions.vec_spherical2cartesian(pos_sph, bvac(pos_sph, psr.r, fv.beq))
        c = SpeedOfLightInVacuum.val # no units hereafter
        eint2 = -cross(v, b) / c
        #println(eint2)
        return eint2
    end


    """
    Calculates Goldreich-Julian charge density
    """
    function GJ_density(pos_sph, r, omega_vec, beq)
        #fv = psr.field_vacuum
        #omega_vec = psr.omega_vec
        b = Functions.vec_spherical2cartesian(pos_sph, bvac(pos_sph, r, beq))
        c = SpeedOfLightInVacuum.val * 1e2 # in cm # no units hereafter
        gj = -dot(omega_vec, b) / (2 * pi * c)
        return gj
    end


    """
    Calculates Goldreich-Julian number charge density
        # TODO change to SI particle / cm^3
        n = (B * Omega) / (2 * pi * c * qe)
        qe = 4.8e-10  # Ładunek elementarny w jednostkach elektrostatycznych w statcoulombs
    """
    function GJ_ndensity(pos_sph, r, omega_vec, beq)
        #fv = psr.field_vacuum
        #omega_vec = psr.omega_vec
        b = Functions.vec_spherical2cartesian(pos_sph, bvac(pos_sph, r, beq))
        c = SpeedOfLightInVacuum.val * 1e2 # in cm no units hereafter
        qe = 4.8032047e-10 # statcoulombs
        gj = -dot(omega_vec, b) / (2 * pi * c * qe)
        return gj
    end




    """
    Electric field in force-free magnetosphere Ruderman & Sutherland 1975 Eq. 2

    all atributes in cartesian coordinates
    r, omega_vec, B - magnetic field

    returns cartesian components of electric field 
    """
    function eff(r, omega_vec, B)
        #pos_sph, rstar, beq, omega)
        c = SpeedOfLightInVacuum.val # no units hereafter
        E = -cross(cross(omega_vec, r), B) / c
        return E
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
    Calculates magnetic and electric fields in force-free (using ForceFree class).
    """
    function calculate_ff!(psr)
        ff = psr.field_forcefree
        ff.beq = beq(psr.p, psr.pdot)

        rs = LinRange(psr.r, ff.rmax, ff.size)
        thetas = LinRange(0, pi, ff.size)
        phis = LinRange(0, 2pi, ff.size)

        for i in 1:ff.size
            for j in 1:ff.size
                for k in 1:ff.size
                    #println(rs[i], " ", thetas[j], " ", phis[k])
                    pos_sph = [rs[i], thetas[j], phis[k]]
                    pos_car = Functions.spherical2cartesian(pos_sph)
                    b_sph = bvac(pos_sph, psr.r, ff.beq)
                    b_car = Functions.vec_spherical2cartesian(pos_sph, b_sph)
                    e_car = eff(pos_car, psr.omega_vec, b_car)
                    push!(ff.locations, pos_car)
                    push!(ff.magnetic, b_car)
                    push!(ff.electric, e_car)
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
                    eint = e_int(pos_sph, psr.r, fv.beq, psr.omega)
                    push!(fv.locs, Functions.spherical2cartesian(pos_sph))
                    push!(fv.eint, Functions.vec_spherical2cartesian(pos_sph, eint))
                end
            end
        end

    end

    """
    Calculates internal electric field of the neutron star 
    """
    function calculate_eint2!(psr)
        fv = psr.field_vacuum
        fv.beq = beq(psr.p, psr.pdot)

        r = psr.r
        thetas = LinRange(0, pi, fv.size)
        phis = LinRange(0, 2pi, fv.size)

        for j in 1:fv.size
            for k in 1:fv.size
                pos_sph = [r, thetas[j], phis[k]]
                eint2 = e_int2(pos_sph, psr)
                #println(typeof(eint2), eint2, pos_sph)
                push!(fv.locs2, Functions.spherical2cartesian(pos_sph))
                push!(fv.eint2, eint2)
            end
        end
    end


    """
    Calculates GJ charge density for the stellar surface 
    """
    function calculate_GJ!(psr; field=nothing, twoD=false)
        if field === nothing
            fv = psr.field_vacuum
        else
            fv = field
        end
        fv.beq = beq(psr.p, psr.pdot)

        r = psr.r
        thetas = LinRange(0.01, pi, fv.size)
        if twoD == false
            phis = LinRange(0, 2pi, fv.size)
        else
            phis = LinRange(0, pi, convert(Int, fv.size / 2))
        end

        for j in 1:size(thetas, 1)
            for k in 1:size(phis, 1)
                pos_sph = [r, thetas[j], phis[k]]
                gj = GJ_ndensity(pos_sph, psr.r, psr.omega_vec, fv.beq)
                #println(gj)
                push!(fv.locs3, Functions.spherical2cartesian(pos_sph))
                push!(fv.gj3, gj)
            end
        end
    end


    """
    Calculates GJ charge density for the heatmap 
    xlims and zlims in kilometers
    """
    function calculate_GJheat!(psr, field_class; xlims=(-30, 30), zlims=(-30, 30), vacuum=true, size=100)
        fi = field_class
        fi.beq = beq(psr.p, psr.pdot)

        xs = LinRange(xlims[1], xlims[2], size)
        zs = LinRange(zlims[1], zlims[2], size)

        for i in 1:size
            for j in 1:size
                x = xs[i] * 1e3 # in meters
                z = zs[j] * 1e3 # in meters
                # TODO fix this
                if (x^2 + z ^2 <= psr.r^2)
                    y = sqrt(psr.r^2 - x^2 - z^2)
                else
                    y = 0
                end
                #y = psr.r

                pos_sph = Functions.cartesian2spherical([x, y, z])
                gj = GJ_ndensity(pos_sph, psr.r, psr.omega_vec, fi.beq)
                #println("$(x/1e3) $(y/1e3) $gj")
                #println(gj)
                if vacuum == true
                    if (pos_sph[1] <= psr.r) #  && (gj < 1e3)
                        push!(fi.xgj, x)
                        push!(fi.zgj, z)
                        push!(fi.gj, gj)
                    end
                else     
                    push!(fi.xgj, x)
                    push!(fi.zgj, z)
                    push!(fi.gj, gj)
                end
            end
        end
    end

end
