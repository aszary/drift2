module Lines
    using LinearAlgebra
    include("functions.jl")
    include("field.jl")
    using .Functions
    using .Field

    export generate_dipole!, generate_sphere

    function generate_sphere(r; size=10)
        u = range(0., stop=2*pi, length=size)
        v = range(0., stop=pi, length=size)
        x = r  * cos.(u) .* sin.(v)'
        y = r  * sin.(u) .* sin.(v)'
        z = r  * ones(size) .* cos.(v)'
        return [x, y, z]
    end


    """
    generate_dipole!(psr; size=100, z=50)

# Arguments

- z: distance from the star's center [in stellar radius]
    """
    function generate_dipole!(psr; size=100, phi_num=20, z=50)
        #phi_num = 10
        # last open magnetic field line
        rlc = psr.r_lc
        tm = theta_max(z, psr)
        thetas = range(0, tm, length=size)
        phis = range(0, 2*pi, length=phi_num)
        for ph in phis
            x = Array{Float64}(undef, size)
            y = Array{Float64}(undef, size)
            z = Array{Float64}(undef, size)
            for (i,th) in enumerate(thetas)
                r = rlc * sin(th)^2 # dipolar magnetic field?
                ca = Functions.spherical2cartesian([r, th, ph])
                x[i] = ca[1]
                y[i] = ca[2]
                z[i] = ca[3]
                #println(th)
            end
            push!(psr.lines, [x, y, z])
        end
    end

    function calculate_polarcap!(psr; phi_num=100)
        theta = Functions.theta_max(1, psr)
        phis = range(0, 2*pi, length=phi_num)
        x = Array{Float64}(undef, phi_num)
        y = Array{Float64}(undef, phi_num)
        z = Array{Float64}(undef, phi_num)
        for (i,ph) in enumerate(phis)
            ca = Functions.spherical2cartesian([psr.r, theta, ph])
            x[i] = ca[1]
            y[i] = ca[2]
            z[i] = ca[3]
        end
        psr.pc = [x, y, z]
    end

    """
    Generates magnetic and electric field lines for vacuum around neutron star

        step in meters
    """
    function generate_vacuum!(psr; step=10, stepsnum=20000, phi=nothing)
        fv = psr.field_vacuum

        # starting points
        r = psr.r
        #thetas = LinRange(0, pi/2, fv.size)
        #thetas = LinRange(0, pi, fv.size+1)[1:end-1]
        thetas = LinRange(0, pi, fv.size)
        if phi === nothing
            phis = LinRange(0, 2pi, fv.size+1)[1:end-1] # get rid of last point
        else
            phis = [phi, phi+pi]
        end

        # TODO add the second half! done? but too many lines?
        
        for i in 1:size(thetas)[1]
            for j in 1:size(phis)[1]
                pos_sph = [r, thetas[i], phis[j]]
                b_sph = Field.bvac(pos_sph, psr.r, fv.beq)
                e_sph = Field.evac(pos_sph, psr.r, fv.beq, psr.omega)
                pos = Functions.spherical2cartesian(pos_sph)                
                b = Functions.vec_spherical2cartesian(pos_sph, b_sph)
                e = Functions.vec_spherical2cartesian(pos_sph, e_sph)
                push!(fv.magnetic_lines, [[pos[1]], [pos[2]], [pos[3]]]) # adding initial position
                ml = fv.magnetic_lines[end]
                posb = copy(pos)
                step = abs(step) # start with positive step
                for k in 1:stepsnum
                    # new fields
                    posb_sph = Functions.cartesian2spherical(posb)
                    if posb_sph[1] < psr.r
                        # going the other direction if needed (e.g. southern hemisphere) or break
                        if size(ml[1], 1) > 2
                            #println("$i $j $k - break")
                            break
                        else
                            step = - step
                            #println("$i $j $k - minus")
                        end
                    end
                    #println(k)
                    b_sph = Field.bvac(posb_sph, psr.r, fv.beq)
                    b = Functions.vec_spherical2cartesian(posb_sph, b_sph)
                    st = b / norm(b) * step
                    posb += st # new position for magnetic
                    push!(ml[1], posb[1])
                    push!(ml[2], posb[2])
                    push!(ml[3], posb[3])
                end
                push!(fv.electric_lines, [[pos[1]], [pos[2]], [pos[3]]]) # adding initial position
                el = fv.electric_lines[end]
                pose = copy(pos)
                step = abs(step) # start with positive step
                for k in 1:stepsnum
                    pose_sph = Functions.cartesian2spherical(pose)
                    if pose_sph[1] < psr.r
                        # going the other direction if needed (e.g. southern hemisphere) or break
                        if size(el[1], 1) > 2
                            break
                        else
                            step = - step
                        end
                    end
                    e_sph = Field.evac(pose_sph, psr.r, fv.beq, psr.omega)
                    e = Functions.vec_spherical2cartesian(pose_sph, e_sph)
                    st = e / norm(e) * step
                    pose += st # new position for magnetic
                    push!(el[1], pose[1])
                    push!(el[2], pose[2])
                    push!(el[3], pose[3])
                end
            end
        end
    end

    """
    Generates magnetic and electric field lines for a force-free case

        step in meters
    """
    function generate_forcefree!(psr; step=10, stepsnum=20000, phi=nothing)
        ff = psr.field_forcefree

        # starting points
        r = psr.r

        thetas = LinRange(0.0001, pi-0.00001, ff.size) # E_ff = 0 for theta=0 and pi
        if phi === nothing
            phis = LinRange(0, 2pi, ff.size+1)[1:end-1] # get rid of last point
        else
            phis = [phi, phi+pi]
        end

        omega_vec = psr.omega_vec

        for i in 1:size(thetas)[1]
            for j in 1:size(phis)[1]
                pos_sph = [r, thetas[i], phis[j]]
                b_sph = Field.bvac(pos_sph, psr.r, ff.beq)
                b = Functions.vec_spherical2cartesian(pos_sph, b_sph)
                pos = Functions.spherical2cartesian(pos_sph)
                e = Field.eff(pos, omega_vec, b)
                push!(ff.magnetic_lines, [[pos[1]], [pos[2]], [pos[3]]]) # adding initial position
                ml = ff.magnetic_lines[end]
                posb = copy(pos)
                step = abs(step) # start with positive step
                for k in 1:stepsnum
                    # new fields
                    posb_sph = Functions.cartesian2spherical(posb)
                    if posb_sph[1] < psr.r
                        # going the other direction if needed (e.g. southern hemisphere) or break
                        if size(ml[1], 1) > 2
                            #println("$i $j $k - break")
                            break
                        else
                            step = - step
                            #println("$i $j $k - minus")
                        end
                    end
                    #println(k)
                    b_sph = Field.bvac(posb_sph, psr.r, ff.beq)
                    b = Functions.vec_spherical2cartesian(posb_sph, b_sph)
                    st = b / norm(b) * step
                    posb += st # new position for magnetic
                    push!(ml[1], posb[1])
                    push!(ml[2], posb[2])
                    push!(ml[3], posb[3])
                end
                push!(ff.electric_lines, [[pos[1]], [pos[2]], [pos[3]]]) # adding initial position
                el = ff.electric_lines[end]
                pose = copy(pos)
                step = abs(step) # start with positive step
                for k in 1:stepsnum
                    pose_sph = Functions.cartesian2spherical(pose)
                    if pose_sph[1] < psr.r
                        # going the other direction if needed (e.g. southern hemisphere) or break
                        if size(el[1], 1) > 2
                            break
                        else
                            step = - step
                        end
                    end
                    b_sph = Field.bvac(pose_sph, psr.r, ff.beq)
                    b = Functions.vec_spherical2cartesian(pose_sph, b_sph)
                    #println(pose)
                    e = Field.eff(pose, omega_vec, b)
                    #println(e)
                    st = e / norm(e) * step
                    pose += st # new position for magnetic
                    push!(el[1], pose[1])
                    push!(el[2], pose[2])
                    push!(el[3], pose[3])
                end
            end
        end
    end

end  # module Lines
