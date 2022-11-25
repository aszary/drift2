module Lines
    include("functions.jl")
    using .Functions

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
                r = rlc * sin(th)^2
                ca = spherical2cartesian([r, th, ph])
                x[i] = ca[1]
                y[i] = ca[2]
                z[i] = ca[3]
                #println(th)
            end
            push!(psr.lines, [x, y, z])
        end
    end


end  # module Lines
