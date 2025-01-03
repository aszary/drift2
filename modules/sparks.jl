module Sparks
    using LinearAlgebra
    using PyCall
    using CoupledFields
    include("functions.jl")
    include("field.jl")
    using .Functions
    using .Field


    function create_grid_sph!(psr, size=200)
        r = psr.r # stellar radius in meters
        theta_max = Functions.theta_max(1, psr)
        thetas = range(0, theta_max, length=size)
        phis = range(0, 2*pi, length=size)
        x = Array{Float64}(undef, size, size)
        y = Array{Float64}(undef, size, size)
        z = Array{Float64}(undef, size, size)
        for (i, th) in enumerate(thetas)
            for (j, ph) in enumerate(phis)
                car = Functions.spherical2cartesian([r, th, ph])
                x[i,j], y[i, j], z[i, j] =  car[1], car[2], car[3]
            end
        end
        psr.grid = [x, y, z]
    end



    """
    # https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
    """
    function create_grid_fib!(psr, size=1000)
        r = psr.r
        theta_max = Functions.theta_max(1, psr)
        rp = r * cos(theta_max) # calculation only at the polar cap
        phi = pi * (3. - sqrt(5.))  # golden angle in radians
        x = []
        y = []
        z = []
        for i in 1:size
            zz = r - ((i-1) / (size-1)) * (r-rp) #  # 2 * r  # z goes from r to -r
            #println(zz)
            radius = sqrt(r^2 - zz * zz)  # radius at z
            theta = phi * i  # golden angle increment
            xx = cos(theta) * radius
            yy = sin(theta) * radius
            push!(x, xx)
            push!(y, yy)
            push!(z, zz)
        end
        x = convert(Array{Float64,1}, x)
        y = convert(Array{Float64,1}, y)
        z = convert(Array{Float64,1}, z)
        psr.grid = [x, y, z]
    end

    """
    simple square grid, not applicable for gradient calculation? or..
    """
    function create_grid_square!(psr, size_=100)
        r = psr.r # stellar radius in meters
        x_min, x_max = extrema(psr.pc[1])
        y_min, y_max = extrema(psr.pc[2])
        z_min = psr.pc[3][1]

        theta_max = Functions.theta_max(1, psr)
        xs = range(x_min, x_max, length=size_)
        ys = range(y_min, y_max, length=size_)
        x = []
        y = []
        z = []
        for (i, xx) in enumerate(xs)
            for (j, yy) in enumerate(ys)
                zz = sqrt(r^2 - xx^2 - yy^2)
                #rr = sqrt(xx^2 + yy^2 + zz^2)
                if zz >= z_min
                    push!(x, xx)
                    push!(y, yy)
                    push!(z, zz)
                end
            end
        end
        x = convert(Array{Float64,1}, x)
        y = convert(Array{Float64,1}, y)
        z = convert(Array{Float64,1}, z)
        psr.grid = [x, y, z]
        psr.grid2 = reduce(hcat, [x, y, z]) # used to calculate gradient... https://stackoverflow.com/questions/70201776/converting-julia-nested-list-to-multidimensional-array
        #println(psr.grid[2][10])
        #println(psr.grid2[10, 2])
        #println(typeof(psr.grid2))
        #println(size(psr.grid2))
    end


    """

    Random sparks for old grid shape (x, y, z)

# Arguments

- min_dist: minimum distant in meters

    """
    function random_sparks_old!(psr; min_dist=30, trials=1000)
        gr = psr.grid
        grid_size = size(gr[1])[1]
        sp = []
        for i in 1:trials
            j = rand(1:grid_size)
            x = gr[1][j]
            y = gr[2][j]
            z = gr[3][j]
            if !([x, y, z] in sp) # not in sparks
                md = 2 * min_dist
                # check distance between sparks
                for s in sp
                    dist = norm(s - [x, y, z])
                    #println(dist)
                    if dist < md
                        md = dist
                    end
                end
                for i in 1:size(psr.pc[1])[1]
                    dist = norm([x, y, z] - [psr.pc[1][i], psr.pc[2][i], psr.pc[3][i]])
                    if dist < md
                        md = dist
                    end
                end
                if md > min_dist
                    push!(sp, [x, y, z])
                end
            end
        end
        #psr.sparks = convert(Array{Float64,1}, sp)
        psr.sparks = sp
        println("Number of sparks added: ", size(sp)[1])
    end

    """
    Electric potential [Filaments]

# Arguments

-r: distance from the spark forming region
    """
    function v(r; a=1)
        return a * log(r)
    end


    function calculate_potential_old!(psr)
        gr = psr.grid
        grid_size = size(gr[1])[1]
        sp = psr.sparks
        spark_num = size(sp)[1]

        vs = Array{Float64}(undef, grid_size)

        for i in 1:grid_size
            vv = 0
            for j in 1:spark_num
                if j != i
                    dist = norm([gr[1][i], gr[2][i], gr[3][i]] - [sp[j][1], sp[j][2], sp[j][3]])
                    if dist != 0
                        vv += v(dist)
                    end
                end
            end
            vs[i] = vv
            #println(vv)
        end

        # gradient does not work!
        #kernelpars = GaussianKP(psr.grid2) # does not work (grid too big)
        #kernelpars = PolynomialKP(psr.grid2) # does not work (grid too big)
        #gradvecfield([0.5 -7], X, Y[:,1:1], kernelpars)
        #grad_v = gradvecfield([0.5 -7], psr.grid2, vs[:,1:1], kernelpars)
        #println(size(grad_v))

        psr.potential = convert(Array{Float64,1}, vs)

    end



    """
    simple square grid, but different data size => applicable for gradient calculation

    # try this https://stackoverflow.com/questions/40338386/calculating-a-3d-gradient-with-unevenly-spaced-points next time?
    """
    function create_grid!(psr; size=50)
        r = psr.r # stellar radius in meters

        x_min, x_max = extrema(psr.pc[1])
        y_min, y_max = extrema(psr.pc[2])
        z_min = psr.pc[3][1]

        xs = range(1.1*x_min, 1.1*x_max, length=size)
        ys = range(1.1*y_min, 1.1*y_max, length=size)
        gr_x = xs
        gr_y = ys
        gr_z = zeros((size, size))
        for (i, xx) in enumerate(xs)
            for (j, yy) in enumerate(ys)
                zz = sqrt(r^2 - xx^2 - yy^2)
                if zz >= z_min
                    gr_z[i, j] = zz
                else
                    gr_z[i, j] = 0  # zero used as magic number if point is outsied of polar cap
                end
            end
        end
        psr.grid = [gr_x, gr_y, gr_z]
    end


    """

    Random sparks for new grid shape (x[size], y[size], z[size, size])
    sparks only by index [i, j] -> location is in grid
# Arguments

- min_dist: minimum distant in meters

    """
    function random_sparks_grid!(psr; min_dist=30, trials=1000)
        gr = psr.grid
        grid_size = size(gr[1])[1]
        sp = []
        for i in 1:trials
            ii = rand(1:grid_size)
            jj = rand(1:grid_size)
            x = gr[1][ii]
            y = gr[2][jj]
            z = gr[3][ii, jj]
            if !([ii, jj] in sp) && (z !=0) # not in sparks # remember to skip z=0
                md = 2 * min_dist
                # check distance between sparks
                for si in sp
                    sx = gr[1][si[1]]
                    sy = gr[2][si[2]]
                    sz = gr[3][si[1], si[2]]
                    dist = norm([sx, sy, sz] - [x, y, z])
                    #println(dist)
                    if dist < md
                        md = dist
                    end
                end
                for i in 1:size(psr.pc[1])[1]
                    dist = norm([x, y, z] - [psr.pc[1][i], psr.pc[2][i], psr.pc[3][i]])
                    if dist < md
                        md = dist
                    end
                end
                if md > min_dist
                    push!(sp, [ii, jj])
                end
            end
        end
        #psr.sparks = convert(Array{Float64,1}, sp)
        psr.sparks = sp
        println("Number of sparks added: ", size(sp)[1])
    end


    """
    Calculates electric potential, electric field and drift velocity for sparks defined by grid indexes...
    """
    function calculate_potential!(psr)
        gr = psr.grid
        grid_size = size(gr[1])[1]
        sp = psr.sparks
        spark_num = size(sp)[1]

        vs = Array{Float64}(undef, grid_size, grid_size)

        v_min =  1e50
        v_max =  -1e50

        for i in 1:grid_size
            for j in 1:grid_size
                vv = 0
                for k in 1:spark_num
                    (ii, jj) = sp[k]
                    if (gr[3][i, j]!=0) # (ii != i) && (jj !=j) && (gr[3][i, j]!=0) # why?
                        sx = gr[1][ii]
                        sy = gr[2][jj]
                        sz = gr[3][ii, jj]
                        dist = norm([gr[1][i], gr[2][j], gr[3][i, j]] - [sx, sy, sz])
                        #vv += v(dist) # nice looking dots (Inf) in the plot, but no
                        if dist != 0
                            vv += v(dist)
                        end
                    end
                end
                if gr[3][i, j] != 0
                    vs[i, j] = vv
                else
                    vs[i, j] = 0
                end
            end
        end

        # set 0 potential (beyound the polar cap) to NaN
        for i in 1:grid_size
            for j in 1:grid_size
                if vs[i, j] == 0
                    vs[i, j] = NaN
                end
            end
        end

        # calculate electric field
        ex = Array{Float64}(undef, grid_size, grid_size)
        ey = Array{Float64}(undef, grid_size, grid_size)

        # julia using diff (bad results!)
        """
        grad_vx = diff(vs, dims=1) # 199x200
        grad_vy = diff(vs, dims=2) # 200x199

        for i in 1:grid_size-1
            for j in 1:grid_size
                ex[i, j] = -grad_vx[i, j]
            end
        end
        # dirty hacks below
        for j in 1:grid_size
            ex[grid_size, j] = -grad_vx[grid_size-1, j]
        end

        for i in 1:grid_size
            for j in 1:grid_size-1
                ey[i, j] = -grad_vy[i, j]
            end
        end
        for j in 1:grid_size
            ey[j, grid_size] = -grad_vy[j, grid_size-1]
        end
        """

        # python gradient calculation
        # TODO find julia solution
        np = pyimport("numpy")
        grad_v2 = np.gradient(vs)
        grad_vx = grad_v2[1]
        grad_vy = grad_v2[2]
        ex = - grad_vx
        ey =  - grad_vy

        # calculate drift velocity
        vdx = Array{Float64}(undef, grid_size, grid_size)
        vdy = Array{Float64}(undef, grid_size, grid_size)
        for i in 1:grid_size
            for j in 1:grid_size
                B = Field.bd(gr[1][i], gr[1][j], psr)
                E = [ex[i, j], ey[i, j], 0]
                #println(B)
                #println(E)
                v = cross(E, B)
                vdx[i, j] = v[1]
                vdy[i, j] = v[2]
            end
        end

        #println(typeof(grad_v2))
        psr.potential = vs
        psr.electric_field = [ex, ey]
        psr.drift_velocity = [vdx, vdy]
    end


    """
    Calculates electric potential, electric field and drift velocity for custom sparks defined by (sx, sy, sz)...
    """
    function calculate_potential_custom!(psr)
        gr = psr.grid
        grid_size = size(gr[1])[1]
        sp = psr.sparks
        spark_num = size(sp)[1]

        vs = Array{Float64}(undef, grid_size, grid_size)

        v_min =  1e50
        v_max =  -1e50

        # TODO start here
        for i in 1:grid_size
            for j in 1:grid_size
                vv = 0
                for k in 1:spark_num
                    sx = sp[k][1]
                    sy = sp[k][2]
                    sz = sp[k][3]
                    dist = norm([gr[1][i], gr[2][j], gr[3][i, j]] - [sx, sy, sz])
                    #vv += v(dist) # nice looking dots (Inf) in the plot, but no
                    if dist != 0
                        vv += v(dist)
                    end
                end
                if gr[3][i, j] != 0
                    vs[i, j] = vv
                else
                    vs[i, j] = 0
                end

            end
        end

        # set 0 potential (beyound the polar cap) to NaN
        for i in 1:grid_size
            for j in 1:grid_size
                if vs[i, j] == 0
                    vs[i, j] = NaN
                end
            end
        end

        # calculate electric field
        ex = Array{Float64}(undef, grid_size, grid_size)
        ey = Array{Float64}(undef, grid_size, grid_size)

        # python gradient calculation
        # TODO find julia solution check calculate_potential!() for inspiration
        np = pyimport("numpy")
        grad_v2 = np.gradient(vs)
        grad_vx = grad_v2[1]
        grad_vy = grad_v2[2]
        ex = - grad_vx
        ey =  - grad_vy

        # calculate drift velocity
        vdx = Array{Float64}(undef, grid_size, grid_size)
        vdy = Array{Float64}(undef, grid_size, grid_size)
        for i in 1:grid_size
            for j in 1:grid_size
                B = Field.bd(gr[1][i], gr[1][j], psr)
                E = [ex[i, j], ey[i, j], 0]
                #println(B)
                #println(E)
                v = cross(E, B)
                vdx[i, j] = v[1]
                vdy[i, j] = v[2]
            end
        end

        #println(typeof(grad_v2))
        psr.potential = vs
        psr.electric_field = [ex, ey]
        psr.drift_velocity = [vdx, vdy]
    end


    """

    Random sparks at the polar cap (grids will be generated later!)

    # Arguments

    - min_dist: minimum distant in meters

    """
    function random_sparks!(psr; min_dist=30, trials=1000)
        sp = []

        # maximum theta
        thm = asin((psr.r_pc - min_dist) / psr.r) # plus distance from pc boundry condition
        for i in 1:trials
            phi = rand() * 2pi
            theta = rand() * thm
            car = Functions.spherical2cartesian([psr.r, theta, phi])
            md = 2 * min_dist
            for s in sp
                dist = norm([s[1], s[2], s[3]] - [car[1], car[2], car[3]])
                if dist < md
                    md = dist
                end
            end
            if md > min_dist
                push!(sp, [car[1], car[2], car[3]])
            end
        end
        psr.sparks = sp
        println("Number of sparks added: ", size(sp)[1])
    end



    """

    Initiate sparks at the polar cap (grids will be generated later!)
    Same distance between sparks at the tracks

    # Arguments

    - rs: list of tracks radius
    - num: number of sparks at the inner track

    """
    function init_sparks1!(psr; rfs=[0.295, 0.5], num=3, center=true)
        sp = []

        # at the center
        if center == true
            car = Functions.spherical2cartesian([psr.r, 0, 0])
            push!(sp, [car[1], car[2], car[3]])
        end

        c1 = nothing
        for (i, rf) in enumerate(rfs)
            r = psr.r_pc * rf
            thm = asin(r / psr.r)
            if i == 1
                c1 = 2 * pi * r
                for phi in range(0, 2pi, length=num+1)[1:num] # really?
                    car = Functions.spherical2cartesian([psr.r, thm, phi])
                    push!(sp, [car[1], car[2], car[3]])
                    #println(phi)
                end
            else
                # calculate track cicumference to adjust track radius
                ci = 2 * pi * r
                num_new = convert(Int, ceil(ci / c1) * num)
                ci_new = ceil(ci / c1) * c1
                r_new = ci_new / (2 * pi)
                println("Radius of track no. $i adjusted to $(r_new/psr.r_pc)")
                thm = asin(r_new / psr.r)
                for phi in range(0, 2pi, length=num_new+1)[1:num_new] # really?
                    car = Functions.spherical2cartesian([psr.r, thm, phi])
                    push!(sp, [car[1], car[2], car[3]])
                    #println(phi)
                end
            end
            #println("$r $c1")
        end
        psr.sparks = sp
        println("Number of sparks added: ", size(sp)[1])
    end


    """

    Initiate sparks at the polar cap (grids will be generated later!)
    Same number of sparks at all tracks

    # Arguments

    - rs: list of tracks radius
    - num: number of sparks at the inner track and all the rest

    """
    function init_sparks2!(psr; rfs=[0.235, 0.71], num=3, offset=true, center=true)
        sp = []
        # at the center
        if center == true
            car = Functions.spherical2cartesian([psr.r, 0, 0])
            push!(sp, [car[1], car[2], car[3]])
        end

        dphi = 2pi / num
        for (i, rf) in enumerate(rfs)
            r = psr.r_pc * rf
            thm = asin(r / psr.r)
            for phi in range(0, 2pi, length=num+1)[1:num] # really?
                if offset == true
                    phi += (i % 2) * dphi / 2
                end
                car = Functions.spherical2cartesian([psr.r, thm, phi])
                push!(sp, [car[1], car[2], car[3]])
                #println(phi)
            end
            #println("$r $c1")
        end
        psr.sparks = sp
        println("Number of sparks added: ", size(sp)[1])
    end



    """
    Initiate sparks at the polar cap (grids will be generated later!)
    "Fibonacci" distribution

    # Arguments

    - num: number of sparks at the whole polar cap
    - rfmax: fraction of the polar cap radius covered by sparks

    """
    function init_sparks3!(psr; num=10, rfmax=0.7)
        sp = []
        # TODO start here
        r = psr.r
        theta_max = rfmax * Functions.theta_max(1, psr)
        rp = r * cos(theta_max) # calculation only at the polar cap
        phi = pi * (3. - sqrt(5.))  # golden angle in radians
        for i in 1:num
            zz = r - ((i-1) / (num-1)) * (r-rp) #  # 2 * r  # z goes from r to -r
            #println(zz)
            radius = sqrt(r^2 - zz * zz)  # radius at z
            theta = phi * i  # golden angle increment
            xx = cos(theta) * radius
            yy = sin(theta) * radius
            push!(sp, [xx, yy, zz])
        end
        psr.sparks = sp
        println("Number of sparks added: ", size(sp)[1])
    end


    """
    Creates grids around sparks to calculate gradient (for simulation)...

    # Arguments
    - prec: grid precision in meters
    - grid_size: odd number for proper spark calculation, even number for proper region calculation
    """
    function create_grids!(psr, prec=0.3, grid_size=5)

        if psr.sparks == nothing
            println("Run init_sparks! frirst..")
            return
        end

        sp = psr.sparks
        gr = []

        for s in sp
            gr_x = range(s[1] - prec*2, s[1]+prec*2, length=grid_size)
            gr_y = range(s[2] - prec*2, s[2]+prec*2, length=grid_size)
            gr_z = zeros((grid_size, grid_size))
            for (i, x) in enumerate(gr_x)
                for (j, y) in enumerate(gr_y)
                    z = sqrt(psr.r^2 - x^2 - y^2)
                    gr_z[i, j] = z
                end
            end
            push!(gr, [gr_x, gr_y, gr_z])
        end
        psr.grid = gr
    end


    """
    Calculates electric potential, electric field and drift velocity for grids around sparks
    """
    function calculate_potentials!(psr)
        grids = psr.grid
        grids_num = size(grids)[1]
        sp = psr.sparks
        spark_num = size(sp)[1]

        potentials = []
        electric_fields = []
        drift_velocities = []
        sparks_velocities = []

        psr.pot_minmax = [1e50, -1e50]

        for ii in 1:grids_num
            gr = grids[ii]
            grid_size = size(gr[1])[1]

            vs = Array{Float64}(undef, grid_size, grid_size)

            for i in 1:grid_size
                for j in 1:grid_size
                    vv = 0
                    for k in 1:spark_num
                        sx = sp[k][1]
                        sy = sp[k][2]
                        sz = sp[k][3]
                        dist = norm([gr[1][i], gr[2][j], gr[3][i, j]] - [sx, sy, sz])
                        #vv += v(dist) # nice looking dots (Inf) in the plot, but no
                        if dist != 0
                            vv += v(dist)
                        end
                    end
                    if vv < psr.pot_minmax[1]
                        psr.pot_minmax[1] = vv
                    end
                    if vv > psr.pot_minmax[2]
                        psr.pot_minmax[2] = vv
                    end
                    vs[i, j] = vv
                    #println(ii, " ", i, " ", j, " ", vs[i, j])
                end
            end
            push!(potentials, vs)

            # calculate electric field
            ex = Array{Float64}(undef, grid_size, grid_size)
            ey = Array{Float64}(undef, grid_size, grid_size)

            # python gradient calculation
            # TODO find julia solution
            np = pyimport("numpy")
            grad_v = np.gradient(vs)
            grad_vx = grad_v[1]
            grad_vy = grad_v[2]
            ex = - grad_vx
            ey = - grad_vy

            # calculate drift velocity
            vdx = Array{Float64}(undef, grid_size, grid_size)
            vdy = Array{Float64}(undef, grid_size, grid_size)
            for i in 1:grid_size
                for j in 1:grid_size
                    B = Field.bd(gr[1][i], gr[1][j], psr)
                    E = [ex[i, j], ey[i, j], 0]
                    #println(B)
                    #println(E)
                    v = cross(E, B)
                    vdx[i, j] = v[1]
                    vdy[i, j] = v[2]
                end
            end
            ind = convert(Int, ceil(grid_size / 2)) # works for odd grid sizes
            push!(sparks_velocities, [vdx[ind, ind], vdy[ind, ind]])
            push!(electric_fields, [ex, ey])
            push!(drift_velocities, [vdx, vdy])
        end
        #println(typeof(grad_v2))
        psr.potential = potentials
        psr.electric_field = electric_fields
        psr.drift_velocity = drift_velocities
        push!(psr.locations, deepcopy(psr.sparks))
        psr.sparks_velocity = sparks_velocities
        push!(psr.sparks_velocities, deepcopy(sparks_velocities))
    end


    function simulate!(psr; steps=1000, speedup=1)
        sp = psr.sparks

        for j in 1:steps
            sv = psr.sparks_velocity
            for (i,s) in enumerate(sp)
                """
                # prevent expanding hack! - not physical! NOPE!!!
                x1 = s[1]
                x2 = s[1] + sv[i][1] * speedup
                y1 = s[2]
                y2 = s[2] + sv[i][2] * speedup
                #println(norm([sv[i][1], sv[i][1]]))
                r = norm([x1, y1])
                phi = angle(x2 + y2*im)
                x2 = r * cos(phi)
                y2 = r * sin(phi)
                s[1] = x2
                s[2] = y2
                #"""
                s[1] = s[1] + sv[i][1] * speedup
                s[2] = s[2] + sv[i][2] * speedup
            end

            create_grids!(psr)
            calculate_potentials!(psr)
        end

    end

end
