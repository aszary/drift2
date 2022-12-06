module Sparks
    using LinearAlgebra
    using PyCall
    include("functions.jl")
    using .Functions


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
    simple square grid, not applicable for gradient calculation
    """
    function create_grid_square!(psr, size=200)
        r = psr.r # stellar radius in meters
        x_min, x_max = extrema(psr.pc[1])
        y_min, y_max = extrema(psr.pc[2])
        z_min = psr.pc[3][1]

        theta_max = Functions.theta_max(1, psr)
        xs = range(x_min, x_max, length=size)
        ys = range(y_min, y_max, length=size)
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

        psr.potential = convert(Array{Float64,1}, vs)

    end



    """
    simple square grid, but different data size => applicable for gradient calculation

    # try this https://stackoverflow.com/questions/40338386/calculating-a-3d-gradient-with-unevenly-spaced-points next time?
    """
    function create_grid!(psr, size=100)
        r = psr.r # stellar radius in meters

        x_min, x_max = extrema(psr.pc[1])
        y_min, y_max = extrema(psr.pc[2])
        z_min = psr.pc[3][1]


        xs = range(x_min, x_max, length=size)
        ys = range(y_min, y_max, length=size)
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
    function random_sparks!(psr; min_dist=30, trials=1000)
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
                        #vv += v(dist) # nice looking dots (Inf) in the plot
                        if dist != 0
                            vv += v(dist)
                        end
                    end
                end
                if vv < v_min
                    v_min = vv
                end
                if vv > v_max
                    v_max = vv
                end
                if gr[3][i, j] != 0
                    vs[i, j] = vv
                else
                    vs[i, j] = 0
                end
            end
        end

        # set 0 potential (beyound the polar cap) to minimum
        #println(v_min)
        #println(v_max)
        for i in 1:grid_size
            for j in 1:grid_size
                if vs[i, j] == 0
                    vs[i, j] = v_max # min or max here depending on v() function!
                end
            end
        end

        # python gradient tests
        #np = pyimport("numpy")
        #grad_v2 = np.gradient(vs)

        # calculate electric field
        grad_vx = diff(vs, dims=1) # 199x200
        grad_vy = diff(vs, dims=2) # 200x199
        ex = Array{Float64}(undef, grid_size, grid_size)
        ey = Array{Float64}(undef, grid_size, grid_size)

        for i in 1:grid_size-1
            for j in 1:grid_size
                ex[i, j] = -grad_vx[i, j]
            end
        end
        # dirty hack
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

        #println(typeof(grad_v2))
        psr.potential = vs
        psr.electric_field = [ex, ey]
    end



end
