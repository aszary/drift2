module Plot
    using GLMakie
    using GeometryBasics
    GLMakie.activate!()
    #using PyPlot
    #using Plots
    #plotlyjs()
    #pyplot()

    # https://discourse.julialang.org/t/determining-xlims/72921/2
    #xlims(ax::Axis) = @lift ($(ax.finallimits).origin[1], $(ax.finallimits).origin[1] + $(ax.finallimits).widths[1])
    #ylims(ax::Axis) = @lift ($(ax.finallimits).origin[2], $(ax.finallimits).origin[2] + $(ax.finallimits).widths[2])


    function plot3d_test(psr)
        fig = figure()
        ax = fig.add_subplot(111, projection="3d")
        plot_surface(psr.sphere[1], psr.sphere[2], psr.sphere[3], alpha=0.3, color="blue") # why error here?
        quiver([0], [0], [0], [psr.magnetic_axis[1]], [psr.magnetic_axis[2]], [psr.magnetic_axis[3]], color="blue")
        quiver([0], [0], [0], [psr.rotation_axis[1]], [psr.rotation_axis[2]], [psr.rotation_axis[3]], color="red")
        #xmami = 0
        #ymami = 0
        #zmami = 0
        for l in psr.lines
            #scatter(l[1], l[2], l[3], marker=".")
            plot(l[1], l[2], l[3], alpha=1)
            #xmami = maximum(l[1]) - minimum(l[1])
            #ymami = maximum(l[2]) - minimum(l[2])
            #zmami = maximum(l[3]) - minimum(l[3])
        end
        #ax.set_aspect("equal", "box") # matplotlib 3.6 required
        #ax.set_box_aspect([xmami, ymami, zmami]) # ? https://stackoverflow.com/questions/8130823/set-matplotlib-3d-plot-aspect-ratio
        #println(psr.lines)
    end

    function plot3d_test2(psr)
        #gr()
        #pyplot()
        p = surface(psr.sphere[1], psr.sphere[2], psr.sphere[3])
        quiver!(p, [0], [0], [0], quiver=([psr.magnetic_axis[1]], [psr.magnetic_axis[2]], [psr.magnetic_axis[3]]))
        quiver!(p, [0], [0], [0], quiver=([psr.rotation_axis[1]], [psr.rotation_axis[2]], [psr.rotation_axis[3]]))

        for l in psr.lines
            plot!(p, l[1], l[2], l[3])
        end
        gui()
        #savefig("myplot.pdf")
    end

    """
    Some Makie tests here, Nice..
    """
    function plot3d_test3(psr)
        #flat = collect(Iterators.flatten(psr.lines[1]))
        #max = maximum(flat)

        #fig = Figure(; resolution=(600, 480))
        #ax1 = Axis3(fig[1, 1]; aspect=(1,1,1), perspectiveness=0.5)

        # plot ns
        #fig, ax1, p = mesh(Sphere(Point3f(0), psr.r), color=:blue, transparency=true) # better camera control Scene, but zlims does not work
        fig, ax1, p = mesh(Sphere(Point3f(0, 0, psr.r), 0.03), color=:blue, transparency=true) # better camera control (Scene), but zlims does not work
        #mesh!(ax1, Sphere(Point3f(0), psr.r), color=:blue, transparency=true)

        #arrows!(ax1, [Point3f(0, 0, 0)], [Point3f(psr.magnetic_axis[1], psr.magnetic_axis[2], psr.magnetic_axis[3])] , color=:red, linewidth =0.03* maximum(psr.magnetic_axis), normalize = false, arrowsize = Vec3f(0.05 * maximum(psr.magnetic_axis), 0.05* maximum(psr.magnetic_axis), 0.1* maximum(psr.magnetic_axis)), transparency=true)
        #arrows!(ax1, [Point3f(0, 0, 0)], [Point3f(psr.rotation_axis[1], psr.rotation_axis[2], psr.rotation_axis[3])] , color=:green, linewidth =0.03* maximum(psr.rotation_axis), normalize = false, arrowsize = Vec3f(0.05 * maximum(psr.rotation_axis), 0.05* maximum(psr.rotation_axis), 0.1* maximum(psr.rotation_axis)), transparency=true)

        # plot polar cap
        #lines!(ax1, psr.pc[1], psr.pc[2], psr.pc[3])

        # plot field lines
        for l in psr.lines
            lines!(ax1, l[1], l[2], l[3])
        end

        # plot grid
        if psr.grid != nothing
            #scatter!(ax1, psr.grid[1], psr.grid[2], psr.grid[3], marker=:diamond, color=:black)
        end


        # plot sparks
        if psr.sparks != nothing
            for sp in psr.sparks
                scatter!(ax1, sp[1], sp[2], sp[3], marker=:xcross, color=:red)
            end
        end

        # plot potential
        if psr.grid != nothing
            #hm = meshscatter!(ax1, psr.grid[1], psr.grid[2], psr.grid[3]; markersize=1.25, color=psr.potential, transparency=false)
            scatter!(ax1, psr.grid[1], psr.grid[2], psr.grid[3], marker=:diamond, color=psr.potential)
            #heatmap!(ax1, psr.grid[1], psr.grid[2], psr.potential, interpolate=false, colorrange=[-155, -135])
            #contourf!(ax1, psr.grid[1], psr.grid[2], psr.potential)
            #Colorbar(fig[1, 4], hm, label="values", height=Relative(0.5))
        end

        #xlims(ax::Axis3) = @lift ($(ax.finallimits).origin[1], $(ax.finallimits).origin[1] + $(ax.finallimits).widths[1])
        #println(xlims(ax1))
        #println(ax1.finallimits)

        #println(methodswith(typeof(ax1)))

        #xlims!(ax1, -max, max)
        #ylims!(ax1, -max, max)
        #zlims!(ax1, 9000, 11000)
        display(fig)
    end


    """
    Test sparks generation and other calculations for old grid/spark generation procedure
    """
    function plot3d(psr)
        #flat = collect(Iterators.flatten(psr.lines[1]))
        #max = maximum(flat)

        #fig = Figure(; resolution=(600, 480))
        #ax1 = Axis3(fig[1, 1]; aspect=(1,1,1), perspectiveness=0.5)

        # plot ns
        #fig, ax1, p = mesh(Sphere(Point3f(0), psr.r), color=:blue, transparency=true) # better camera control Scene, but zlims does not work
        fig, ax1, p = mesh(Sphere(Point3f(0, 0, psr.r), 0.03), color=:blue, transparency=true) # better camera control (Scene), but zlims does not work
        #mesh!(ax1, Sphere(Point3f(0), psr.r), color=:blue, transparency=true)

        #arrows!(ax1, [Point3f(0, 0, 0)], [Point3f(psr.magnetic_axis[1], psr.magnetic_axis[2], psr.magnetic_axis[3])] , color=:red, linewidth =0.03* maximum(psr.magnetic_axis), normalize = false, arrowsize = Vec3f(0.05 * maximum(psr.magnetic_axis), 0.05* maximum(psr.magnetic_axis), 0.1* maximum(psr.magnetic_axis)), transparency=true)
        #arrows!(ax1, [Point3f(0, 0, 0)], [Point3f(psr.rotation_axis[1], psr.rotation_axis[2], psr.rotation_axis[3])] , color=:green, linewidth =0.03* maximum(psr.rotation_axis), normalize = false, arrowsize = Vec3f(0.05 * maximum(psr.rotation_axis), 0.05* maximum(psr.rotation_axis), 0.1* maximum(psr.rotation_axis)), transparency=true)

        # plot polar cap
        #lines!(ax1, psr.pc[1], psr.pc[2], psr.pc[3])

        # plot field lines
        for l in psr.lines
            lines!(ax1, l[1], l[2], l[3])
        end

        # plot grid
        if psr.grid != nothing
            #scatter!(ax1, psr.grid[1], psr.grid[2], psr.grid[3], marker=:diamond, color=:black)
        end


        # plot sparks
        if psr.sparks != nothing
            for sp in psr.sparks
                scatter!(ax1, sp[1], sp[2], sp[3], marker=:xcross, color=:red)
            end
        end

        # plot potential
        if psr.grid != nothing
            #hm = meshscatter!(ax1, psr.grid[1], psr.grid[2], psr.grid[3]; markersize=1.25, color=psr.potential, transparency=false)
            scatter!(ax1, psr.grid[1], psr.grid[2], psr.grid[3], marker=:diamond, color=psr.potential)
            #heatmap!(ax1, psr.grid[1], psr.grid[2], psr.potential, interpolate=false, colorrange=[-155, -135])
            #contourf!(ax1, psr.grid[1], psr.grid[2], psr.potential)
            #Colorbar(fig[1, 4], hm, label="values", height=Relative(0.5))
        end

        #xlims(ax::Axis3) = @lift ($(ax.finallimits).origin[1], $(ax.finallimits).origin[1] + $(ax.finallimits).widths[1])
        #println(xlims(ax1))
        #println(ax1.finallimits)

        #println(methodswith(typeof(ax1)))

        #xlims!(ax1, -max, max)
        #ylims!(ax1, -max, max)
        #zlims!(ax1, 9000, 11000)
        display(fig)
    end


    function potential3D(psr)
        gr = psr.grid
        grid_size = size(gr[1])[1]

        # data for potential plotting
        x = Array{Float64}(undef, grid_size * grid_size)
        y = Array{Float64}(undef, grid_size * grid_size)
        z = Array{Float64}(undef, grid_size * grid_size)
        v = Array{Float64}(undef, grid_size * grid_size)
        ex = Array{Float64}(undef, grid_size * grid_size)
        ey = Array{Float64}(undef, grid_size * grid_size)

        ind = 0
        for i in 1:grid_size
            for j in 1:grid_size
                ind += 1
                x[ind] = gr[1][i]
                y[ind] = gr[2][j]
                z[ind] = gr[3][i,j]
                v[ind] = psr.potential[i, j]
                ex[ind] = psr.electric_field[1][i, j]
                ey[ind] = psr.electric_field[2][i, j]
            end
        end

        #fig = Figure(; resolution=(600, 480))
        #ax1 = Axis3(fig[1, 1]; aspect=(1,1,1), perspectiveness=0.5)
        #fig, ax1, p = mesh(Sphere(Point3f(0, 0, psr.r), 0.03), color=:blue, transparency=true)
        fig, ax1, p = mesh(Sphere(Point3f(0, 0, 0), 0.03), color=:blue, transparency=true)

        # plot polar cap
        #lines!(ax1, psr.pc[1], psr.pc[2], psr.pc[3])

        # plot sparks
        if psr.sparks != nothing
            for (i, j) in psr.sparks
                #scatter!(ax1, gr[1][i], gr[2][j], gr[3][i, j], marker=:xcross, color=:red)
            end
        end

        #ze = zeros(size(z))

        heatmap!(ax1, x, y, v, interpolate=false) #, colorrange=[-155, -135])
        #hm = meshscatter!(ax1, x, y, ze; markersize=1.25, color=v, transparency=false)
        arrows!(ax1, x, y, ex, ey)

        display(fig)
    end


    function potential2D(psr)
        gr = psr.grid
        grid_size = size(gr[1])[1]

        # data for potential plotting
        x = Array{Float64}(undef, grid_size * grid_size)
        y = Array{Float64}(undef, grid_size * grid_size)
        z = Array{Float64}(undef, grid_size * grid_size)
        v = Array{Float64}(undef, grid_size * grid_size)
        ex = Array{Float64}(undef, grid_size * grid_size)
        ey = Array{Float64}(undef, grid_size * grid_size)
        vx = Array{Float64}(undef, grid_size * grid_size)
        vy = Array{Float64}(undef, grid_size * grid_size)

        ind = 0
        for i in 1:grid_size
            for j in 1:grid_size
                ind += 1
                x[ind] = gr[1][i]
                y[ind] = gr[2][j]
                z[ind] = gr[3][i,j]
                v[ind] = psr.potential[i, j]
                ex[ind] = psr.electric_field[1][i, j]
                ey[ind] = psr.electric_field[2][i, j]
                vx[ind] = psr.drift_velocity[1][i, j]
                vy[ind] = psr.drift_velocity[2][i, j]
            end
        end

        fig, ax1, p = heatmap(x, y, v, interpolate=false) #, colorrange=[-155, -135])
        #hm = meshscatter!(ax1, x, y, ze; markersize=1.25, color=v, transparency=false)
        #arrows!(x, y, ex, ey, color=:white)
        arrows!(x, y, vx, vy, color=:white)

        display(fig)
    end

    """
    Converts grid_x[size], grid_y[size], grid_z[size, size] to x[size*size], y[size*size], z[size*size]
    gets potential, electric field and drift velocity
    """
    function get_gridxyz(gr, potential, electric_field, drift_velocity)
        # TODO add potentials poltting
        grid_size = size(gr[1])[1]

        x = Array{Float64}(undef, grid_size * grid_size)
        y = Array{Float64}(undef, grid_size * grid_size)
        z = Array{Float64}(undef, grid_size * grid_size)
        v = Array{Float64}(undef, grid_size * grid_size)
        ex = Array{Float64}(undef, grid_size * grid_size)
        ey = Array{Float64}(undef, grid_size * grid_size)
        vx = Array{Float64}(undef, grid_size * grid_size)
        vy = Array{Float64}(undef, grid_size * grid_size)

        ind = 0
        for i in 1:grid_size
            for j in 1:grid_size
                ind += 1
                x[ind] = gr[1][i]
                y[ind] = gr[2][j]
                z[ind] = gr[3][i,j]
                v[ind] = potential[i, j]
                ex[ind] = electric_field[1][i, j]
                ey[ind] = electric_field[2][i, j]
                vx[ind] = drift_velocity[1][i, j]
                vy[ind] = drift_velocity[2][i, j]
            end
        end
        return (x, y, z, v, ex, ey, vx, vy)
    end


    function sparks(psr)

        #fig = Figure(; resolution=(600, 480))
        #ax1 = Axis3(fig[1, 1]; aspect=(1,1,1), perspectiveness=0.5)
        fig, ax1, p = mesh(Sphere(Point3f(0, 0, psr.r), 0.03), color=:blue, transparency=true)
        #fig, ax1, p = mesh(Sphere(Point3f(0, 0, 0), 0.03), color=:blue, transparency=true)

        # plot polar cap
        lines!(ax1, psr.pc[1], psr.pc[2], psr.pc[3])

        # plot sparks
        if psr.sparks != nothing
            for s in psr.sparks
                scatter!(ax1, s[1], s[2], s[3], marker=:xcross, color=:red)
            end
        end

        # plot grids
        for (i, gr) in enumerate(psr.grid)
            x, y, z, v, ex, ey, vx, vy = get_gridxyz(gr, psr.potential[i], psr.electric_field[i], psr.drift_velocity[i])
            scatter!(ax1, x, y, z, marker=:diamond, color=v, colorrange=(psr.pot_minmax[1], psr.pot_minmax[2]))
            arrows!(ax1, x, y, z, vx, vy, zeros(length(vx)), color=:black) # drift velocity
            #arrows!(ax1, )
        end
        display(fig)
    end


    function simulation(psr)

        #fig = Figure(; resolution=(600, 480))
        #ax1 = Axis3(fig[1, 1]; aspect=(1,1,1), perspectiveness=0.5)
        fig, ax1, p = mesh(Sphere(Point3f(0, 0, psr.r), 0.03), color=:blue, transparency=true)
        #fig, ax1, p = mesh(Sphere(Point3f(0, 0, 0), 0.03), color=:blue, transparency=true)

        # plot polar cap
        lines!(ax1, psr.pc[1], psr.pc[2], psr.pc[3])

        for sp in psr.locations[1]
            scatter!(ax1, sp[1], sp[2], sp[3], marker=:diamond)
        end

        for (i, sv) in enumerate(psr.sparks_velocity)
            x = psr.locations[1][i][1]
            y = psr.locations[1][i][2]
            z = psr.locations[1][i][3]
            #println(x, " ",  y, " ", z)
            #println(sv[1], " ",  sv[2], " ", z)
            arrows!(ax1, [x], [y], [z], [sv[1]], [sv[2]], [0], color=:red) # drift velocity
        end



        display(fig)
    end




end  # module Plot
