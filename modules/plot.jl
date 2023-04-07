module Plot
    using Glob
    using GLMakie
    using CairoMakie
    using LinearAlgebra
    #using GeometryBasics # no more Point3f?
    #using Meshes
    GLMakie.activate!()
    #CairoMakie.activate!()
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


    function potential2Dv2(psr)
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

        #fig = Figure(; resolution=(600, 600))
        #ax = Axis(fig[1, 1]; aspect=(1,1))
        #heatmap!(fig, x, y, v, interpolate=false) #, colorrange=[-155, -135])

        fig, ax1, p = heatmap(x, y, v, interpolate=false) #, colorrange=[-155, -135])
        resize!(fig, (700, 700)) # changes resolution
        #hm = meshscatter!(ax1, x, y, ze; markersize=1.25, color=v, transparency=false)

        #arrows!(x, y, ex, ey, color=:white)
        #arrows!(x, y, vx, vy, color=:white)

        # get last file
        filename = get_newfilename("output", "potential2D_", "png")
        println("Filename: $filename")

        save(filename, fig)
        #save("output/potential2D.svg",fig) # does not work
        #save("output/potential2D.pdf",fig) # does not work
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


    function simulation2d(psr)

        # plot polar cap
        fig, ax, pl = lines(psr.pc[1], psr.pc[2])
        #println(fieldnames(typeof(fig))) # HERE useful!
        ax.aspect = DataAspect()

        observables = []

        for (i, sp) in enumerate(psr.locations[1])
            #println(i)
            o1 = Observable(sp[1])
            o2 = Observable(sp[2])
            push!(observables, [o1, o2])
            scatter!(o1, o2, marker=:diamond) #, color=colors[i])
        end

        #=
        # plot grids
        for (i, gr) in enumerate(psr.grid)
            x, y, z, v, ex, ey, vx, vy = get_gridxyz(gr, psr.potential[i], psr.electric_field[i], psr.drift_velocity[i])
            scatter!(x, y, marker=:diamond, color=v, colorrange=(psr.pot_minmax[1], psr.pot_minmax[2]))
            arrows!(x, y, vx, vy, color=:black) # drift velocity
            #arrows!(ax1, )
        end
        =#

        #display(fig)

        # TODO start here

        # get last file
        filename = get_newfilename("output", "sparks_", "mp4")
        println("Filename: $filename")

        # save animation
        record(fig, filename, 1:length(psr.locations), framerate=600) do i
        # animation
        #display(fig)
        #for (i, loc) in enumerate(psr.locations)
        #   sleep(0.01)
            loc = psr.locations[i]
            for (j, obs) in  enumerate(observables)
                obs[1][] = loc[j][1]
                obs[2][] = loc[j][2]
            end
        end

    end

    """
    Returns new file name incremented by +1
    """
    function get_newfilename(dir, filestart, ext="mp4")
        # get last file
        files = glob("$filestart*.$ext", "$dir")
        if length(files) == 0
                return "$dir/$(filestart)1.$ext"
        end
        files = replace.(files, "$dir/$filestart"=>"", ".$ext"=>"")
        nums = parse.(Int, files)
        num = maximum(nums) + 1
        return "$dir/$filestart$num.$ext"
    end


    function simulation3d(psr)

        #fig = Figure(; resolution=(600, 480))
        #ax1 = Axis3(fig[1, 1]; aspect=(1,1,1), perspectiveness=0.5)
        fig, ax1, p = mesh(Sphere(Point3f(0, 0, psr.r), 0.03), color=:blue, transparency=true)
        #fig, ax1, p = mesh(Sphere(Point3f(0, 0, 0), 0.03), color=:blue, transparency=true)

        # plot polar cap
        lines!(ax1, psr.pc[1], psr.pc[2], psr.pc[3])

        observables = []
        for (i, sp) in enumerate(psr.locations[1])
            #println(i)
            o1 = Observable(sp[1])
            o2 = Observable(sp[2])
            o3 = Observable(sp[3])
            push!(observables, [o1, o2, o3])
            scatter!(ax1, o1, o2, o3, marker=:diamond, color=:black)
        end
        for (i, sv) in enumerate(psr.sparks_velocity)
            x = psr.locations[1][i][1]
            y = psr.locations[1][i][2]
            z = psr.locations[1][i][3]
            #println(x, " ",  y, " ", z)
            #println(sv[1], " ",  sv[2], " ", z)
            arrows!(ax1, [x], [y], [z], [sv[1]], [sv[2]], [0], color=:red) # drift velocity
        end


        # plot grids
        for (i, gr) in enumerate(psr.grid)
            x, y, z, v, ex, ey, vx, vy = get_gridxyz(gr, psr.potential[i], psr.electric_field[i], psr.drift_velocity[i])
            scatter!(ax1, x, y, z, marker=:diamond, color=v, colorrange=(psr.pot_minmax[1], psr.pot_minmax[2]))
            arrows!(ax1, x, y, z, vx, vy, zeros(length(vx)), color=:black) # drift velocity
            #arrows!(ax1, )
        end

        display(fig)

        # animation
        for (i,loc) in enumerate(psr.locations)
            for (j, obs) in  enumerate(observables)
                obs[1][] = loc[j][1]
                obs[2][] = loc[j][2]
                obs[3][] = loc[j][3]
            end

            println(i)
            sleep(0.01)
        end
    end


    function field3d(psr, field_class)
        # Field.Vacuum class
        fi = field_class

        #println(fi.electric[1])
        # normalize fields to stellar radius
        for i in 1:size(fi.magnetic, 1)
            fi.magnetic[i] =  fi.magnetic[i] / fi.beq * 0.5 * psr.r
            fi.electric[i] =  fi.electric[i] / fi.beq * 0.5 * psr.r
        end
        # interior and internal electrif fields disabled
        #=
        for i in 1:size(fi.eint, 1)
            fi.eint[i] = fi.eint[i] / fi.beq * 0.5 * psr.r
        end
        for i in 1:size(fi.eint2, 1)
            fi.eint2[i] = fi.eint2[i] / fi.beq * 0.5 * psr.r
        end
        =#

        x, y, z = [], [], []
        for i in 1:size(fi.gj3, 1)
            #println(fi.gj3[i])
            push!(x, fi.locs3[i][1])
            push!(y, fi.locs3[i][2])
            push!(z, fi.locs3[i][3])
        end

        #println(fi.electric[1])

        #println(norm(fi.magnetic[1]))

        #fig = Figure(; resolution=(600, 480))
        #ax1 = Axis3(fig[1, 1]; aspect=(1,1,1), perspectiveness=0.5)
        fig, ax1, p = mesh(Sphere(Point3f(0, 0, 0), psr.r), color=:white, transparency=true)
        
        # magnetic and electric fields
        #arrows!(ax1, Point3f.(fi.locations), Vec3.(fi.magnetic), color=:blue, arrowsize=Vec3f(0.1*psr.r, 0.1*psr.r, 0.2*psr.r), linewidth=0.05*psr.r)
        #arrows!(ax1, Point3f.(fi.locations), Vec3.(fi.electric), color=:red, arrowsize=Vec3f(0.1*psr.r, 0.1*psr.r, 0.2*psr.r), linewidth=0.05*psr.r)
        
        # interior electric field
        #arrows!(ax1, Point3f.(fi.locs), Vec3.(fi.eint), color=:indianred, arrowsize=Vec3f(0.1*psr.r, 0.1*psr.r, 0.2*psr.r), linewidth=0.05*psr.r)
        # internal electric field
        #arrows!(ax1, Point3f.(fi.locs2), Vec3.(fi.eint2), color=:orange, arrowsize=Vec3f(0.1*psr.r, 0.1*psr.r, 0.2*psr.r), linewidth=0.05*psr.r)

        # plot field lines
        for l in fi.magnetic_lines
            lines!(ax1, convert(Array{Float64}, l[1]), convert(Array{Float64}, l[2]), convert(Array{Float64}, l[3]))
        end

        for l in fi.electric_lines
            lines!(ax1, convert(Array{Float64}, l[1]), convert(Array{Float64}, l[2]), convert(Array{Float64}, l[3]), color=:red, linewidth=5)
        end

        meshscatter!(ax1, convert(Array{Float64},x) , convert(Array{Float64},y) , convert(Array{Float64},z), markersize = 510.5, color = convert(Array{Float64},fi.gj3), colormap=:seismic)
        #=
        for i in 1:size(fi.gj3, 1)
            if fi.gj3[i] > 0
                meshscatter!(ax1, x[i] ,y[i] , z[i], markersize = 510.5, color=:red)
            else
                meshscatter!(ax1, x[i] ,y[i] , z[i], markersize = 510.5, color=:blue)
            end
        end
        =#
        display(fig)
    end

    function vacuum2d(psr)
        # Field.Vacuum class
        fv = psr.field_vacuum

        # normalize fields to stellar radius [in kilometers]
        for i in 1:size(fv.magnetic, 1)
            fv.magnetic[i] =  fv.magnetic[i] / fv.beq * 0.5 * psr.r/1e3
            fv.electric[i] =  fv.electric[i] / fv.beq * 0.5 * psr.r/1e3
        end

        x, y = [], []
        for i in 1:size(fv.gj3, 1)
            #println(fv.gj3[i])
            push!(x, fv.locs3[i][1])
            push!(y, fv.locs3[i][3])
        end

        CairoMakie.activate!()
        # Figure size
        size_inches = (11/2.54, 11/2.54) # 11cm x 11cm
        size_pt = 72 .* size_inches
        #println(size_pt)
        fig = Figure(resolution = size_pt, fontsize = 8)

        #fig = Figure(resolution=(900, 600))
        ax = Axis(fig[1, 1]; aspect=DataAspect())

        # charges
        heatmap!(ax, convert(Array{Float64},fv.xgj)/1e3 , convert(Array{Float64},fv.zgj)/1e3, convert(Array{Float64},fv.gj))

        # plot field lines
        for (i, l) in enumerate(fv.magnetic_lines)
            if i ==1
                lines!(ax, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:blue, linewidth=0.3, label="magnetic lines")
            else
                lines!(ax, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:blue, linewidth=0.3)
            end
        end

        for (i, l) in enumerate(fv.electric_lines)
            if i == 1
                lines!(ax, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:red, linewidth=0.3, label="electric lines")
            else
                lines!(ax, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:red, linewidth=0.3)
            end

        end
        #lines!(ax, Circle(Point2f(0, 0), psr.r), color=:grey, linewidth=3)
        # neutron star
        #mesh!(ax, Circle(Point2f(0, 0), psr.r/1e3), color=(:grey, 0.4))

        # charges # obsolete
        #=
        for i in 1:size(fv.gj3, 1)
            if fv.gj3[i] > 0
                scatter!(ax, fv.locs3[i][1]/1e3, fv.locs3[i][3]/1e3, color=:indianred, markersize=3)
            else
                scatter!(ax, fv.locs3[i][1]/1e3, fv.locs3[i][3]/1e3, color=:royalblue, markersize=3)
            end
        end
        =#

        #Legend(fig, ax)
        #axislegend(ax) # TODO add later
        # in kilometers now
        #xlims!(ax, -90, 90)
        #ylims!(ax, -60, 60)
        xlims!(ax, -30, 30)
        ylims!(ax, -20, 20)
        #tightlimits!(ax)

        filename = "output/vacuum2d.pdf"
        println(filename)
        save(filename, fig, pt_per_unit = 1)
        #display(fig)
    end


    function field2d(psr, field_class)
        fi = field_class

        # normalize fields to stellar radius [in kilometers]
        for i in 1:size(fi.magnetic, 1)
            fi.magnetic[i] =  fi.magnetic[i] / fi.beq * 0.5 * psr.r/1e3
            fi.electric[i] =  fi.electric[i] / fi.beq * 0.5 * psr.r/1e3
        end

        CairoMakie.activate!()
        # Figure size
        size_inches = (11/2.54, 11/2.54) # 11cm x 11cm
        size_pt = 72 .* size_inches
        #println(size_pt)
        fig = Figure(resolution = size_pt, fontsize = 8)

        #fig = Figure(resolution=(900, 600))
        ax = Axis(fig[1, 1]; aspect=DataAspect())

        # charges
        heatmap!(ax, convert(Array{Float64},fi.xgj)/1e3 , convert(Array{Float64},fi.zgj)/1e3, convert(Array{Float64},fi.gj), transparency=true)

        # plot field lines
        for (i, l) in enumerate(fi.magnetic_lines)
            if i ==1
                lines!(ax, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:blue, linewidth=0.3, label="magnetic lines")
            else
                lines!(ax, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:blue, linewidth=0.3)
            end
        end

        for (i, l) in enumerate(fi.electric_lines)
            if i == 1
                lines!(ax, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:red, linewidth=0.3, label="electric lines")
            else
                lines!(ax, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:red, linewidth=0.3)
            end

        end
        #lines!(ax, Circle(Point2f(0, 0), psr.r), color=:grey, linewidth=3)
        # neutron star
        #mesh!(ax, Circle(Point2f(0, 0), psr.r/1e3), color=(:grey, 0.4))

        # charges
        #=
        for i in 1:size(fi.gj3, 1)
            if fi.gj3[i] > 0
                scatter!(ax, fi.locs3[i][1]/1e3, fi.locs3[i][3]/1e3, color=:indianred, markersize=3)
            else
                scatter!(ax, fi.locs3[i][1]/1e3, fi.locs3[i][3]/1e3, color=:royalblue, markersize=3)
            end
        end
        =#

        #Legend(fig, ax)
        #axislegend(ax) # TODO add later
        # in kilometers now
        #xlims!(ax, -90, 90)
        #ylims!(ax, -60, 60)
        xlims!(ax, -30, 30)
        ylims!(ax, -20, 20)
        #tightlimits!(ax)

        filename = "output/field2d.pdf"
        println(filename)
        save(filename, fig, pt_per_unit = 1)
        #display(fig)
    end

end  # module Plot
