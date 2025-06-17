module Plot
    using Glob
    using GLMakie
    using CairoMakie
    using LinearAlgebra
    using ColorSchemes
    using GeometryBasics # no more Point3f?
    #using Meshes
    GLMakie.activate!()
    #CairoMakie.activate!()
    #using PyPlot
    #using Plots
    #plotlyjs()
    #pyplot()

    include("functions.jl")
    include("field.jl")
    using .Functions
    using .Field


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
        #println(fieldnames(typeof(fig))) # HERE useful!, dump is better
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
        heatmap!(ax, convert(Array{Float64},fi.xgj)/1e3 , convert(Array{Float64},fi.zgj)/1e3, convert(Array{Float64},fi.gj))

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

    function fields0(psr, field_vacuum, field_forcefree)

        fv = field_vacuum
        ff = field_forcefree

        # normalize fields to stellar radius [in kilometers]
        for i in 1:size(fv.magnetic, 1)
            fv.magnetic[i] =  fv.magnetic[i] / fv.beq * 0.5 * psr.r/1e3
            fv.electric[i] =  fv.electric[i] / fv.beq * 0.5 * psr.r/1e3
        end

        # normalize fields to stellar radius [in kilometers]
        for i in 1:size(ff.magnetic, 1)
            ff.magnetic[i] =  ff.magnetic[i] / ff.beq * 0.5 * psr.r/1e3
            ff.electric[i] =  ff.electric[i] / ff.beq * 0.5 * psr.r/1e3
        end


        CairoMakie.activate!()
        # Figure size
        size_inches = (17/2.54, 6.5/2.54) # 17cm x 6.5cm
        size_pt = 72 .* size_inches
        #println(size_pt)
        fig = Figure(resolution = size_pt, fontsize = 8)

        ax1 = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="x (km)", ylabel="z (km)") #, xminorticksvisible=true, yminorticksvisible=true)
       # charges
        heatmap!(ax1, convert(Array{Float64},fv.xgj)/1e3 , convert(Array{Float64},fv.zgj)/1e3, convert(Array{Float64},fv.gj))
        #heatmap!(ax1, convert(Array{Float64},fv.xgj)/1e3 , convert(Array{Float64},fv.zgj)/1e3, log10.(abs.(convert(Array{Float64},fv.gj)))) # log10
        # plot field lines
        for (i, l) in enumerate(fv.magnetic_lines)
            if i ==1
                lines!(ax1, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:royalblue1, linewidth=0.7, label="magnetic lines")
            else
                lines!(ax1, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:royalblue1, linewidth=0.7)
            end
        end
        for (i, l) in enumerate(fv.electric_lines)
            if i == 1
                lines!(ax1, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:indianred1, linewidth=0.7, label="electric lines")
            else
                lines!(ax1, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:indianred1, linewidth=0.7)
            end

        end
        lines!(ax1, Circle(Point2f(0, 0), psr.r/1e3), color=:palegreen1, linewidth=1)
        # neutron star
        #mesh!(ax1, Circle(Point2f(0, 0), psr.r/1e3), color=(:grey, 0.4))
        # magnetic axis
        arrows!(ax1, [0], [19], [0], [3], linewidth=0.77, arrowsize=3, color=:black, transparency=true)
        text!(ax1, 0.5, 25.7, text=L"\mathbf{\mu}", fontsize=8)
        # rotation axis
        arrows!(ax1, [0], [25], [0], [3], linewidth=0.77, arrowsize=3, color=:black, transparency=true)
        text!(ax1, 0.5, 18.3, text=L"\mathbf{\Omega}", fontsize=8)
        tightlimits!(ax1)
        hidexdecorations!(ax1, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        hideydecorations!(ax1, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        text!(ax1, -27, 25, text="a)", fontsize=9)
        xlims!(ax1, -30, 30)
        ylims!(ax1, -30, 30)
        

        ax2 = Axis(fig[1, 2]; aspect=DataAspect(), xlabel="x (km)", ylabel="z (km)") #, xminorticksvisible=true, yminorticksvisible=true)
        # charges
        hm = heatmap!(ax2, convert(Array{Float64},ff.xgj)/1e3 , convert(Array{Float64},ff.zgj)/1e3, convert(Array{Float64},ff.gj)) #, colorrange=[-6732.297190412704, 6732.297190412704])
        # TODO play with log scale?
        #hm = heatmap!(ax2, convert(Array{Float64},ff.xgj)/1e3 , convert(Array{Float64},ff.zgj)/1e3, log10.(abs.(convert(Array{Float64},ff.gj)))) #, colorrange=[-6732.297190412704, 6732.297190412704])
        cb = Colorbar(fig[1, 3], hm, label=L"$n_{GJ}$", vertical=true, height=Relative(0.7))
        #println(dump(cb,maxdepth=1))
        # plot field lines
        for (i, l) in enumerate(ff.magnetic_lines)
            if i ==1
                lines!(ax2, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:royalblue1, linewidth=0.7, label="magnetic lines")
            else
                lines!(ax2, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:royalblue1, linewidth=0.7)
            end
        end
        for (i, l) in enumerate(ff.electric_lines)
            if i == 1
                lines!(ax2, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:indianred1, linewidth=0.7, label="electric lines")
            else
                lines!(ax2, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:indianred1, linewidth=0.7)
            end

        end
        # neutron star
        lines!(ax2, Circle(Point2f(0, 0), psr.r/1e3), color=:palegreen1, linewidth=1)
        #mesh!(ax2, Circle(Point2f(0, 0), psr.r/1e3), color=(:grey, 0.4))
        # magnetic axis
        arrows!(ax2, [0], [19], [0], [3], linewidth=0.7, arrowsize=3, color=:black, transparency=true)
        text!(ax2, 0.5, 25.7, text=L"\mathbf{\mu}", fontsize=8)
        # rotation axis
        arrows!(ax2, [0], [25], [0], [3], linewidth=0.7, arrowsize=3, color=:black, transparency=true)
        text!(ax2, 0.5, 18.3, text=L"\mathbf{\Omega}", fontsize=8)
        tightlimits!(ax2)
        hidexdecorations!(ax2, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        hideydecorations!(ax2, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        text!(ax2, -27, 25, text="b)", fontsize=9)
        xlims!(ax2, -30, 30)
        ylims!(ax2, -30, 30)


        filename = "output/fields0.pdf"
        println(filename)
        save(filename, fig, pt_per_unit = 1)
        #display(fig)
    end


    function fields(psr, field_vacuum, field_forcefree)

        fv = field_vacuum
        ff = field_forcefree

        # normalize fields to stellar radius [in kilometers]
        for i in 1:size(fv.magnetic, 1)
            fv.magnetic[i] =  fv.magnetic[i] / fv.beq * 0.5 * psr.r/1e3
            fv.electric[i] =  fv.electric[i] / fv.beq * 0.5 * psr.r/1e3
        end

        # normalize fields to stellar radius [in kilometers]
        for i in 1:size(ff.magnetic, 1)
            ff.magnetic[i] =  ff.magnetic[i] / ff.beq * 0.5 * psr.r/1e3
            ff.electric[i] =  ff.electric[i] / ff.beq * 0.5 * psr.r/1e3
        end


        # to play with log scale in heatmaps
        xv, zv, gjv = [], [], [], []
        x2v, z2v, gj2v = [], [], [], []
        for i in 1:size(fv.gj, 1)
            if fv.gj[i] > 0
                push!(xv, fv.xgj[i])
                push!(zv, fv.zgj[i])
                push!(gjv, log10(fv.gj[i]))
            else
                push!(x2v, fv.xgj[i])
                push!(z2v, fv.zgj[i])
                push!(gj2v, log10(abs(fv.gj[i])))
            end
        end

        # to play with log scale in heatmaps
        x, z, gj = [], [], [], []
        x2, z2, gj2 = [], [], [], []
        for i in 1:size(ff.gj, 1)
            if ff.gj[i] > 0
                push!(x, ff.xgj[i])
                push!(z, ff.zgj[i])
                push!(gj, log10(ff.gj[i]))
            else
                push!(x2, ff.xgj[i])
                push!(z2, ff.zgj[i])
                push!(gj2, log10(abs(ff.gj[i])))
            end
        end
        cmap1 = ColorSchemes.Reds
        cmap2 = ColorSchemes.Blues

        CairoMakie.activate!()
        # Figure size
        size_inches = (17/2.54, 6.5/2.54) # 17cm x 6.5cm
        size_pt = 72 .* size_inches
        println(size_pt)
        fig = Figure(resolution = size_pt, fontsize = 8, figure_padding=(0, 5, 0, 5))

        ax1 = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="x (km)", ylabel="z (km)") #, xminorticksvisible=true, yminorticksvisible=true)
        # charges
        #heatmap!(ax1, convert(Array{Float64},fv.xgj)/1e3 , convert(Array{Float64},fv.zgj)/1e3, convert(Array{Float64},fv.gj))
        #heatmap!(ax1, convert(Array{Float64},fv.xgj)/1e3 , convert(Array{Float64},fv.zgj)/1e3, log10.(abs.(convert(Array{Float64},fv.gj)))) # log10
        hm = heatmap!(ax1, convert(Array{Float64}, xv)/1e3 , convert(Array{Float64},zv)/1e3, convert(Array{Float64},gjv); colormap=cmap2, colorrange=[9, 11])
        hm2 = heatmap!(ax1, convert(Array{Float64}, x2v)/1e3 , convert(Array{Float64},z2v)/1e3, convert(Array{Float64},gj2v); colormap=cmap1, colorrange=[9, 11])
        # plot field lines
        for (i, l) in enumerate(fv.magnetic_lines)
            if i ==1
                lines!(ax1, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:royalblue1, linewidth=0.7, label="magnetic lines")
            else
                lines!(ax1, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:royalblue1, linewidth=0.7)
            end
        end
        for (i, l) in enumerate(fv.electric_lines)
            if i == 1
                lines!(ax1, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:indianred1, linewidth=0.7, label="electric lines")
            else
                lines!(ax1, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:indianred1, linewidth=0.7)
            end
        end
        lines!(ax1, Circle(Point2f(0, 0), psr.r/1e3), color=:black, linewidth=1)
        # neutron star
        #mesh!(ax1, Circle(Point2f(0, 0), psr.r/1e3), color=(:grey, 0.4))
        # magnetic axis
        arrows!(ax1, [0], [19], [0], [3], linewidth=0.77, arrowsize=3, color=:black, transparency=true)
        text!(ax1, 0.5, 19, text=L"\mathbf{\mu}", fontsize=8)
        # rotation axis
        arrows!(ax1, [0], [-19], [0], [-3], linewidth=0.77, arrowsize=3, color=:black, transparency=true)
        text!(ax1, 1, -22.7, text=L"\mathbf{\Omega}", fontsize=8)
        tightlimits!(ax1)
        hidexdecorations!(ax1, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        hideydecorations!(ax1, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        text!(ax1, -27, 23, text="a)", fontsize=9)
        xlims!(ax1, -30, 30)
        ylims!(ax1, -30, 30)
        

        ax2 = Axis(fig[1, 4]; aspect=DataAspect(), xlabel="x (km)", ylabel="z (km)") #, xminorticksvisible=true, yminorticksvisible=true)
        # charges
        #hm = heatmap!(ax2, convert(Array{Float64},ff.xgj)/1e3 , convert(Array{Float64},ff.zgj)/1e3, convert(Array{Float64},ff.gj)) #, colorrange=[-6732.297190412704, 6732.297190412704])
        # TODO play with log scale?
        #last_color = cmap1.colors[end]
        #cmap2.colors[end] = last_color
        # :hsv
        hm = heatmap!(ax2, convert(Array{Float64}, x)/1e3 , convert(Array{Float64},z)/1e3, convert(Array{Float64},gj); colormap=cmap2, colorrange=[9, 11])
        hm2 = heatmap!(ax2, convert(Array{Float64}, x2)/1e3 , convert(Array{Float64},z2)/1e3, convert(Array{Float64},gj2); colormap=cmap1, colorrange=[9, 11])
        #hm = heatmap!(ax2, convert(Array{Float64},ff.xgj)/1e3 , convert(Array{Float64},ff.zgj)/1e3, log10.(abs.(convert(Array{Float64},ff.gj)))) #, colorrange=[-6732.297190412704, 6732.297190412704])
        #println(dump(cb,maxdepth=1))
        # plot field lines
        for (i, l) in enumerate(ff.magnetic_lines)
            if i ==1
                lines!(ax2, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:royalblue1, linewidth=0.7, label="magnetic lines")
            else
                lines!(ax2, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:royalblue1, linewidth=0.7)
            end
        end
        for (i, l) in enumerate(ff.electric_lines)
            if i == 1
                lines!(ax2, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:indianred1, linewidth=0.7, label="electric lines")
            else
                lines!(ax2, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:indianred1, linewidth=0.7)
            end

        end
        # neutron star
        lines!(ax2, Circle(Point2f(0, 0), psr.r/1e3), color=:black, linewidth=1)
        #mesh!(ax2, Circle(Point2f(0, 0), psr.r/1e3), color=(:grey, 0.4))
        # magnetic axis
        arrows!(ax2, [0], [19], [0], [3], linewidth=0.7, arrowsize=3, color=:black, transparency=true)
        text!(ax2, 0.5, 19, text=L"\mathbf{\mu}", fontsize=8)
        # rotation axis
        arrows!(ax2, [0], [-19], [0], [-3], linewidth=0.7, arrowsize=3, color=:black, transparency=true)
        text!(ax2, 1, -22.7, text=L"\mathbf{\Omega}", fontsize=8)
        tightlimits!(ax2)
        hidexdecorations!(ax2, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        hideydecorations!(ax2, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        text!(ax2, -27, 23, text="b)", fontsize=9)
        xlims!(ax2, -30, 30)
        ylims!(ax2, -30, 30)

        cb = Colorbar(fig[1,2], hm, label=L"$\log_{10}{(p^+)}$", vertical=true, width=Relative(0.1))
        cb = Colorbar(fig[1,3], hm2, label=L"$\log_{10}{(e^-)}$", vertical=true, width=Relative(0.1))
        colgap!(fig.layout, 1, 0)
        colgap!(fig.layout, 2, 0)
        colgap!(fig.layout, 3, 0)
        colsize!(fig.layout, 2, Relative(0.1)) # does not work?
        colsize!(fig.layout, 3, Relative(0.1)) # does not work?

        filename = "output/fields.pdf"
        println(filename)
        save(filename, fig, pt_per_unit = 1)
        #save(replace(filename, ".pdf"=>".png"), fig, pt_per_unit=1) # low-res
        # pdftoppm fields.pdf fields -png -r 300
        #display(fig)

    end


    function polar_cap_obsolete(psr)
        spark_radius = 20 # in meters

        ff = psr.field_forcefree
        # normalize fields to stellar radius [in kilometers]
        for i in 1:size(ff.magnetic, 1)
            ff.magnetic[i] =  ff.magnetic[i] / ff.beq * 0.5 * psr.r/1e3
            ff.electric[i] =  ff.electric[i] / ff.beq * 0.5 * psr.r/1e3
        end

        # spark num to plot
        sn = 8

        sx = psr.sparks[sn][1]
        sy = psr.sparks[sn][2]

        pnum = 6 # points at which fields are calculated
        cov = 0.7 # points coverage
        vlength = 0.15 #  drift velocity lenght

        points = Array{Float64}(undef, pnum, 3)

        dx = 2 * cov * spark_radius / (pnum -1)
        for i in 1:pnum
            points[i, 1] = sx - cov * spark_radius + dx * (i - 1)
            points[i, 2] = sy
            points[i, 3] = sqrt(psr.r ^2 - points[i,1]^2 - points[i,2]^2)
        end

        # drift velocity for the LBC model
        vdl = Array{Float64}(undef, pnum, 3)
        # electric field in the LBC model
        el = Array{Float64}(undef, pnum, 3)
        # electric field in the MC model
        em = Array{Float64}(undef, pnum, 3)
        # drift velocity for the MC model
        vdm = Array{Float64}(undef, pnum, 3)
       
        for i in 1:pnum
            # electric fields
            el[i, 1] = - points[i, 1]
            el[i, 2] = - points[i, 2]
            el[i, 3] = 0 
            el[i,:] = el[i,:] / norm(el[i,:]) * cov * spark_radius # normalize the length
            em[i, 1] = sx - points[i, 1]
            em[i, 2] = sy - points[i, 2]
            em[i, 3] = 0 
            em[i,:] = em[i,:] / norm(em[i,:]) * vlength * 0.5 * spark_radius # normalize the length
            # magneitc field
            r = Functions.cartesian2spherical(points[i,:])
            b_sph = Field.dipole(1, r[2]) # r_theta used    
            b_car = Functions.vec_spherical2cartesian(r, b_sph)
            # drift velocities
            v = cross(el[i, :], b_car)
            vdl[i, :] = v / norm(v) * vlength * spark_radius # normalize
            vm = cross(em[i, :], b_car)
            vdm[i, :] = vm / norm(vm) * 2 * vlength * spark_radius # normalize
        end

        # Generate curved arrow (for the top panel)
        ome = Array{Float64}(undef,2, 100)
        ph = range(0.8*pi, 2.1*pi, length=100)
        for i in 1:100
            ome[1, i] = 35 * cos(ph[i]) 
            ome[2, i] = 35 * sin(ph[i])
        end

        # generate civulation line
        cilx = collect(range(sx - 1.2* spark_radius, sx+1.2*spark_radius, length=100) )
        cily = fill(sy, 100)


        CairoMakie.activate!()

        # Figure size
        size_inches = (17/2.54, 11/2.54) # 17cm x 11cm
        size_pt = 72 .* size_inches
        #println(size_pt)
        fig = Figure(resolution=size_pt, fontsize=8, figure_padding=(1, 2, 0, 0)) # left, right, bottom, top

        # TODO:
        # - charges between sparks..
        # - reorder velocities top-down...
        top = Axis(fig[1, 2]; aspect=DataAspect(), xlabel="x (m)", ylabel="y (m)", xminorticksvisible=true, yminorticksvisible=true, xaxisposition=:top)
        hidexdecorations!(top, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        hideydecorations!(top, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        # plot polar cap boundry
        lines!(top, psr.pc[1], psr.pc[2], color=:green)
        # plot sparks
        if psr.sparks != nothing
            for sp in psr.sparks
                lines!(top, Circle(Point2f(sp[1], sp[2]), spark_radius), color=:grey, linewidth=1)
                #scatter!(top, sp[1], sp[2], sp[3], marker=:xcross, color=:red)
            end
        end
        lines!(top, cilx, cily, linewidth=0.7, linestyle=:dash, color=:palegreen1)
        lines!(top, ome[1,:], ome[2, :], color=:black, linewidth=0.7)
        arrows!(top, [ome[1, 1]], [ome[2, 1]], [ome[1, 1]-ome[1, 2]], [ome[2, 1]-ome[2, 2]], arrowsize=5, color=:black)
        text!(top, -20, -65, text=L"\mathbf{\Omega}", fontsize=8)

        lef = Axis(fig[2, 1]; aspect=DataAspect(), xlabel="x (m)", ylabel="y (m)", xminorticksvisible=true, yminorticksvisible=true)
        hidexdecorations!(lef, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        hideydecorations!(lef, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        lines!(lef, Circle(Point2f(psr.sparks[sn][1], psr.sparks[sn][2]), spark_radius), color=:grey, linewidth=1)
        lines!(lef, cilx, cily, linewidth=0.7, linestyle=:dash, color=:palegreen1)
        for i in 1:pnum
            arrows!(lef, [points[i, 1]], [points[i, 2]], [el[i, 1]], [el[i, 2]], arrowsize=5, color=:indianred)
            arrows!(lef, [points[i, 1]], [points[i, 2]], [vdl[i, 1]], [vdl[i, 2]], arrowsize=5, color=:black)
        end
        text!(lef, sx-2.1, sy+5, text=L"\mathbf{E_{\perp}^{\prime}}", fontsize=8, color=:indianred)
        text!(lef, sx-2.1, sy-7, text=L"\mathbf{v_{d}^{\prime}}", fontsize=8, color=:black)


        mid = Axis(fig[2, 2]; aspect=DataAspect(), xlabel="x (m)", ylabel="y (m)", xminorticksvisible=true, yminorticksvisible=true)
        hidexdecorations!(mid, label=true, ticklabels=true, ticks=true, grid=true, minorgrid=true, minorticks=true)
        hideydecorations!(mid, label=true, ticklabels=true, ticks=true, grid=true, minorgrid=true, minorticks=true)
        hidespines!(mid)
        # neutron star
        lines!(mid, Circle(Point2f(0, 0), psr.r/1e3), color=:palegreen1, linewidth=0.5)
        # plot field lines
        for (i, l) in enumerate(ff.magnetic_lines)
            lines!(mid, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:black, linewidth=0.1)
        end
        xlims!(mid, -30, 30)
        ylims!(mid, -30, 70)
        for l in psr.lines
            #plot!(mid, l[1], l[3])
        end

        rig = Axis(fig[2, 3]; aspect=DataAspect(), xlabel="x (m)", ylabel="y (m)", xminorticksvisible=true, yminorticksvisible=true, yaxisposition=:right)
        hidexdecorations!(rig, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        hideydecorations!(rig, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)

        lines!(rig, Circle(Point2f(psr.sparks[sn][1], psr.sparks[sn][2]), spark_radius), color=:grey, linewidth=1)
        lines!(rig, cilx, cily, linewidth=0.7, linestyle=:dash, color=:palegreen1)
        for i in 1:pnum
            arrows!(rig, [points[i, 1]], [points[i, 2]], [em[i, 1]], [em[i, 2]], arrowsize=5, color=:indianred)
            arrows!(rig, [points[i, 1]], [points[i, 2]], [vdm[i, 1]], [vdm[i, 2]], arrowsize=5, color=:black)
        end
        text!(rig, sx-3.1, sy+1, text=L"\mathbf{E_{\perp}^{\prime}}", fontsize=8, color=:indianred)
        text!(rig, sx-2.1, sy-7, text=L"\mathbf{v_{d}^{\prime}}", fontsize=8, color=:black)

        filename = "output/polar_cap_obsolete.pdf"
        filename2 = replace(filename, ".pdf"=>".svg")
        println(filename)
        save(filename, fig, pt_per_unit = 1)
        #println(filename2)
        #save(filename2, fig, pt_per_unit = 1) # slow
        #display(fig)

    end

    """
    anti-aligned case
    """
    function polar_cap(psr)
        spark_radius = 20 # in meters

        ff = psr.field_forcefree
        # normalize fields to stellar radius [in kilometers]
        for i in 1:size(ff.magnetic, 1)
            ff.magnetic[i] =  ff.magnetic[i] / ff.beq * 0.5 * psr.r/1e3
            ff.electric[i] =  ff.electric[i] / ff.beq * 0.5 * psr.r/1e3
        end

        # spark num to plot
        sn = 8

        sx = psr.sparks[sn][1]
        sy = psr.sparks[sn][2]

        pnum = 6 # points at which fields are calculated
        cov = 0.7 # points coverage
        vlength = 0.15 #  drift velocity lenght

        points = Array{Float64}(undef, pnum, 3)

        dy = 2 * cov * spark_radius / (pnum -1)
        for i in 1:pnum
            points[i, 1] = sx
            points[i, 2] = sy - cov * spark_radius + dy * (i - 1)
            points[i, 3] = sqrt(psr.r ^2 - points[i,1]^2 - points[i,2]^2)
        end

        # drift velocity for the LBC model
        vdl = Array{Float64}(undef, pnum, 3)
        # electric field in the LBC model
        el = Array{Float64}(undef, pnum, 3)
        # electric field in the MC model
        em = Array{Float64}(undef, pnum, 3)
        # drift velocity for the MC model
        vdm = Array{Float64}(undef, pnum, 3)
       
        for i in 1:pnum
            # electric fields
            el[i, 1] = - points[i, 1]
            el[i, 2] = - points[i, 2]
            el[i, 3] = 0 
            el[i,:] = el[i,:] / norm(el[i,:]) * vlength * spark_radius # normalize the length
            em[i, 1] = sx - points[i, 1]
            em[i, 2] = sy - points[i, 2]
            em[i, 3] = 0 
            em[i,:] = em[i,:] / norm(em[i,:]) * vlength * 0.5 * spark_radius # normalize the length
            # magneitc field
            r = Functions.cartesian2spherical(points[i,:])
            b_sph = Field.dipole(1, r[2]) # r_theta used    
            b_car = Functions.vec_spherical2cartesian(r, b_sph)
            # drift velocities
            v = cross(el[i, :], b_car)
            vdl[i, :] = v / norm(v) * 2* vlength * spark_radius # normalize
            vm = cross(em[i, :], b_car)
            vdm[i, :] = vm / norm(vm) * 2 * vlength * spark_radius # normalize
        end

        # Generate curved arrow (for the top panel)
        ome = Array{Float64}(undef,2, 100)
        ph = range(0.8*pi, 2.1*pi, length=100)
        for i in 1:100
            ome[1, i] = 35 * cos(ph[i]) 
            ome[2, i] = 35 * sin(ph[i])
        end

        # generate circulation line
        cilx = fill(sx, 100)
        cily = collect(range(sy - 1.1* spark_radius, sy+1.1*spark_radius, length=100) )

        # axis limits
        xyl = (-170, 170)
        xl = (sx-30, sx+30)
        yl = [sy-25, sy+25]


        CairoMakie.activate!()

        # Figure size
        size_inches = (17/2.54, 11/2.54) # 17cm x 11cm
        size_pt = 72 .* size_inches
        #println(size_pt)
        fig = Figure(resolution=size_pt, fontsize=8, figure_padding=(1, 2, 0, 0)) # left, right, bottom, top

        # TODO:
        # - charges between sparks..
        # - reorder velocities top-down...
        top = Axis(fig[1, 2]; aspect=DataAspect(), xlabel="x (m)", ylabel="y (m)", xminorticksvisible=true, yminorticksvisible=true, xaxisposition=:top)
        hidexdecorations!(top, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        hideydecorations!(top, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        # plot polar cap boundry
        mesh!(top, Rect(xyl[1], xyl[1], xyl[2]-xyl[1], xyl[2]-xyl[1]), color=(:grey93))
        lines!(top, psr.pc[1], psr.pc[2], color=:green)
        # plot sparks
        if psr.sparks != nothing
            for sp in psr.sparks
                mesh!(top, Circle(Point2f(sp[1], sp[2]), spark_radius), color=:white)
                lines!(top, Circle(Point2f(sp[1], sp[2]), spark_radius), color=:grey, linewidth=1)

                #scatter!(top, sp[1], sp[2], sp[3], marker=:xcross, color=:red)
            end
        end
        lines!(top, cilx, cily, linewidth=0.7, linestyle=:dash, color=:palegreen1)
        lines!(top, ome[1,:], ome[2, :], color=:black, linewidth=0.7)
        arrows!(top, [ome[1, 1]], [ome[2, 1]], [ome[1, 1]-ome[1, 2]], [ome[2, 1]-ome[2, 2]], arrowsize=5, color=:black)
        text!(top, 10, 20, text=L"\mathbf{\Omega}", fontsize=8)
        xlims!(xyl[1], xyl[2])
        ylims!(xyl[1], xyl[2])


        lef = Axis(fig[2, 1]; aspect=DataAspect(), xlabel="x (m)", ylabel="y (m)", xminorticksvisible=true, yminorticksvisible=true)
        hidexdecorations!(lef, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        hideydecorations!(lef, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        mesh!(lef, Rect(xl[1], yl[1], xl[2]-xl[1], yl[2]-yl[1]), color=(:grey93))
        mesh!(lef, Circle(Point2f(sx, sy), spark_radius), color=(:white))
        lines!(lef, Circle(Point2f(psr.sparks[sn][1], psr.sparks[sn][2]), spark_radius), color=:grey, linewidth=1)
        lines!(lef, cilx, cily, linewidth=0.7, linestyle=:dash, color=:palegreen1)
        for i in 1:pnum
            arrows!(lef, [points[i, 1]], [points[i, 2]], [el[i, 1]], [el[i, 2]], arrowsize=5, color=:indianred)
            arrows!(lef, [points[i, 1]], [points[i, 2]], [vdl[i, 1]], [vdl[i, 2]], arrowsize=5, color=:black)
        end
        text!(lef, sx-7, sy+5, text=L"\mathbf{E_{\perp}^{\prime}}", fontsize=8, color=:indianred)
        text!(lef, sx+7, sy-7, text=L"\mathbf{v_{d}^{\prime}}", fontsize=8, color=:black)
        # magnetic field
        mesh!(lef, Circle(Point2f(sx-spark_radius/3, sy-spark_radius/3), spark_radius/65), color=:orange)
        lines!(lef, Circle(Point2f(sx-spark_radius/3, sy-spark_radius/3), spark_radius/20), color=:orange, linewidth=0.5)
        text!(lef, sx-spark_radius/3-1, sy-spark_radius/3+1, text=L"\mathbf{B}", fontsize=8, color=:orange)
        xlims!(xl[1], xl[2])
        ylims!(yl[1], yl[2])


        mid = Axis(fig[2, 2]; aspect=DataAspect(), xlabel="x (m)", ylabel="y (m)", xminorticksvisible=true, yminorticksvisible=true)
        hidexdecorations!(mid, label=true, ticklabels=true, ticks=true, grid=true, minorgrid=true, minorticks=true)
        hideydecorations!(mid, label=true, ticklabels=true, ticks=true, grid=true, minorgrid=true, minorticks=true)
        hidespines!(mid)
        # neutron star
        #lines!(mid, Circle(Point2f(0, 0), psr.r/1e3), color=:palegreen1, linewidth=0.5)
        lines!(mid, Circle(Point2f(0, 0), psr.r/1e3), color=:black, linewidth=0.7)
        # plot field lines
        for (i, l) in enumerate(ff.magnetic_lines)
            lines!(mid, convert(Array{Float64}, l[1])/1e3, convert(Array{Float64}, l[3])/1e3, color=:black, linewidth=0.1)
        end
        xlims!(mid, -30, 30)
        ylims!(mid, -30, 70)
        for l in psr.lines
            #plot!(mid, l[1], l[3])
        end
        arrows!(mid, [0], [19], [0], [3], linewidth=0.77, arrowsize=3, color=:black, transparency=true)
        text!(mid, 1.7, 18.3, text=L"\mathbf{\mu}", fontsize=8)
        # rotation axis
        arrows!(mid, [0], [-19], [0], [-3], linewidth=0.77, arrowsize=3, color=:black, transparency=true)
        text!(mid, 1.7, -23, text=L"\mathbf{\Omega}", fontsize=6)


        rig = Axis(fig[2, 3]; aspect=DataAspect(), xlabel="x (m)", ylabel="y (m)", xminorticksvisible=true, yminorticksvisible=true, yaxisposition=:right)
        hidexdecorations!(rig, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        hideydecorations!(rig, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=true, minorticks=false)
        mesh!(rig, Rect(xl[1], yl[1], xl[2]-xl[1], yl[2]-yl[1]), color=(:grey93))
        mesh!(rig, Circle(Point2f(sx, sy), spark_radius), color=(:white))
        lines!(rig, Circle(Point2f(psr.sparks[sn][1], psr.sparks[sn][2]), spark_radius), color=:grey, linewidth=1)
        lines!(rig, cilx, cily, linewidth=0.7, linestyle=:dash, color=:palegreen1)
        for i in 1:pnum
            arrows!(rig, [points[i, 1]], [points[i, 2]], [em[i, 1]], [em[i, 2]], arrowsize=5, color=:indianred)
            arrows!(rig, [points[i, 1]], [points[i, 2]], [vdm[i, 1]], [vdm[i, 2]], arrowsize=5, color=:black)
        end
        text!(rig, sx+3, sy+5, text=L"\mathbf{E_{\perp}^{\prime}}", fontsize=8, color=:indianred)
        text!(rig, sx+7, sy-7, text=L"\mathbf{v_{d}^{\prime}}", fontsize=8, color=:black)
        # magnetic field
        mesh!(rig, Circle(Point2f(sx-spark_radius/3, sy-spark_radius/3), spark_radius/65), color=:orange)
        lines!(rig, Circle(Point2f(sx-spark_radius/3, sy-spark_radius/3), spark_radius/20), color=:orange, linewidth=0.5)
        text!(rig, sx-spark_radius/3-1, sy-spark_radius/3+1, text=L"\mathbf{B}", fontsize=8, color=:orange)
        xlims!(xl[1], xl[2])
        ylims!(yl[1], yl[2])


        filename = "output/polar_cap.pdf"
        filename2 = replace(filename, ".pdf"=>".svg")
        println(filename)
        save(filename, fig, pt_per_unit = 1)
        #println(filename2)
        #save(filename2, fig, pt_per_unit = 1) # slow # not vector graphic anymore!
        #display(fig)

    end
      function test(psr)
         # calculate electric potential
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
                # przygotuj siatkę dla heatmapy
        heatmap_x = reshape(x, grid_size, grid_size)
        heatmap_y = reshape(y, grid_size, grid_size)
        heatmap_v = reshape(v, grid_size, grid_size)
        #z_plane = fill(psr.r, size(heatmap_x))  # wysokość - na powierzchni pulsara

        z_plane = similar(heatmap_x)  # ta sama wielkość i typ

        for i in 1:size(heatmap_x, 1), j in 1:size(heatmap_x, 2)
            x_val = heatmap_x[i, j]
            y_val = heatmap_y[i, j]
            val = psr.r^2 - x_val^2 - y_val^2
            if val >= 0
                z_plane[i, j] = sqrt(val)   # górna półkula
            else
                z_plane[i, j] = NaN         # poza sferą
            end
        end
        # rysuj heatmapę i pulsar
        #fig, ax1, p = mesh(Sphere(Point3f(0, 0, 0,), psr.r), color=:blue, transparency=true)
        fig = Figure()
        ax1 = Axis3(fig[1, 1])
        surface!(
            ax1,
            heatmap_x,
            heatmap_y,
            z_plane;
            color=heatmap_v,
            colormap=:plasma,
            shading=false
        )
        
        # kula pulsara
        #mesh!(ax1, Sphere(Point3f0(0, 0, 0), psr.r), color=:blue, transparency=true)

        # linia bieguna (opcjonalnie)
        
        lines!(ax1, psr.pc[1], psr.pc[2], psr.pc[3], color=:white)
        if psr.sparks != nothing
            for sp in psr.sparks
                 i, j = sp
                x_sp = gr[1][i]
                y_sp = gr[2][j]
                z_sp = gr[3][i, j]
                scatter!(ax1, [x_sp], [y_sp], [z_sp], marker=:xcross, color=:red)
            end
        end
        Colorbar(fig[1, 2], colormap=:plasma, limits=extrema(heatmap_v), label="Potencjał")

        display(fig)
            

        return

        #fig, ax1, p = mesh(Sphere(Point3f(0, 0, 0), psr.r), color=:blue, transparency=true)
        # plot polar cap
        #lines!(ax1, psr.pc[1], psr.pc[2], psr.pc[3])

        
        # plot sparks
        if psr.sparks != nothing
            for (i, j) in psr.sparks
                #scatter!(ax1, gr[1][i], gr[2][j], gr[3][i, j], marker=:xcross, color=:red)
            end
        end

        #ze = zeros(size(z))

        #heatmap!(ax1, x, y, v, interpolate=false) #, colorrange=[-155, -135])
        #hm = meshscatter!(ax1, x, y, ze; markersize=1.25, color=v, transparency=false)
        #arrows!(ax1, x, y, ex, ey)


        display(fig)

        #=
        #fv = field_vacuum
        #ff = field_forcefree
        ff = psr.field_forcefree
        fig, ax1, p = mesh(Sphere(Point3f(0, 0, 0,), psr.r), color=:blue, transparency=true)
        
        #mid = Axis(fig[2, 2]; aspect=DataAspect(), xlabel="x (m)", ylabel="y (m)", xminorticksvisible=true, yminorticksvisible=true)
        #hidexdecorations!(mid, label=true, ticklabels=true, ticks=true, grid=true, minorgrid=true, minorticks=true)
        #hideydecorations!(mid, label=true, ticklabels=true, ticks=true, grid=true, minorgrid=true, minorticks=true)
        #hidespines!(mid)

        lines!(ax1, psr.pc[1], psr.pc[2], psr.pc[3])

        for (i, l) in enumerate(ff.magnetic_lines)
            lines!(ax1, convert(Array{Float64}, l[1]), convert(Array{Float64}, l[2]), convert(Array{Float64}, l[3]), color=:black, linewidth=0.1)
        end
        cam = cameracontrols(ax1.scene)
        cam.eyeposition[] = Point3f0(1, 1, 100)  # zmiana widoku- działa zoom pod warunkiem, że sie zzoomuje od góry 
        #cam.lookat[] = Point3f0(0, 0, 0)             # look at the origin

        display(fig)
        
        =#
    end

    

end  # module Plot
