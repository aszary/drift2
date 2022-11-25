module Plot
    using GLMakie
    using GeometryBasics
    GLMakie.activate!()
    #using PyPlot
    #using Plots
    #plotlyjs()
    #pyplot()

    function plot3d(psr)
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

    function plot3d2(psr)
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

    function plot3d3(psr)
        flat = collect(Iterators.flatten(psr.lines[1]))
        max = maximum(flat)

        # TODO fix the arrows!

        #fig = Figure(; resolution=(600, 480))
        #ax1 = Axis3(fig[1, 1]; aspect=(1,1,1), perspectiveness=0.5)
        #mesh!(ax1, Sphere(Point3f(0), psr.r), color=:blue, transparency=true)
        fig, ax1, p = mesh(Sphere(Point3f(0), psr.r), color=:blue, transparency=true) # better camera control

        arrows!(ax1, [Point3f(0, 0, 0)], [Point3f(psr.magnetic_axis[1], psr.magnetic_axis[2], psr.magnetic_axis[3])] , color=:red, linewidth =0.03* maximum(psr.magnetic_axis), normalize = false, arrowsize = Vec3f(0.05 * maximum(psr.magnetic_axis), 0.05* maximum(psr.magnetic_axis), 0.1* maximum(psr.magnetic_axis)), transparency=true)
        arrows!(ax1, [Point3f(0, 0, 0)], [Point3f(psr.rotation_axis[1], psr.rotation_axis[2], psr.rotation_axis[3])] , color=:green, linewidth =0.03* maximum(psr.rotation_axis), normalize = false, arrowsize = Vec3f(0.05 * maximum(psr.rotation_axis), 0.05* maximum(psr.rotation_axis), 0.1* maximum(psr.rotation_axis)), transparency=true)
        println(size(psr.lines))
        for l in psr.lines
            lines!(ax1, l[1], l[2], l[3])
        end
        #xlims!(ax1, -max, max)
        #ylims!(ax1, -max, max)
        #zlims!(ax1, -max, max)
        #arrows!(ax1, [Point3f(0, 0, 0)], [Point3f(100000, 0, 0)], linewidth = 1000.05, arrowsize = Vec3f(0.1, 0.1, 0.2), normalize = false, color = :red)
        #scatter!(ax1, x, y, z; markersize=15)
        display(fig)
    end



end  # module Plot
