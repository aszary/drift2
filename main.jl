module Drift2
    include("modules/functions.jl")
    include("modules/plot.jl")
    include("modules/lines.jl")
    using .Functions
    using .Lines
    using .Plot

    # define pulsar parameters
    mutable struct Pulsar
        p # pulsar period in [s]
        r # neutron star radius in [m]
        #beta_frac # fraction of r_pc covered by sparks path?
        r_pc # polar cap radius
        r_lc # light cylinder radius
        magnetic_axis # in cartesian coordinates
        rotation_axis # in cartesian coordinates
        lines # magnetic lines (for ploting)
        function Pulsar(p, r)
            r_pc = rdp(p, r)
            r_lc = rlc(p)
            #sphere = generate_sphere(r) # GLMakie the King :D
            return new(p, r, r_pc, r_lc, [0, 0, 2*r], [0, 0, 1.5*r], [])
        end
    end


    function main()

        # initialize pulsar instance
        psr = Pulsar(1, 10e3) # period 1 s, radius 10 km

        #psr.lines = [1,2,3]
        #println(psr.lines)

        generate_dipole!(psr)
        #Plot.plot3d(psr)
        #Plot.plot3d2(psr)
        Plot.plot3d3(psr)

        #println(psr.sphere)

        #println(spherical2cartesian([3, 0, 0]))
        println("Bye")
    end
end

d2 = Drift2.main()
