module Drift2
    include("modules/functions.jl")
    include("modules/plot.jl")
    include("modules/lines.jl")
    include("modules/sparks.jl")
    using .Functions
    using .Lines
    using .Plot
    using .Sparks

    # define pulsar parameters
    mutable struct Pulsar
        p # pulsar period in [s]
        r # neutron star radius in [m]
        #beta_frac # fraction of r_pc covered by sparks path?
        r_pc # polar cap radius [in m]
        r_lc # light cylinder radius [in m]
        magnetic_axis # in cartesian coordinates
        rotation_axis # in cartesian coordinates
        lines # magnetic lines (for ploting)
        pc # polar cap boundry x, y, z
        grid
        sparks
        potential
        electric_field
        function Pulsar(p, r)
            r_pc = rdp(p, r)
            r_lc = rlc(p)
            #sphere = generate_sphere(r) # GLMakie is the King :D
            return new(p, r, r_pc, r_lc, [0, 0, 2*r], [0, 0, 1.5*r], [], nothing, nothing, nothing, nothing, nothing)
        end
    end


    function test_plots()

        # initialize pulsar instance
        psr = Pulsar(1, 10e3) # period 1 s, radius 10 km

        #Lines.generate_dipole!(psr)
        Lines.calculate_polarcap!(psr)
        #Sparks.create_grid_sph!(psr) # does not work for plotting any more
        #Sparks.create_grid_fib!(psr) # nice but not very useful
        Sparks.create_grid_square!(psr)
        Sparks.random_sparks_old!(psr)
        Sparks.calculate_potential_old!(psr)
        Plot.plot3d_test3(psr)




    end


    function main()

        #test_plots()
        #return

        # initialize pulsar instance
        psr = Pulsar(1, 10e3) # period 1 s, radius 10 km

        #Lines.generate_dipole!(psr)
        Lines.calculate_polarcap!(psr)
        #Sparks.create_grid_sph!(psr) # does not work for plotting any more
        #Sparks.create_grid_fib!(psr) # nice but not very useful
        Sparks.create_grid!(psr)
        Sparks.random_sparks!(psr)
        Sparks.calculate_potential!(psr)
        Plot.potential(psr)


        println("Bye")
    end
end

d2 = Drift2.main()
