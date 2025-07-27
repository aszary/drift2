module Drift2
    include("modules/functions.jl")
    include("modules/plot.jl")
    include("modules/lines.jl")
    include("modules/sparks.jl")
    include("modules/field.jl")
    using .Functions
    using .Lines
    using .Plot
    using .Sparks
    using .Field

    # define pulsar parameters
    mutable struct Pulsar
        p # pulsar period in [s]
        pdot # pulsar period derivative in [s/s]
        r # neutron star radius in [m] # TODO change?
        #beta_frac # fraction of r_pc covered by sparks path?
        r_pc # polar cap radius [in m]
        r_lc # light cylinder radius [in m]
        magnetic_axis # in cartesian coordinates
        rotation_axis # in cartesian coordinates
        omega
        omega_vec # Omega vector in cartesian components
        lines # magnetic lines (for ploting)
        pc # polar cap boundry x, y, z
        grid
        grid2
        sparks # single sparks locations
        locations # for the simulation
        sparks_velocity # single spark valocity
        sparks_velocities # for the simulation
        potential
        pot_minmax
        electric_field # at the polar cap
        drift_velocity # at the polar cap
        field_vacuum # magnetic and electric fields calculations
        field_forcefree # magnetic and electric fields calculations
        """
        align = 1 -> align rotator
        align = -1 -> anti-align rotator
        """
        function Pulsar(p, pdot, r; align=1, magnetic_axis=[0, 0, 2*r], rotation_axis=[0, 0, 1.5*r])
            r_pc = rdp(p, r)
            r_lc = rlc(p)
            omega = 2 * pi / p # SI units
            omega_vec = [0, 0, align * omega] # align rotator
            #sphere = generate_sphere(r) # GLMakie is the King :D
            return new(p, pdot, r, r_pc, r_lc, magnetic_axis, rotation_axis, omega, omega_vec, [], nothing, nothing, nothing, nothing, [], nothing, [], nothing, nothing, nothing, nothing, Field.Vacuum(), Field.ForceFree())
        end
    end


    function gradient3D_old()
        # initialize pulsar instance
        psr = Pulsar(1, 1e-15, 10e3) # period 1 s, radius 10 km

        #Lines.generate_dipole!(psr)
        Lines.calculate_polarcap!(psr)
        #Sparks.create_grid_sph!(psr) # does not work for plotting any more
        #Sparks.create_grid_fib!(psr) # nice but not very useful
        Sparks.create_grid_square!(psr)
        Sparks.random_sparks_old!(psr)
        Sparks.calculate_potential_old!(psr)
        Plot.plot3d_test3(psr)
    end


    function sparks_fullgrid()
        # initialize pulsar instance
        psr = Pulsar(1, 1e-15, 10e3) # period 1 s, radius 10 km

        #Lines.generate_dipole!(psr)
        Lines.calculate_polarcap!(psr)
        #Sparks.create_grid_sph!(psr) # does not work for plotting any more
        #Sparks.create_grid_fib!(psr) # nice but not very useful
        Sparks.create_grid!(psr)
        Sparks.random_sparks_grid!(psr)
        Sparks.calculate_potential!(psr)
        #Plot.potential3D(psr)
        Plot.potential2D(psr)
    end


    function sparks_smallgrids()
        # initialize pulsar instance
        psr = Pulsar(1, 1e-15, 10e3) # period 1 s, radius 10 km

        #Lines.generate_dipole!(psr)
        Lines.calculate_polarcap!(psr)
        #Sparks.random_sparks!(psr)
        #Sparks.init_sparks1!(psr)
        #Sparks.init_sparks2!(psr)
        Sparks.init_sparks3!(psr)

        Sparks.create_grid!(psr; size=100)
        Sparks.calculate_potential_custom!(psr)
        #Plot.potential2D(psr)
        Plot.potential2Dv2(psr)

        # simulation
        Sparks.create_grids!(psr)
        Sparks.calculate_potentials!(psr)
        #Plot.sparks(psr)
        Sparks.simulate!(psr)
        Plot.simulation2d(psr)
        Plot.simulation3d(psr)
    end


    function fields()

        #psr = Pulsar(1, 1e-15, 10e3) # period 1 s, pdot 1e-15, radius 10 km # aligned rotator
        #psr = Pulsar(1, 1e-15, 10e3; align=-1) # period 1 s,  pdot 1e-15, radius 10 km # anti-aligned rotator
        psr = Pulsar(1, 1e-15, 10e3; magnetic_axis = [1*10e3, 0, 0], rotation_axis = [0,0,1.5*10e3]) # period 1 s,  pdot 1e-15, radius 10 km # anti-aligned rotator

        #Lines.generate_dipole!(psr)
        #Field.calculate_eint!(psr) # useless?
        #Field.calculate_eint2!(psr) # useless?

        #Field.calculate_vac!(psr)
        #Lines.generate_vacuum!(psr)
        #Field.calculate_GJ!(psr)
        #Plot.field3d(psr, psr.field_vacuum)

        #Lines.generate_vacuum!(psr; phi=0) # phi=0 for 2d plot
        ##Field.calculate_GJ!(psr; twoD=true) # obsolete
        #Field.calculate_GJheat!(psr, psr.field_vacuum; size=1000)
        #Plot.vacuum2d(psr)
        
        # 3D here
        Field.calculate_ff!(psr)
        Lines.generate_forcefree!(psr)
        Field.calculate_GJ!(psr, field=psr.field_forcefree)
        Plot.field3d(psr, psr.field_forcefree)
        return
        
        #Lines.generate_forcefree!(psr; phi=0)
        #Field.calculate_GJheat!(psr, psr.field_forcefree; vacuum=false)
        #Plot.field2d(psr, psr.field_forcefree)

        # in the paper
        Field.calculate_vac!(psr)
        Lines.generate_vacuum!(psr; phi=0) # phi=0 for 2d plot
        Field.calculate_GJheat!(psr, psr.field_vacuum; size=1000)
        Field.calculate_ff!(psr)
        Lines.generate_forcefree!(psr; phi=0)
        Field.calculate_GJheat!(psr, psr.field_forcefree; vacuum=false, size=1000)
        #Plot.fields0(psr, psr.field_vacuum, psr.field_forcefree) # aligned, but change psr first!
        Plot.fields(psr, psr.field_vacuum, psr.field_forcefree) # anti-aligned

    end


    function polar_cap()

        psr = Pulsar(1, 1e-15, 10e3) # period 1 s, radius 10 km

        Lines.generate_dipole!(psr)
        Lines.calculate_polarcap!(psr)
        
        Sparks.init_sparks1!(psr, rfs=[90/psr.r_pc], num=8) # circle radius: 90 meters 
        
        Field.calculate_ff!(psr)
        Lines.generate_forcefree!(psr; phi=0)

        #Plot.polar_cap_obsolete(psr)
        Plot.polar_cap(psr) # in the paper
    end

    function test_andrzej()
        psr = Pulsar(1, 1e-15, 10e3) # period 1 s, radius 10 km
    
        Lines.calculate_polarcap!(psr)
    
        # calculate field
        Field.calculate_ff!(psr)
        Lines.generate_forcefree!(psr; phi=0)
    
        Sparks.create_grid!(psr, size=10)
        Sparks.random_sparks_grid!(psr)
        Sparks.calculate_potential!(psr)
        #Plot.potential3D(psr)
    
        #Plot.plot3d(psr)
        Plot.test_andrzej(psr)
    
    end
    


    function main()
        test_andrzej()
        #gradient3D_old()
        #sparks_fullgrid()
        #sparks_smallgrids()
        #fields() # plot in the paper
        #polar_cap() # plot in the paper
        return

        psr = Pulsar(1, 1e-15, 10e3) # period 1 s, radius 10 km

        Lines.generate_dipole!(psr)
        Lines.calculate_polarcap!(psr)
        
        Sparks.init_sparks1!(psr, rfs=[90/psr.r_pc], num=8) # circle radius: 90 meters 
        
        Field.calculate_ff!(psr)
        Lines.generate_forcefree!(psr; phi=0)

        #Plot.polar_cap_obsolete(psr)
        Plot.polar_cap(psr)


        println("Bye")
    end
end

d2 = Drift2.main()
