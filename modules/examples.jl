using GLMakie, LinearAlgebra

n = 10

θ = [0;(0.5:n-0.5)/n;1]
φ = [(0:2n-2)*2/(2n-1);2]
x = [cospi(φ)*sinpi(θ) for θ in θ, φ in φ]
y = [sinpi(φ)*sinpi(θ) for θ in θ, φ in φ]
z = [cospi(θ) for θ in θ, φ in φ]

f2(x, y, z) = cross([x, y, z], [1.0, 0.0, 0.0])
tans = f2.(vec(x), vec(y), vec(z))
u = [a[1] for a in tans]
v = [a[2] for a in tans]
w = [a[3] for a in tans]

#scene = Scene();
        fig = Figure(; resolution=(600, 480))
        ax1 = Axis3(fig[1, 1]; aspect=(1,1,1), perspectiveness=0.5)
surface!(ax1, x, y, z)
arr = GLMakie.arrows!(ax1, 
    vec(x), vec(y), vec(z), u, v, w;
    arrowsize = 0.1, linecolor = (:gray, 0.7), linewidth = 0.02, lengthscale = 0.1
)
display(fig)