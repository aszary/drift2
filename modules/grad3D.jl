using CoupledFields
using GLMakie
using GeometryBasics
GLMakie.activate!()

sz = 1000

g(x,y,z) = x .* exp.(-x.^2 - y.^2 - z.^2)
X = -2 .+ 4*rand(sz, 3)
Y = g.(X[:,1], X[:,2], X[:,3])

"""
println(size(X))
println(size(Y))
println(length(Y))
println(size(Y[:,1:1]))
println(typeof(Y))
#println(Y[:,1:1])
#println(Y)
"""

kernelpars = GaussianKP(X)
∇g = gradvecfield([0.5 -7], X, Y[:,1:1], kernelpars)
println(size(∇g[1]))

ngx = Array{Float64}(undef, sz)
ngy = Array{Float64}(undef, sz)
ngz = Array{Float64}(undef, sz)

for i in 1:sz
        ngx[i] = ∇g[i][1]
        ngy[i] = ∇g[i][2]
        ngz[i] = ∇g[i][3]
end


fig, ax1, p = mesh(Sphere(Point3f(0, 0, 0), 0.03), color=:blue, transparency=true) # better camera control (Scene), but zlims does not work
#println(X[:, 1])
#fig = Figure(; resolution=(600, 480))
#ax1 = Axis3(fig[1, 1]; aspect=(1,1,1), perspectiveness=0.5)
#scatter!(ax1,X[:, 1] ,X[:, 2] ,X[:, 3] , marker=:diamond, color=:black)
hm = meshscatter!(ax1, X[:, 1], X[:, 2], X[:, 3]; markersize=0.05, color=Y, transparency=false)
arrows!(ax1, X[:, 1], X[:, 2], X[:, 3],  ngx, ngy, ngz, color=:white, arrowsize=Vec3f(0.08, 0.09, 0.1), linewidth=0.05)
display(fig)
