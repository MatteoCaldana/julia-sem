include("core.jl")
include("mesh.jl")

using PProf       # @pprof
using ProfileView # @profview

deg = 1
nel = 128
mesh = CartesianMesh([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], [1, 1, 1])
dof_map, dof_support = distribute_dof(mesh, deg)

for i = [2, 4, 8, 16, 32, 64]
  mesh = CartesianMesh([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], i*[1, 1, 1])
  @time distribute_dof(mesh, deg)
end
# @pprof distribute_dof(mesh, deg)
# read(stdin, Char)