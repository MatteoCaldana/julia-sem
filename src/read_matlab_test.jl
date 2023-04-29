using MAT 
using LinearAlgebra

include("mesh.jl")

function test_linear_solver()
  nel = 4
  deg = 3
  vars = matread("../test-data/3d_stiff_deg"*repr(deg)*"_nel"*repr(nel)*".mat")

  A = vars["A"]
  f = vars["f"]
  u = vars["un"]

  u_julia = A \ f
  e = u - u_julia

  println(norm(e, 2))
  println(norm(e, Inf))
end

function test_dof_map(nel, deg)
  vars = matread("../test-data/dofmap_idx.deg"*repr(deg)*".nel"*repr(nel)*".mat")
  dof_map_mat = Int64.(vars["idx"])

  vars = matread("../test-data/projec_idx.deg"*repr(deg)*".nel"*repr(nel)*".mat")
  dof_map_fine_mat = Int64.(vars["idx_f"])


  mesh = CartesianMesh([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], nel*[1, 1, 1])
  dof_map, _ = distribute_dof_v2(mesh, deg)
  dof_map_fine = projection_element_mapping(mesh, deg)

  err = maximum(abs.(dof_map_mat - dof_map))
  err_f = maximum(abs.(dof_map_fine_mat - dof_map_fine))
  println(err, " ", err_f)
end

for nel = [1, 2, 4, 8]
  for deg = [1, 2, 3, 4]
    test_dof_map(nel, deg)
  end
end
