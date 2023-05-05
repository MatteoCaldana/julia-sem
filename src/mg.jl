include("mesh.jl")
include("core.jl")

struct MGLevel
  refel::FESpace
  mesh::CartesianMesh

  has_coarse::Bool
  coarse::Union{MGLevel, Missing}

  dof_map::Matrix{Int64}
  dof_map_to_fine::Matrix{Int64}
  dof_support::Matrix{Float64}
  boundary_idxs::Vector{Int64}
  boundary_values::Vector{Float64}

  Ph::Matrix{Float64}
  Pp::Matrix{Float64}

  ndofs::Int64
  ndofs_fine::Int64

  function MGLevel(deg, mesh, bc_fun)
    fe = FESpace(deg)

    has_coarse = sum(mesh.ns .% 2) == 0
    if has_coarse
      coarse_mesh = CartesianMesh(mesh.corner_a, mesh.corner_b, mesh.ns .÷ 2)
      coarse = MGLevel(deg, coarse_mesh, bc_fun)
    else
      coarse = missing
    end
    dof_map, dof_support = distribute_dof_v2(mesh, deg)
    dof_map_to_fine = projection_element_mapping(mesh, deg)
    boundary_idxs = find_bc(mesh, dof_support)
    boundary_values = bc_fun(dof_support[:, boundary_idxs])

    ndofs = prod(mesh.ns .* deg .+ 1)
    ndofs_fine = prod(2 .* mesh.ns .* deg .+ 1)

    return new(
      fe, 
      mesh, 
      has_coarse, coarse, 
      dof_map, dof_map_to_fine, dof_support, 
      boundary_idxs, boundary_values,
      reduce(kron, Iterators.repeated(fe.Ph, mesh.dim)),
      reduce(kron, Iterators.repeated(fe.Pp, mesh.dim)),
      ndofs, ndofs_fine
    )
  end
end

# from this level to the finer one
# TODO: check this matches homg
function interpolation_vmult(lvl::MGLevel, x::Vector{Float64})::Vector{Float64}
  y = zeros(Float64, lvl.ndofs_fine)

  for ie in axes(lvl.dof_map_to_fine, 2)
    for j in axes(lvl.dof_map_to_fine, 1)
      y[lvl.dof_map_to_fine[:, ie]] += lvl.Ph * x[lvl.dof_map[:, ie]] 
    end 
  end

  return y
end

# TODO: check this is equal to the transpose of the above
function projection_vmult(lvl::MGLevel, x::Vector{Float64})::Vector{Float64}
  y = zeros(Float64, lvl.ndofs)

  for ie in axes(lvl.dof_map, 2)
    for j in axes(lvl.dof_map, 1)
      y[lvl.dof_map[:, ie]] += transpose(lvl.Ph) * x[lvl.dof_map_to_fine[:, ie]] 
    end 
  end

  return y
end

nel = 4
mesh = CartesianMesh([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], nel*[1, 1, 1])
mglvl = MGLevel(1, mesh, xyz -> zeros(size(xyz, 2)))

interpolation_vmult(mglvl, ones(mglvl.ndofs))