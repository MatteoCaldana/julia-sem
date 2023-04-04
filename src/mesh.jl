include("core.jl")

using LinearAlgebra
using BenchmarkTools

abstract type Mesh end

struct CartesianMesh <: Mesh
  nelem::Int64
  nvert::Int64
  dim::Int64
  ns::Vector{Int64}
  nprod::Vector{Int64}
  nprodv::Vector{Int64}
  corner_a::Vector{Float64}
  corner_b::Vector{Float64}
  jd::Vector{Float64}

  function CartesianMesh(a::Vector{Float64}, b::Vector{Float64}, ns::Vector{Int64})
    @assert length(a) == length(b)
    @assert length(a) == length(ns)
    dim = length(ns)
    nprod = ones(Int64, dim + 1)
    nprodv = ones(Int64, dim + 1)
    nvert = 1
    for i = 1:dim
      nprod[i+1] = nprod[i] * ns[i]
      nprodv[i+1] = nprodv[i] * (ns[i] + 1)
      nvert *= ns[i] + 1
    end
    h = (b - a) ./ ns
    return new(nprod[end], nprodv[end], dim, ns, nprod, nprodv, a, b, h ./ 2)
  end
end

function tuple_idx(prod::Vector{Int64}, id::Int64)
  dim = length(prod) - 1
  tid = Vector{Int64}(undef, dim)
  for j = 1:dim
    @inbounds tid[j] = div((id - 1) % prod[j+1], prod[j])
  end
  return tid
end

function flat_idx(prod::Vector{Int64}, id::Vector{Int64})
  fid = id[1]
  for j = 2:length(prod)-1
    @inbounds fid += id[j] * prod[j]
  end
  return fid + 1
end

# give the index of the mesh element, the type of the component, i.e.
#  .1 == vertex
#  .2 == edge
#  .3 == face
#  ...  
# and the id of the component in the local enumeration 
# return the global id of the component for the mesh
#       y
#       |
#       |
#       ---------x
#      /
#     /
#    z
# first three indexes indicates the cartesian move of the element 
# the fourth indicates the new id
const EDGE_ID_TABLE = [
  [0, 0, 0, 0];; # x bottom far
  [0, 0, 0, 1];; # y left far 
  [1, 0, 0, 1];; # y right far
  [0, 1, 0, 0];; # x top far
  [0, 0, 0, 2];; # z left bottomp
  [1, 0, 0, 2];; # z right bottom
  [0, 1, 0, 2];; # z left top
  [1, 1, 0, 2];; # z right top
  [0, 0, 1, 0];; # x bottom near
  [0, 0, 1, 1];; # y left near
  [1, 0, 1, 1];; # y right near
  [0, 1, 1, 0]   # x top near
]
const POW2 = [1, 2, 4]
function get_component_id(mesh::CartesianMesh, eid::Int64, type::Int64, id::Int64)
  id1 = id - 1
  if type == 4 #"element"
    return eid
  else
    idx = tuple_idx(mesh.nprod, eid)
    if type == 3 # "face"
      # face with id is:
      # 1: xy far
      # 2: xz bottom
      # 3: yz left
      # 4: yz right
      # 5: xz top
      # 6: xy near
      # @assert id <= 6
      if id > 3
        idx[id-3] += 1
        id1 = 5 - id1
      end
      eid = flat_idx(mesh.nprodv, idx)
      return 3 * eid + id1
    elseif type == 2 # "edge"
      # @assert id <= 12
      idx += EDGE_ID_TABLE[1:3, id]
      eid = flat_idx(mesh.nprodv, idx)
      return 3 * eid + EDGE_ID_TABLE[4, id]
    elseif type == 1 #"vertex"
      # @assert id <= 8
      for i = 1:3
        @inbounds idx[i] += (id1 & POW2[i]) > 0
      end
      return flat_idx(mesh.nprodv, idx)
    end
  end
end

function in_range(v::Vector{Int64}, mmax::Int64)
  for j = 1:length(v)
    if v[j] < 0 || v[j] >= mmax
      return false
    end
  end
  return true
end

function find_neighbours(id::Int64, dim::Int64, degp::Int64)
  prod = degp .^ [0:dim;]
  idx = tuple_idx(prod, id)

  deltas = hcat(I(dim), -I(dim))
  neighbours = []

  for i in axes(deltas, 2)
    candidate_neighbour = idx + deltas[:, i]
    if in_range(candidate_neighbour, degp)
      push!(neighbours, flat_idx(prod, candidate_neighbour))
    end
  end

  return neighbours
end

function get_reference_dof(deg::Int64, dim::Int64)::Matrix{Int64}
  if deg == 1
    return hcat(ones(Int64, 8), [1:8;])
  end
  degp = deg + 1
  ndofs = degp^dim

  # 2^dim   == vertexes
  # 2^dim-1 == edges
  # 2^dim-2 == faces
  # 2^dim-3 == cubes
  # ...
  dof_1d = [2, ones(Int64, deg - 1)..., 2]
  cart_prod = map(prod, Base.product(ntuple(x -> dof_1d, dim)...))
  cart_prod_f = (dim - 62) .+ leading_zeros.(reshape(cart_prod, ndofs))

  # id for each component (vertex, edge, ...), different components have different ids
  ids = zeros(Int64, ndofs)
  curr_max_ids = ones(Int64, dim + 1)
  for i = 1:ndofs
    neighbours = find_neighbours(i, dim, degp)
    found = false
    for j in eachindex(neighbours)
      n = neighbours[j]
      if cart_prod_f[n] == cart_prod_f[i] && ids[n] != 0
        found = true
        ids[i] = ids[n]
        break
      end
    end
    if !found
      ids[i] = curr_max_ids[cart_prod_f[i]]
      curr_max_ids[cart_prod_f[i]] += 1
    end
  end

  return hcat(cart_prod_f, ids)
end

function element_transform(mesh::CartesianMesh, id::Int64, pt::Vector{Float64})::Vector{Float64}
  el_idx = tuple_idx(mesh.nprod, id)
  return mesh.jd .* (pt + 2 * el_idx .+ 1.0) + mesh.corner_a
end

function distribute_dof(grid::Mesh, deg::Int64)::Tuple{Matrix{Int64},Matrix{Float64}}
  degp = deg + 1
  reference_dof = get_reference_dof(deg, grid.dim)
  rx, _ = xwlgl(degp)
  rx = Float64.(rx)

  dof_map = zeros(Int64, size(reference_dof, 1), grid.nelem)
  dof_support = Vector{Vector{Float64}}(undef, 0)
  # lookup table to check if already enumerated
  table = Array{Vector{Int64}}(undef, grid.dim + 1, 4 * grid.nvert + 10)
  curr_dof_id = 1

  for i in axes(dof_map, 2)
    for j in axes(dof_map, 1)
      type = reference_dof[j, 1]
      id = reference_dof[j, 2]
      comp = get_component_id(grid, i, type, id)

      if dof_map[j, i] == 0
        if isassigned(table, type, comp) && length(table[type, comp]) == (deg - 1)^(type - 1)
          idx = 1
          for k in axes(dof_map, 1)
            if (reference_dof[k, 1] == type) && (reference_dof[k, 2] == id)
              dof_map[k, i] = table[type, comp][idx]
              idx += 1
            end
          end
        else
          if !isassigned(table, type, comp)
            table[type, comp] = zeros(Int64, 0)
          end
          dof_map[j, i] = curr_dof_id
          push!(table[type, comp], curr_dof_id)
          push!(dof_support, element_transform(grid, i, rx[tuple_idx(degp .^ [0:grid.dim;], j).+1]))
          curr_dof_id += 1
        end
      end
    end
  end
  return dof_map, hcat(dof_support...)
end

function find_bc(mesh::CartesianMesh, pts)
  EPS = 1e-15
  is_bc = sum((abs.(pts .- mesh.corner_a) .< EPS) + (abs.(pts .- mesh.corner_b) .< EPS), dims=1)
  return findall(x -> x > 0, vec(is_bc))
end

# mesh = CartesianMesh([-1., -1., -1.], [1., 1., 1.], [2, 2, 2])
# dof_map, dof_support = distribute_dof(mesh, 2)
# display(dof_map')
# println(size(dof_support), dof_support)
