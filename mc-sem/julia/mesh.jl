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
  dh::Vector{Float64}

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
    return new(nprod[end], nprodv[end], dim, ns, nprod, nprodv, a, b, h)
  end
end

function tuple_idx(prod::Vector{Int64}, id::Int64)
  dim = length(prod) - 1
  tid = Vector{Int64}(undef, dim)
  for j = 1:dim
    tid[j] = div((id - 1) % prod[j+1], prod[j])
  end
  return tid
end

function flat_idx(prod::Vector{Int64}, id::Vector{Int64})
  fid = id[1]
  for j = 2:length(prod)-1
    fid += id[j] * prod[j]
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
function get_component_id(mesh::CartesianMesh, eid::Int64, type::Int64, id::Int64)
  id1 = id - 1
  if type == 4 #"element"
    return eid
  else
    idx = tuple_idx(mesh.nprod, eid)
    if type == 3 # "face"
      @assert id <= 6
      if id > 3
        idx[id-3] += 1
        id1 -= 3
      end
      eid = flat_idx(mesh.nprodv, idx)
      return 3 * eid + id1
    elseif type == 2 # "edge"
      @assert id <= 12
      TABLE = [
        [0, 0, 0, 0];;
        [0, 0, 0, 1];;
        [1, 0, 0, 1];;
        [0, 1, 0, 0];;
        [0, 0, 0, 2];;
        [1, 0, 0, 2];;
        [0, 1, 0, 2];;
        [1, 1, 0, 2];;
        [0, 0, 1, 0];;
        [0, 0, 1, 1];;
        [1, 0, 1, 1];;
        [0, 1, 1, 0]
      ]
      idx += TABLE[1:3, id]
      eid = flat_idx(mesh.nprodv, idx)
      return 3 * eid + TABLE[4, id]
    elseif type == 1 #"vertex"
      @assert id <= 8
      idx += [id1 & 1, (id1 & 2) > 0, (id1 & 4) > 0]
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

  for i = 1:size(deltas, 2)
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
  cart_prod = map(prod, Base.product(ntuple(x->dof_1d, dim)...))
  cart_prod_f = (dim - 62) .+ leading_zeros.(reshape(cart_prod, ndofs))

  # id for each component (vertex, edge, ...), different components have different ids
  ids = zeros(Int64, ndofs)
  curr_max_ids = ones(Int64, dim+1)
  for i = 1:ndofs
    neighbours = find_neighbours(i, dim, degp)
    found = false
    for j = 1:length(neighbours)
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
  return mesh.dh .* ( pt ./ 2 + el_idx .+ 0.5) + mesh.corner_a
end

function distribute_dof(grid::Mesh, deg::Int64)::Tuple{Matrix{Int64},Matrix{Float64}}
  degp = deg + 1
  reference_dof = get_reference_dof(deg, grid.dim)
  rx, _ = xwlgl(degp)
  rx = Float64.(rx)

  dof_map = zeros(Int64, grid.nelem, size(reference_dof, 1))
  dof_support = Vector{Vector{Float64}}(undef, 0)
  # lookup table to check if already enumerated
  table = Array{Vector{Int64}}(undef, grid.dim + 1, 4 * grid.nvert + 10)
  curr_dof_id = 1

  for i = 1:size(dof_map, 1)
    for j = 1:size(dof_map, 2)
      type = reference_dof[j, 1]
      id = reference_dof[j, 2]
      comp = get_component_id(grid, i, type, id)

      if dof_map[i, j] == 0
        if isassigned(table, type, comp) && length(table[type, comp]) == (deg - 1)^(type - 1)
          idx = 1
          for k = 1:size(dof_map, 2)
            if (reference_dof[k, 1] == type) && (reference_dof[k, 2] == id)
              dof_map[i, k] = table[type, comp][idx]
              idx += 1
            end
          end
        else
          if !isassigned(table, type, comp)
            table[type, comp] = zeros(Int64, 0)
          end
          dof_map[i, j] = curr_dof_id
          push!(table[type, comp], curr_dof_id)
          push!(dof_support, element_transform(grid, i, rx[tuple_idx(degp.^[0:grid.dim;], j) .+ 1]))
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

#get_reference_dof(3, 3)
# mesh = CartesianMesh([-1., -1., -1.], [1., 1., 1.], [1, 1, 2])
# dof_map, dof_support = distribute_dof(mesh, 3)
# println(dof_map)
# println(size(dof_support), dof_support)
