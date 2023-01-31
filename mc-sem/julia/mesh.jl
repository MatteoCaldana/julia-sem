# # # flattend meshgrid of hyper-rectangle with opposite vertexes a and b
# # # on each direction i it is subsivided in ns[i] - 1 elements
# # function flattened_meshgrid(a::Vector{Float64}, b::Vector{Float64}, ns::Vector{Int64})::Matrix{Float64}
# #   @assert length(a) == length(b)
# #   @assert length(a) == length(ns)
# #   nprod = ones(Int64, length(ns) + 1)
# #   dt = Array{Float64}(undef, size(ns))
# #   for i = 1:length(ns)
# #     nprod[i + 1] = nprod[i] * ns[i]
# #     dt[i] = (b[i] - a[i]) / (ns[i] - 1)
# #   end
# #   mgrid = Array{Float64}(undef, size(ns, 1), nprod[end])
# #   @inbounds for i=1:nprod[end]
# #     mgrid[1, i] = ((i - 1) % nprod[2]) * dt[1] + a[1]
# #     @inbounds for j = 2:size(ns, 1)
# #       mgrid[j, i] = div(i - 1, nprod[j]) * dt[j] + a[j]
# #     end
# #   end

# #   return mgrid
# # end

# # function hyper_rectangle_edges(ns::Vector{Int64})::Matrix{Int64}
# #   nprod = ones(Int64, length(ns) + 1)
# #   for i = 1:length(ns)
# #     nprod[i + 1] = nprod[i] * ns[i]
# #   end

# #   nedges = 0
# #   for dim = 1:length(ns)
# #     nedges += (ns[dim] - 1) * div(nprod[end], ns[dim])
# #   end

# #   edges = Array{Int64}(undef, 2, nedges)
# #   filled = 0
# #   @inbounds for dim = 1:length(ns)
# #     a = vcat([1:nprod[dim]:(nprod[dim+1] - nprod[dim]);]', [1+nprod[dim]:nprod[dim]:nprod[dim+1];]')
# #     @inbounds for i = 1:div(nprod[end], ns[dim])
# #       edges[:, filled + (i - 1)*size(a, 2) + 1: filled + i*size(a, 2)] = a .+ (((i - 1) % nprod[dim]) + div(i - 1,  nprod[dim]) * nprod[dim+1])
# #     end
# #     filled += (ns[dim] - 1) * div(nprod[end], ns[dim])
# #   end
# #   return edges
# # end

# # function hyper_rectangle_faces(ns::Vector{Int64})::Matrix{Int64}
# #   nprod = ones(Int64, length(ns) + 1)
# #   for i = 1:length(ns)
# #     nprod[i + 1] = nprod[i] * ns[i]
# #   end

# #   nfaces = 0
# #   for dim = 1:length(ns)
# #     dim2 = dim % length(ns) + 1
# #     nfaces += (ns[dim] - 1) * (ns[dim2] - 1) * div(nprod[end], ns[dim]*ns[dim2])
# #   end

# #   faces = Array{Int64}(undef, 4, nfaces)
# #   filled = 0
# #   @inbounds for dim = 1:length(ns)
# #     a = vcat([1:nprod[dim]:(nprod[dim+1] - nprod[dim]);]', [1+nprod[dim]:nprod[dim]:nprod[dim+1];]')
# #     @inbounds for i = 1:div(nprod[end], ns[dim])
# #       edges[:, filled + (i - 1)*size(a, 2) + 1: filled + i*size(a, 2)] = a .+ (((i - 1) % nprod[dim]) + div(i - 1,  nprod[dim]) * nprod[dim+1])
# #     end
# #     filled += (ns[dim] - 1) * div(nprod[end], ns[dim])
# #   end
# #   return edges
# # end

#println(flattened_meshgrid([0., 0., 0.], [1., 1., 1.], [3, 4, 5]))
#println(hyper_rectangle_edges([3, 4, 5]))
#hyper_rectangle_faces([3, 4, 5])

###############################################################################
###############################################################################

using LinearAlgebra

abstract type Mesh
end

struct CartesianMesh <: Mesh
  nelem::Int64
  nvert::Int64
  dim::Int64
  ns::Vector{Int64}
  nprod::Vector{Int64}
  nprodv::Vector{Int64}

  function CartesianMesh(ns::Vector{Int64})
    dim = length(ns)
    nprod = ones(Int64, dim + 1)
    nprodv = ones(Int64, dim + 1)
    nvert = 1
    for i = 1:dim
      nprod[i + 1] = nprod[i] * ns[i]
      nprodv[i + 1] = nprodv[i] * (ns[i] + 1)
      nvert *= ns[i] + 1
    end
    return new(nprod[end], nprodv[end], dim, ns, nprod, nprodv)
  end
end

function tuple_idx(prod::Vector{Int64}, id::Int64)
  dim = length(prod) - 1
  tid = Vector{Int64}(undef, dim)
  for j = 1:dim
    tid[j] = div((id - 1) % prod[j + 1], prod[j])
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
  println("eid: ", eid)
  println("type: ", type)
  println("id: ", id)

  id1 = id - 1
  if type == 4 #"element"
    return eid
  else 
    idx = tuple_idx(mesh.nprod, eid)
    println("idx: ", idx)
    if type == 3 # "face"
      @assert id <= 6
      if id > 3
        idx[id - 3] += 1
      end
      eid = flat_idx(mesh.nprodv, idx)
      return 3*eid + id1
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
        [0, 1, 1, 0];; 
      ]
      println(size(TABLE))
      println("tab:", TABLE[1:3, id], id)
      idx += TABLE[1:3, id]
      eid = flat_idx(mesh.nprodv, idx)
      println("idx: ", idx)
      println("eid: ", eid)
      return 3*eid + TABLE[4, id]
    elseif type == 1 #"vertex"
      @assert id <= 8
      idx += [id1 & 1, (id1 & 2) > 0, (id1 & 4) > 0]
      println("idx: ", idx)
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
  prod = degp.^[0:dim;]
  idx = tuple_idx(prod, id)

  deltas = hcat(I(dim), -I(dim))
  neighbours = []

  for i = 1:size(deltas, 2)
    candidate_neighbour = idx + deltas[:, i]
    if in_range(candidate_neighbour, degp)
      #println(id, idx, candidate_neighbour, flat_idx(prod, candidate_neighbour))
      push!(neighbours, flat_idx(prod, candidate_neighbour))
    end
  end

  return neighbours
end

function get_reference_dof(deg::Int64, dim::Int64)::Matrix{Int64}
  degp = deg + 1
  ndofs = degp^dim
  dof_1d = ones(Int64, degp)
  dof_1d[1] = 2
  dof_1d[end] = 2

  # 2^dim   == vertexes
  # 2^dim-1 == edges
  # 2^dim-2 == faces
  # 2^dim-3 == cubes
  # ...
  cart_prod = 1
  for i = 1:dim
    shape = ones(Int64, dim)
    shape[i] = degp
    cart_prod = cart_prod .* reshape(dof_1d, tuple(shape...))
  end

  cart_prod_f = reshape(cart_prod, ndofs)
  # id for each component (vertex, edge, ...), different components have different ids
  ids = zeros(Int64, ndofs)

  curr_max_ids = ones(Int64, 2^dim) # need just dim+1 but too lazy to compute log
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

  return hcat((dim - 62) .+ leading_zeros.(cart_prod_f), ids)
end

function distribute_dof(grid::Mesh, deg::Int64)::Matrix{Int64}
  degp = deg + 1
  reference_dof = get_reference_dof(deg, grid.dim)

  dof_map = zeros(Int64, grid.nelem, size(reference_dof, 1))
  # lookup table to check if already enumerated
  table = Array{Vector{Int64}}(undef, grid.dim + 1, 4*grid.nvert+10)
  curr_dof_id = 1

  for i = 1:size(dof_map, 1)
    for j = 1:size(dof_map, 2)
      type = reference_dof[j, 1]
      id = reference_dof[j, 2]
      comp = get_component_id(grid, i, type, id)
      println("comp: ", comp)
      
      if dof_map[i, j] == 0
        println("not zero")
        if isassigned(table, type, comp) && length(table[type, comp]) == (deg - 1)^(type - 1)
          println("well def")
          println(type, " ", id, " ", (deg - 1)^(type - 1))
          
          println(table[type, comp])
          idx = 1
          for k = 1:size(dof_map, 2)
            if (reference_dof[k, 1] == type) && (reference_dof[k, 2] == id)
              dof_map[i, k] = table[type, comp][idx]
              idx += 1
            end
          end
        else
          if !isassigned(table, type, comp)
            println("not at all")
            table[type, comp] = zeros(Int64, 0)
          end
          dof_map[i, j] = curr_dof_id
          push!(table[type, comp], curr_dof_id)
          curr_dof_id += 1
        end
      end
      println(dof_map[i, :])
      println("----------------------------------------") 
    end
    println("========================================")
    println("========================================")
  end
  println(size(dof_map))
  return dof_map
end


#get_reference_dof(3, 3)
mesh = CartesianMesh([3, 1, 1])
dof_map = distribute_dof(mesh, 3)
println(dof_map)
