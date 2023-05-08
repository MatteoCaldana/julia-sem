include("mesh.jl")
include("core.jl")
include("matrix_free_utils.jl")

using SparseArrays 

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
  P::SparseMatrixCSC{Float64, Int64}

  A::Matrix{Float64}
  Ad::Vector{Float64}
  M::Vector{Float64}

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

    Ph = reduce(kron, Iterators.repeated(fe.Ph, mesh.dim))
    Pp = reduce(kron, Iterators.repeated(fe.Pp, mesh.dim))

    A = Float64.(stiff_3d_sp(fe.w_128, fe.D_128, mesh.jd[1], mesh.jd[2], mesh.jd[3]))
    M = Float64.(reshape(map(prod, Base.product(ntuple(x -> fe.w_128, mesh.dim)...)), length(fe.w)^mesh.dim))
    M *= mesh.jd[1] * mesh.jd[2] * mesh.jd[3]

    ndofs = prod(mesh.ns .* deg .+ 1)
    ndofs_fine = prod(2 .* mesh.ns .* deg .+ 1)

    NP_c = (deg + 1)^mesh.dim
    NP_f = (2*deg + 1)^mesh.dim
    ne = prod(mesh.ns)
    NPNP = NP_c * NP_f;

    I = zeros(Int64, ne * NPNP);
    J = zeros(Int64, ne * NPNP);
    val = zeros(Float64, ne * NPNP);

    for ie=1:ne
      idx_c = dof_map[:, ie]
      idx_f = dof_map_to_fine[:, ie]

      ind1 = repeat(idx_f,NP_c,1);
      ind2 = reshape(repeat(idx_c',NP_f,1),NPNP,1);
      st = (ie-1)*NPNP+1;
      en = ie*NPNP;
      I[st:en] = ind1;
      J[st:en] = ind2;
      val[st:en] = Ph[:];
    end

    return new(
      fe, 
      mesh, 
      has_coarse, coarse, 
      dof_map, dof_map_to_fine, dof_support, 
      boundary_idxs, boundary_values,
      Ph, Pp,
      sparse(I,J,val,ndofs_fine,ndofs, (x,y)->max(x, y)),
      A, diag(A), M,
      ndofs, ndofs_fine
    )
  end
end

# from this level to the finer one
function interpolation_vmult(lvl::MGLevel, x::Vector{Float64})::Vector{Float64}
  y = Vector{Float64}(undef, lvl.ndofs_fine)
  for ie in axes(lvl.dof_map_to_fine, 2)
    y[lvl.dof_map_to_fine[:, ie]] = lvl.Ph * x[lvl.dof_map[:, ie]] 
  end
  return y
end

function projection_vmult(lvl::MGLevel, x::Vector{Float64})::Vector{Float64}
  y = zeros(lvl.ndofs)
  for ie in axes(lvl.dof_map, 2)
    # TODO: is this worth it? 
    #       would it be better to materialize the projection/interpolation
    #       considering this extra cost and the fact this is a dirty trick?
    y[lvl.dof_map[:, ie]] += Float64.(y[lvl.dof_map[:, ie]] .== 0) .* (transpose(lvl.Ph) * x[lvl.dof_map_to_fine[:, ie]]) 
  end
  return y
end

function test_transpose()
  nel = 2
  for deg = [1, 2, 3]
    mesh = CartesianMesh([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], nel*[1, 1, 1])
    bc = xyz -> zeros(size(xyz, 2))
    mglvl = MGLevel(deg, mesh, bc)

    lvl = mglvl
    while true
      I = materialize_linear_map(x -> interpolation_vmult(lvl, x), lvl.ndofs_fine, lvl.ndofs)
      P = materialize_linear_map(x -> projection_vmult(lvl, x),    lvl.ndofs, lvl.ndofs_fine)
      println("deg: ", deg, " ndofs:", lvl.ndofs, " ", all(P.==I'))
      lvl.has_coarse || break
      lvl = lvl.coarse
    end
  end
end

#test_transpose()