include("core.jl")
include("mesh.jl")

using MAT 


struct FESpace
  n::Int64
  np::Int64
  x::Vector{Float64}
  w::Vector{Float64}
  d::Matrix{Float64}

  x_128::Vector{Float128}
  w_128::Vector{Float128}
  d_128::Matrix{Float128}

  function FESpace(n::Int64)
    x, w = xwlgl(n + 1)
    d = derlgl(x)
    return new(n, n + 1, Float64.(x), Float64.(w), Float64.(d), x, w, d)
  end
end

function assemble(fe::FESpace, mesh::Mesh, dof_map)
  Aloc = Float64.(stiff_3d_sp(fe.w_128, fe.d_128, mesh.dh[1] / 2, mesh.dh[2] / 2, mesh.dh[3] / 2))
  Mloc = Float64.(reshape(map(prod, Base.product(ntuple(x -> fe.w_128, mesh.dim)...)), length(fe.w)^mesh.dim))
  Mloc *= mesh.dh[1] * mesh.dh[2] * mesh.dh[3] / 8
  ndof = dof_map[end, end]
  A = zeros(Float64, ndof, ndof)
  M = zeros(Float64, ndof)
  for i = 1:size(dof_map, 1)
    A[dof_map[i, :], dof_map[i, :]] += Aloc
    M[dof_map[i, :]] += Mloc
  end
  return A, M
end

function apply_dirichlet_bc(A, x, b, bc, sol, dof_support)
  for idx = 1:length(bc)
    i = bc[idx]
    A[i, :] .= 0
    A[i, i] = 1
    v = sol(dof_support[:, i:i])[1]
    x[i] = v
    b[i] = v
  end
end

function cg(Avmult, b, x, tol)
  r = b - Avmult(x)
  p = r
  rsold = r' * r
  for i = 1:length(b)
    Ap = Avmult(p)
    alpha = rsold / (p' * Ap)
    x += alpha * p
    r -= alpha * Ap
    rsnew = r' * r
    if sqrt(rsnew) < tol
      rsold = rsnew
      break
    end
    p = r + (rsnew / rsold) * p
    rsold = rsnew
  end
  println("Residual:", rsold)
  return x
end

function compute_err(u, un, dun, fe, mesh, dof_map)
  dim = Int64(round(log(size(dof_map, 2)) / log(length(fe.w)), RoundNearest))
  ww = map(prod, Base.product(ntuple(x -> fe.w, dim)...))
  ww = reshape(ww, length(ww))
  l2_err = 0
  h1_err = 0
  for i = 1:size(dof_map, 1)
    li = dof_map[i, :]
    # there should be a jacobian of the element here, 
    # so this is the norm up to a constant (if all elements have the same jacobian)
    l2_err += sum((u[li] - un[li]) .^ 2 .* ww)

    u_loc = reshape(u[li], ntuple(x -> length(fe.w), dim))
    ux = Array{Float64}(undef, size(u_loc))
    uy = Array{Float64}(undef, size(u_loc))
    uz = Array{Float64}(undef, size(u_loc))
    for i = 1:length(fe.w)
      for j = 1:length(fe.w)
        ux[:, i, j] = fe.d * u_loc[:, i, j] / mesh.dh[1] / 2
        uy[i, :, j] = fe.d * u_loc[i, :, j] / mesh.dh[2] / 2
        uz[i, j, :] = fe.d * u_loc[i, j, :] / mesh.dh[3] / 2
      end
    end
    ux = reshape(ux, length(ux))
    uy = reshape(uy, length(uy))
    uz = reshape(uz, length(uz))
    h1_err += sum(((ux - dun[li, 1]) .^ 2 + (uy - dun[li, 2]) .^ 2 + (uz - dun[li, 3]) .^ 2) .* ww)
  end
  return sqrt(l2_err), sqrt(h1_err)
end

# function force(xyz)
#   return 6 * xyz[1, :] - (xyz[2, :] .^ 2 + xyz[3, :] .^ 2) .* sin.(xyz[2, :] .* xyz[3, :])
# end

# function sol(xyz)
#   return sin.(xyz[2, :] .* xyz[3, :]) + xyz[1, :] .^ 3
# end

function force(xyz)
  return 12 * π * π * sin.(2 * π * xyz[1, :]) .* cos.(2 * π * xyz[2, :]) .* sin.(2 * π * xyz[3, :])
end

function sol(xyz)
  return sin.(2 * π * xyz[1, :]) .* cos.(2 * π * xyz[2, :]) .* sin.(2 * π * xyz[3, :])
end

# for deg = 1:4
#   for nel = [1, 2, 4, 8, 16, 32, 64]
#     println("Basis deg: ", deg)
#     println("#elements: ", nel)
#     println("Building FE space")
#     fespace = FESpace(deg)
#     println("Building mesh")
#     mesh = CartesianMesh([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], nel*[1, 1, 1])
#     println("Distributing dof")
#     dof_map, dof_support = distribute_dof(mesh, deg)
#     println("Problem has ", dof_map[end, end], " dofs")
#     println("Assembling system")
#     A = assemble(fespace, mesh, dof_map)
#     f = force(dof_support)
#     println("Applying BC")
#     x0 = zeros(Float64, size(f))
#     bc = find_bc(mesh, dof_support)
#     apply_dirichlet_bc(A, x0, f, bc, sol, dof_support)
#     println("Solving system")
#     u = A \ f

#     println("Compunting error")
#     u_ex = sol(dof_support)
#     du_ex = zeros(Float64, size(u_ex, 1), 3)
#     println("Res:",  norm(A*u - f) / norm(f))
#     println("Err:",  norm(u - u_ex))
#     err_l2, err_h1 = compute_err(u, u_ex, du_ex, fespace, mesh, dof_map)
#     println("Err L2:", err_l2)
#     println("Err H1:", err_h1)
#     println("-------------------------------------")
#   end
# end

function loo(x, y)
  return maximum(abs.(x - y))  
end

for deg = 1:4
  for nel = 1:4
    println("Degree: ", deg)
    println("#Elements: ", nel)
    vars = matread("test-data/3d_stiff_deg"*repr(deg)*"_nel"*repr(nel)*"_v2.mat")

    fespace = FESpace(deg)
    mesh = CartesianMesh([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], nel * [1, 1, 1])
    dof_map, dof_support = distribute_dof(mesh, deg) 

    println(loo(dof_support', vars["xyz"]))
    println(loo(dof_map', vars["nov"]))

    A, M = assemble(fespace, mesh, dof_map)
    f = force(dof_support) .* M

    println("Test A: ", loo(A, vars["A_org"]))
    println("Test f: ", loo(f, vars["f_org"]))

    bc = find_bc(mesh, dof_support)
    x0 = zeros(Float64, size(f))
    apply_dirichlet_bc(A, x0, f, bc, sol, dof_support)
    println("Test w bc A: ", loo(A, vars["A"]))
    println("Test w bc f: ", loo(f, vars["f"]))

    u = A \ f
    println("Test u: ", loo(u, vars["un"]))
    u_ex = sol(dof_support)

    println("Res:",  norm(A*u - f) / norm(f))
    println("Err:",  norm(u - u_ex) / norm(u_ex), " (normalized)")

    println("Err       :",  norm(u - u_ex))
    println("Err MATLAB:",  norm(vars["un"] - vars["u_ex"]))
    println("------------------------------------------------")
  end
end

# b = @benchmark cg(x->A*x, f, x0, 1e-10)
# print_benckmark(b)