include("core.jl")
include("mesh.jl")
include("csr.jl")

using MAT 
#using MKL

function assemble_local(fe::FESpace, mesh::Mesh)
  Aloc = Float64.(stiff_3d_sp(fe.w_128, fe.D_128, mesh.jd[1], mesh.jd[2], mesh.jd[3]))
  Mloc = Float64.(reshape(map(prod, Base.product(ntuple(x -> fe.w_128, mesh.dim)...)), length(fe.w)^mesh.dim))
  Mloc *= mesh.jd[1] * mesh.jd[2] * mesh.jd[3]
  return Aloc, Mloc
end

function cg(Avmult, b, x, tol)
  r = b - Avmult(x)
  p = r
  rsold = r' * r
  if sqrt(rsold) < tol
    println("No iterations needed")
    return 0, x
  end

  i = 1
  while i <= length(b)
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
    i += 1
  end
  println("Residual:   ", sqrt(rsold))
  println("Iterations: ", i)
  return i, x
end

function compute_errors(u, un, dun, fe, mesh, dof_map)
  dim = Int64(round(log(size(dof_map, 1)) / log(length(fe.w)), RoundNearest))
  ww = map(prod, Base.product(ntuple(x -> fe.w, dim)...))
  ww = reshape(ww, length(ww))
  l2_err = 0
  h1_err = 0
  for li in eachcol(dof_map)
    j = mesh.jd[1] * mesh.jd[2] * mesh.jd[3]
    l2_err += j * sum((u[li] - un[li]) .^ 2 .* ww)

    u_loc = reshape(un[li], ntuple(x -> length(fe.w), dim))
    ux = Array{Float64}(undef, size(u_loc))
    uy = Array{Float64}(undef, size(u_loc))
    uz = Array{Float64}(undef, size(u_loc))
    for i = 1:length(fe.w)
      for j = 1:length(fe.w)
        ux[:, i, j] = fe.D * u_loc[:, i, j] / mesh.jd[1]
        uy[i, :, j] = fe.D * u_loc[i, :, j] / mesh.jd[2]
        uz[i, j, :] = fe.D * u_loc[i, j, :] / mesh.jd[3]
      end
    end
    ux = reshape(ux, length(ux))
    uy = reshape(uy, length(uy))
    uz = reshape(uz, length(uz))
    h1_err += j * sum(((ux - dun[1][li]) .^ 2 + (uy - dun[2][li]) .^ 2 + (uz - dun[3][li]) .^ 2) .* ww)  
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

function dsol(xyz)
  pi2 = 2 * π
  sinxyz = sin.(2 * π * xyz)
  cosxyz = cos.(2 * π * xyz)
  return pi2*cosxyz[1, :].*cosxyz[2, :].*sinxyz[3, :], 
        -pi2*sinxyz[1, :].*sinxyz[2, :].*sinxyz[3, :], 
         pi2*sinxyz[1, :].*cosxyz[2, :].*cosxyz[3, :]
end

function vmult(Aloc, x, dof_map, bc)
  y = zeros(Float64, size(x))
  for i in axes(dof_map, 2)
    for j in axes(dof_map, 1)
      @inbounds jj = dof_map[j, i]
      @inbounds @simd for k = Aloc.row_ptr[j] + 1 : Aloc.row_ptr[j+1]
        @inbounds y[jj] += Aloc.val[k] * x[dof_map[Aloc.col_ind[k], i]]
      end
    end 
  end
  @simd for i in eachindex(bc)
    @inbounds y[bc[i]] = x[bc[i]]
  end
  return y
end

degs = [1:8;]
nels = 2 .^ [0:8;]
err_table = zeros(Float64, length(nels), length(degs), 4)
elapsed_times = zeros(Float64, 3, 0)

for i in eachindex(nels)
  for j in eachindex(degs)
    nel = nels[i]
    deg = degs[j]
    println("Basis deg: ", deg)
    println("#elements: ", nel, " (per dimension)")
    ndof = (nel * deg + 1)^3
    println("Problem has ", ndof, " dofs")

    if ndof > 1000000
      println("Skip")
      println("-------------------------------------")
      continue
    end

    println("Building FE space")
    fespace = FESpace(deg)
    println("Building mesh")
    mesh = CartesianMesh([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], nel*[1, 1, 1])
    println("Distributing dof")
    dof_map, dof_support = distribute_dof_v2(mesh, deg)
    
    println("Assembling system")
    Aloc, Mloc = assemble_local(fespace, mesh)
    Adloc = diag(Aloc)
    M = zeros(Float64, ndof)
    Ad = zeros(Float64, ndof)
    for idx in eachcol(dof_map)
      M[idx] += Mloc
      Ad[idx] += Adloc
    end
    f = force(dof_support) .* M
    println("Evaluating exact solution")
    u_ex = sol(dof_support)
    du_ex = zeros(Float64, size(u_ex, 1), 3)
    println("Applying BC")
    x0 = zeros(Float64, size(f))
    bc = find_bc(mesh, dof_support)
    f[bc] = u_ex[bc]
    x0[bc] = u_ex[bc]

    println("Solving system with CG")
    #Aloc = sparse(Aloc) # WARNING: it is CSC, with is suboptimal for vmult
    Aloc = CSR(Aloc)
    t = @elapsed it, u = cg(x->vmult(Aloc, x, dof_map, bc), f, x0, 1e-8)
    println("Computing error")
    println("Res: ",  norm(vmult(Aloc, u, dof_map, bc) - f) / norm(f))
    println("Err: ",  norm(u - u_ex))
    du_ex = dsol(dof_support)
    err_l2, err_h1 = compute_errors(u, u_ex, du_ex, fespace, mesh, dof_map)
    println("Err L2:        ", err_l2)
    println("Err H1 (semi): ", err_h1)

    err_table[i, j, 1] = norm(u - u_ex)
    err_table[i, j, 2] = norm(u - u_ex, Inf)
    err_table[i, j, 3] = err_l2
    err_table[i, j, 4] = err_h1

    global elapsed_times
    elapsed_times = hcat(elapsed_times, [t, deg, ndof])
    display(elapsed_times)
    println("-------------------------------------")
  end
end

#display(err_table)

conv_table = log.(err_table[2:end, :, :]./err_table[1:end-1, :, :])./log.(nels[1:end-1]./nels[2:end]) 
display(conv_table[:, :, 2:4])
