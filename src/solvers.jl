function cg(Avmult, b, x, tol)
  r = b - Avmult(x)
  p = r
  rsold = dot(r, r)
  r0 = sqrt(rsold)
  
  if r0 < tol
    println("No iterations needed")
    return 0, x
  end

  i = 1
  while i <= length(b)
    Ap = Avmult(p)
    alpha = rsold / (p' * Ap)
    x += alpha * p
    r -= alpha * Ap
    rsnew = dot(r, r)
    if sqrt(rsnew) / r0 < tol
      rsold = rsnew
      break
    end
    p = r + (rsnew / rsold) * p
    rsold = rsnew
    i += 1
  end
  println("Residual:   ", sqrt(rsold) / r0, "(r0=", r0, ")")
  println("Iterations: ", i)
  return i, x
end

function pcg(Avmult, b, Pvmult, x, tol)
  r = b - Avmult(x)
  z = Pvmult(r)
  p = z
  rzold = r' * z
  r0 = norm(r)

  if r0 < tol
    println("No iterations needed")
    return 0, x
  end

  rn = norm(r)
  i = 1
  while i <= length(b)
    Ap = Avmult(p)
    alpha = rzold / (p' * Ap)
    x += alpha * p
    r -= alpha * Ap
    rznew = r' * z
    if rn / r0 < tol
      break
    end
    rn = norm(r)
    z = Pvmult(r)
    p = z + (rznew / rzold) * p
    rzold = rznew
    i += 1
  end
  println("Residual:   ", rn / r0, "(r0=", r0, ")")
  println("Iterations: ", i)
  return i, x
end

function jacobi(Avmult, 
  b::Vector{Float64}, 
  x::Vector{Float64}, 
  nit::Int64, 
  tol::Float64,
  invD_omega::Vector{Float64})

  r = b - Avmult(x)
  r0 = norm(r)

  if r0 < tol
    println("No iterations needed")
    return 0, x
  end

  i = 1
  while i <= nit
    r = b - Avmult(x);
    if norm(r) / r0 < tol
      break
    end
    x += invD_omega .* r;
    i += 1
  end
  println("Residual:   ", norm(r) / r0, "(r0=", r0, ")")
  println("Iterations: ", i)
  return i, x
end

function solve_mg(grid, vmult_generic, get_diag, num_vcyc, v1, v2, b, x, tol)
  Avmult = x->vmult_generic(grid.A, x, grid.dof_map, grid.boundary_idxs)
  r = Avmult(x) - b;
  r0 = norm(r);
  println("Initial residual is: ", r0)
  for i=1:num_vcyc
      println("Vcycle iter: ", i)
      x = vcycle(grid, vmult_generic, get_diag, v1, v2, b, x, i);
      r = Avmult(x) - b;
      println("|res| = ", norm(r))
      println("======================================================")
      println("======================================================")
      if (norm(r)/r0 < tol)
          iter = i;
          rr = norm(r)/r0;
          return iter, x;
      end
  end
  iter = num_vcyc;
  rr = norm(r)/r0;
  return iter, x
end

using MAT
include("matrix_free_utils.jl")

function fix_subnormal(x)
  return ifelse.(abs.(x) .< 2.0 * eps(Float64), 0.0, x)
end

function vcycle(grid, vmult, get_diag, v1, v2, b, x, i)::Vector{Float64}
  Avmult = x->vmult(grid.A, x, grid.dof_map, grid.boundary_idxs)

  if i == 1
    if grid.has_coarse
      Pmat = matread("../test-data/P."*repr(grid.ndofs)*".mat")["PP"]
      P = Matrix(grid.coarse.P)
      println("Error P:", maximum(abs.(P - Pmat)))

      Rmat = matread("../test-data/R."*repr(grid.coarse.ndofs)*".mat")["RR"]
      R = Matrix(grid.coarse.P')
      println("Error R:", maximum(abs.(R - Rmat)))
    end

    Amat = matread("../test-data/K."*repr(grid.ndofs)*".mat")["KK"]
    A = materialize_linear_map(Avmult, grid.ndofs, grid.ndofs)
    println("Error A:", maximum(abs.(A - Amat)))
    display(diag(A))
    display(diag(Amat))
  end

  bmat = matread("../test-data/rhs."*repr(grid.ndofs)*".it"*repr(i)*".mat")["rhs"]
  println("Error rhs:", maximum(abs.(b - bmat)))

  xmat = matread("../test-data/u."*repr(grid.ndofs)*".0.it"*repr(i)*".mat")["u"]
  println("Error x (0):", maximum(abs.(x - xmat)))


  # 0. handle coarses level
  if ( !grid.has_coarse )
    #TODO: if len(bc_values) == len(x) return bc_values
    return cg(Avmult, b, x, 1e-15)[2] #discard number of iterations
  end

  # TODO: save once for all invD_omega
  invD_omega = 1 ./ get_diag(grid.Ad, grid.dof_map, grid.boundary_idxs)
  if i == 1

    JiDmat = matread("../test-data/JiD."*repr(grid.ndofs)*".mat")["JiD"]
    display(grid.A)
    display(grid.Ad)
    println("Error JiD:", maximum(abs.(invD_omega - JiDmat)))
    display(invD_omega)
    display(JiDmat)
  end

  # 1. pre-smooth
  for it = 1:v1
    x += invD_omega .* (b - Avmult(x));
  end

  xmat = matread("../test-data/u."*repr(grid.ndofs)*".1.it"*repr(i)*".mat")["u"]
  println("Error x (1):", maximum(abs.(x - xmat)))

  # 2. compute residual
  res = Avmult(x) - b;

  resmat = matread("../test-data/res."*repr(grid.ndofs)*".it"*repr(i)*".mat")["res"]
  println("Error res:", maximum(abs.(res - resmat)))

  # 3. restrict
  res_coarse = grid.coarse.P' * res;
  res_coarsemat = matread("../test-data/res_coarse."*repr(grid.coarse.ndofs)*".raw.it"*repr(i)*".mat")["res_coarse"]
  println("Error res_coarse raw:", maximum(abs.(res_coarse - res_coarsemat)))

  res_coarse[grid.coarse.boundary_idxs] = grid.coarse.boundary_values;

  bdy_idxmat = matread("../test-data/bdy_idx."*repr(size(grid.coarse.boundary_idxs, 1))*".it"*repr(i)*".mat")["bdy_idx"]
  println("Error bdy_idx:", maximum(abs.(grid.coarse.boundary_idxs - bdy_idxmat)))

  res_coarsemat = matread("../test-data/res_coarse."*repr(grid.coarse.ndofs)*".bdy.it"*repr(i)*".mat")["res_coarse"]
  println("Error res_coarse bdy:", maximum(abs.(res_coarse - res_coarsemat)))

  # 4. recurse
  x_corr_coarse = vcycle(grid.coarse, vmult, get_diag, v1, v2, res_coarse, zeros(size(res_coarse)), i);

  u_corr_coarsemat = matread("../test-data/u_corr_coarse."*repr(grid.coarse.ndofs)*".it"*repr(i)*".mat")["u_corr_coarse"]
  println("Error u_corr_coarse:", maximum(abs.(x_corr_coarse - u_corr_coarsemat)))

  # 5. prolong and correct
  x -= grid.coarse.P * x_corr_coarse;

  xmat = matread("../test-data/u."*repr(grid.ndofs)*".5.it"*repr(i)*".mat")["u"]
  println("Error u (5):", maximum(abs.(x - xmat)))

  # 6. post-smooth
  for it = 1:v2
    r = invD_omega .* (Avmult(x) - b);
    x -= r;
  end

  xmat = matread("../test-data/u."*repr(grid.ndofs)*".6.it"*repr(i)*".mat")["u"]
  println("Error u (6):", maximum(abs.(x - xmat)))
  return x
end

function chebyshev(Avmult, b, x, nit, eig_max, eig_min)
  # adjust the eigenvalues to hit the upper spectrum
  # beta = eig_max;
  # alpha = 0.25*eig_max;% (grid.eig_min + grid.eig_max)/2;

  # delta = (beta - alpha)/2;
  # theta = (beta + alpha)/2;
  # s1 = theta/delta;
  # rhok = 1/s1;

  # d = zeros(size(u));

  # res = -grid.residual ( rhs, u );
  # d = res/theta.* grid.jacobi_invdiag;
  # u = u + d;

  # for iter = 2:v
  #     rhokp1 = 1/ (2*s1 - rhok);
  #     d1 = rhokp1 * rhok;
  #     d2 = 2*rhokp1 / delta;
  #     rhok = rhokp1;
  #     res = -grid.residual ( rhs, u );
  #     d = d1 * d + d2 * res.*grid.jacobi_invdiag;
  #     u = u + d;
  # end
  return
end

function power()
end