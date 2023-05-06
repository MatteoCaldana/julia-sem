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
  Avmult = x->vmult_generic(x, grid.dof_map, grid.boundary_idxs)
  r = Avmult(x) - b;
  r0 = norm(r);
  println("Initial residual is: ", r0)
  for i=1:num_vcyc
      x = vcycle(grid, vmult_generic, get_diag, v1, v2, b, x);
      r = Avmult(x) - b;
      println("|res| = ", norm(r))
      if (norm(r)/r0 < tol)
          iter = i;
          rr = norm(r)/r0;
          return;
      end
  end
  println("------------------------------------------");
  iter = num_vcyc;
  rr = norm(r)/r0;
  return iter, x
end

function vcycle(grid, vmult_generic, get_diag, v1, v2, b, x)::Vector{Float64}
  Avmult = x->vmult_generic(x, grid.dof_map, grid.boundary_idxs)
  if ( !grid.has_coarse )
    #TODO: if len(bc_values) == len(x) return bc_values
    return cg(Avmult, b, x, 1e-15)[2] #discard number of iterations
  end
  # TODO: cache invD_omega
  invD_omega = 1 ./ get_diag(grid.dof_map)
  for it = 1:v1
    x += invD_omega .* (b - Avmult(x));
  end
  res = Avmult(x) - b;
  res_coarse = projection_vmult(grid.coarse, res);
  res_coarse[grid.coarse.boundary_idxs] = grid.coarse.boundary_values;
  x_corr_coarse = vcycle(grid.coarse, vmult_generic, get_diag, v1, v2, res_coarse, zeros(size(res_coarse)));
  x -= interpolation_vmult(grid.coarse, x_corr_coarse);
  for it = 1:v2
    x += invD_omega .* (b - Avmult(x));
  end
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