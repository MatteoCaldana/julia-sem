using DoubleFloats  # Double64
using Quadmath      # Float128
using BenchmarkTools
using SparseArrays

function pnleg(x::Array{Float128,1}, n::Int64)::Array{Float128,1}
  if n == 0
    return ones(Float128, size(x))
  else
    p1 = ones(Float128, size(x))
    p2, p3 = x, x
    for k = 1:n-1
      p3 = (Float128(2 * k + 1) * x .* p2 - k * p1) / Float128(k + 1)
      p1, p2 = p2, p3
    end
    return p3
  end
end

function pnleg1(x::Array{Float128,1}, n::Int64)::Tuple{Vector{Float128},Vector{Float128}}
  if n == 0
    return zeros(Float128, size(x)), ones(Float128, size(x))
  elseif n == 1
    return ones(Float128, size(x)), x
  else
    p1 = x
    p2 = x
    p11 = x
    p21 = x
    p3 = x
    p31 = x
    p1 = ones(Float128, size(x))
    p2 = x
    p11 = zeros(Float128, size(x))
    p21 = ones(Float128, size(x))
    for k = 1:n-1
      k2p1 = Float128(2 * k + 1)
      kp1 = Float128(1) / Float128(k + 1)
      p3 = (k2p1 * x .* p2 - k * p1) * kp1
      p31 = (k2p1 * (x .* p21 + p2) - k * p11) * kp1
      p11, p21 = p21, p31
      p1, p2 = p2, p3
    end
    return p31, p3
  end
end

function pnleg2(x::Array{Float128,1}, n::Int64)::Tuple{Vector{Float128},Vector{Float128},Vector{Float128}}
  if n <= 1
    return zeros(Float128, size(x)), zeros(Float128, size(x)), zeros(Float128, size(x))
  elseif n == 2
    return zeros(Float128, size(x)), zeros(Float128, size(x)), ones(Float128, size(x))
  else
    p1 = ones(Float128, size(x))
    p2 = x
    p11 = zeros(Float128, size(x))
    p21 = ones(Float128, size(x))
    p12 = zeros(Float128, size(x))
    p22 = zeros(Float128, size(x))
    p3 = x
    p31 = x
    p32 = x
    for k = 1:n-1
      k2p1 = Float128(2 * k + 1)
      kp1 = Float128(1) / Float128(k + 1)
      p3 = (k2p1 * x .* p2 - k * p1) * kp1
      p31 = (k2p1 * (x .* p21 + p2) - k * p11) * kp1
      p32 = (k2p1 * (x .* p22 + p21 * 2) - k * p12) * kp1
      p1, p2 = p2, p3
      p11, p21 = p21, p31
      p12, p22 = p22, p32
    end
    return p32, p31, p3
  end
end


function ixlgl(n::Int64)::Vector{Float128}
  x = ones(Float128, n - 2)
  for i = 1:n-2
    x[i] = cos(Float128(pi) * Float128(i) / Float128(n - 1))
  end
  for i = 1:40
    p2, p1, p = pnleg2(x, n - 1)
    dx = p1 ./ p2
    dx = ifelse.(isnan.(dx), Float128(0), dx)
    x -= dx
    if maximum(abs.(dx)) < Float128(10.0 * eps(Float128))
      break
    end
  end
  return sort(x)
end

function jacobi_eval(x::Array{Float128,1}, n::Int64, alpha::Float128, beta::Float128)::Tuple{Array{Float128,1},Array{Float128,1}}
  apb = alpha + beta
  ab2 = alpha * alpha - beta * beta

  if n == 0
    return ones(Float128, size(x)), zeros(Float128, size(x))
  else
    p1 = ones(Float128, size(x))
    p2 = zeros(Float128, size(x))
    pd1 = zeros(Float128, size(x))
    pd2 = zeros(Float128, size(x))
    p = (alpha - beta .+ (apb + 2) * x) * 0.5
    pd = 0.5 * (apb + 2) * ones(Float128, size(x))
    for k = 1:n-1
      k2 = k * 2
      k2ab = k2 + alpha + beta
      k1 = k + 1
      k2 = k * 2
      k2ab = k2 + alpha + beta
      k2ab1 = k2ab + 1
      k2ab2 = k2ab1 + 1
      p2, p1 = p1, p
      pd2, pd1 = pd1, pd
      a1 = 2 * k1 * (k1 + apb) * k2ab
      a21 = k2ab1 * ab2
      a22 = k2ab2 * k2ab1 * k2ab
      a3 = 2 * (k + alpha) * (k + beta) * k2ab2
      p = ((a21 .+ a22 * x) .* p1 - a3 * p2) / a1
      pd = (a22 * p1 + (a21 .+ a22 * x) .* pd1 - a3 * pd2) / a1
    end
    return p, pd
  end
end

function jacobi_roots(n::Int64, alpha::Float128, beta::Float128)::Array{Float128,1}
  TOL = Float128(10.0 * eps(Float128))
  KMAX = 40
  @assert n > 0 "polinomial degree is not positive"

  x = zeros(Float128, n)
  x0 = [cos(Float128(pi) / Float128((2 * n)))]
  x1 = x0
  for j = 1:n
    diff = TOL + 1
    kiter = 0
    while kiter <= KMAX && diff >= TOL
      p, pd = jacobi_eval(x0, n, alpha, beta)
      ss = j > 1 ? sum(Float128(1) ./ (x0[1] .- x[1:j-1])) : 0
      x1 = x0 .- p[1] / (pd[1] - ss * p[1])
      diff = abs(x1[1] - x0[1])
      kiter = kiter + 1
      x0 = x1
    end
    x0 = (x1 .+ cos((2 * (j + 1) - 1) * pi / (2 * n))) / 2.0
    x[j] = x1[1]
  end
  return sort(x)
end


function xwlgl(np::Int64)::Tuple{Array{Float128,1},Array{Float128,1}}
  @assert np > 1 "the number of quadrature nodes is smaller than 2"

  if np == 2
    return [-1, 1], [1, 1]
  else
    x = zeros(Float128, np)
    w = zeros(Float128, np)
    n = np - 1
    x[1] = -1
    x[np] = 1
    x[2:n] = jacobi_roots(n - 1, Float128(1), Float128(1))
    coef = Float128(2) / Float128(n * np)
    w = coef ./ (pnleg(x, n) .^ 2)
  end
  return x, w
end


function derlgl(x::Array{Float128,1})::Array{Float128,2}
  np = size(x)[1]
  d = zeros(Float128, np, np)
  n = np - 1
  for j = 1:np
    lnxj = pnleg([x[j]], n)[1]
    for i = 1:np
      if i != j
        lnxi = pnleg([x[i]], n)[1]
        d[i, j] = lnxi / ((x[i] - x[j]) * lnxj)
      end
    end
  end
  d[1, 1] = -0.25 * n * np
  d[np, np] = -d[1, 1]
  return d
end


function vandermonde(x::Vector{Float128})::Tuple{Matrix{Float128},Matrix{Float128}}
  np = size(x)[1]
  V = ones(Float128, np, np)
  for i = 1:np-1
    V[:, i+1] = V[:, i] .* x
  end
  return V, inv(V)
end

function derlgl_v2(Vx::Matrix{Float128}, Vinv::Matrix{Float128})::Matrix{Float128}
  np = size(Vinv)[1]
  return Vx[:, 1:end-1] * (Vinv[2:end, :] .* range(1, np - 1))
end


function test()
  for n = 3:16
    @time x, w = xwlgl(n)
    @time x2 = ixlgl(n)

    @time d = derlgl(x)
    @time V, Vinv = vandermonde(x)
    @time d2 = derlgl_v2(V, Vinv)

    z = pnleg1(x, n - 1)[1][2:end-1]
    println("Jacobi deflation residual:          ", maximum(abs.(z)))
    println("Jacobi deflation vs Newton-Raphson: ", maximum(abs.(x[2:end-1] - x2)))


    println("Exact D vs Vandermonde:             ", maximum(abs.(d2 - d)))
    println("=============================================================================")
  end
end

function stiff_2d_sp(w, d, jx, jy)
  np = length(w)
  nn = np * np
  jyx = jy / jx
  jxy = jx / jy
  ww = w .* w'
  ww = ww[:]
  A = zeros(nn, nn)
  for i = 1:np
    inde = (i-1)*np+1:i*np
    A[inde, inde] += d' * (ww[inde] .* d) * jyx
    inde = i:np:nn
    A[inde, inde] += d' * (ww[inde] .* d) * jxy
  end
  return A
end

function stiff_3d_sp(w, d, jx, jy, jz)
  np = length(w)
  mn = np * np
  nn = np * np * np
  j = jx * jy * jz
  jyx = jy / jx
  jxy = jx / jy
  s = d' * (w .* d)
  A = zeros(nn, nn)

  for li = 1:np
    il = (li - 1) * mn
    for ki = 1:np
      ik = (ki - 1) * np
      coefx = w[ki] * w[li] * jz * jy / jx
      for i = 1:np
        ii = il + ik + i
        coefy = w[i] * w[li] * jz * jx / jy
        coefz = w[i] * w[ki] * jy * jx / jz

        i1 = il + ik + 1
        i2 = il + ik + np
        A[ii, i1:i2] = A[ii, i1:i2] + s[i, 1:np] * coefx

        i1 = il + i
        i2 = il + mn
        A[ii, i1:np:i2] = A[ii, i1:np:i2] + s[ki, 1:np] * coefy

        i1 = ik + i
        A[ii, i1:mn:nn] = A[ii, i1:mn:nn] + s[li, 1:np] * coefz
      end
    end
  end

  return A
end

function to_float64_fix_subnormal(x)
  return ifelse.(abs.(x) .< Float128(2.0 * eps(Float64)), 0.0, Float64.(x))
end

struct FESpace
  n::Int64
  np::Int64

  x::Vector{Float64}        # 1D reference coordinates of the interpolation nodes
  w::Vector{Float64}        # 1D weights for gll quadrature
  D::Matrix{Float64}        # 1D derivative of Lagrange interpolants at the interpolation nodes
                            # the same of (derlgl, derlgl_v2, transpose(V \ gradV))
  V::Matrix{Float64}        # 1D Vandermonde matrix of Legendre polynomials
  gradV::Matrix{Float64}    #    and their derivative (Nrp x Nrp)
  Ph::Matrix{Float64}       # interpolation from this element to its 4/8 children
  Pp::Matrix{Float64}       # interpolation from this element to its 2p version
  M::Matrix{Float64}        # exact 1D Mass matrix (Nrp x Nrp) at gll

  x_128::Vector{Float128}
  w_128::Vector{Float128}
  D_128::Matrix{Float128}
  V_128::Matrix{Float128}
  gradV_128::Matrix{Float128}
  Ph_128::Matrix{Float128}
  Pp_128::Matrix{Float128}
  M_128::Matrix{Float128}

  function FESpace(n::Int64)
    np = n + 1
    x, w = xwlgl(np)
    d = derlgl(x)

    x_hby2  = [Float128(0.5)*(x .- Float128(1))...; Float128(0.5)*(x[2:end] .+ Float128(1))...]
    x_2p, _ = xwlgl(2*n + 1)

    V     = zeros(Float128, np, np)
    gradV = zeros(Float128, np, np)
    Vph   = zeros(Float128, np, 2*n+1)
    Vpp   = zeros(Float128, np, 2*n+1)
    for i=1:np
      p, pd = jacobi_eval(x, i-1, Float128(0.0), Float128(0.0))
      V[i,:]     = p
      gradV[i,:] = pd
      Vph[i,:]   = jacobi_eval(x_hby2, i-1, Float128(0.0), Float128(0.0))[1]
      Vpp[i,:]   = jacobi_eval(x_2p,   i-1, Float128(0.0), Float128(0.0))[1]
    end

    #D = transpose(V \ gradV);

    p_h_1d = transpose(V \ Vph);
    p_p_1d = transpose(V \ Vpp);

    iV    = inv(V);
    M     = iV * iV';
    invM  = inv(M);


    return new(n, np, 
      to_float64_fix_subnormal(x), 
      to_float64_fix_subnormal(w), 
      to_float64_fix_subnormal(d),
      to_float64_fix_subnormal(V),
      to_float64_fix_subnormal(gradV),
      to_float64_fix_subnormal(p_h_1d),
      to_float64_fix_subnormal(p_p_1d),
      to_float64_fix_subnormal(M),
      x, 
      w, 
      d,
      V,
      gradV,
      p_h_1d,
      p_p_1d,
      M,    
    )
  end
end

function print_benckmark(b)
  io = IOBuffer()
  show(io, "text/plain", b)
  s = String(take!(io))
  println(s)
end


function bench()
  for np = 2:9
    println("========================================================")
    println("n: ", np - 1)
    println("--------------------------------------------------------")
    x, w = xwlgl(np)
    d = derlgl(x)

    A = Float64.(stiff_2d_sp(w, d, 1.0, 1.0))
    u = rand(size(A, 1))

    b = @benchmark v = 2 * d
    print_benckmark(b)
    println("========================================================")
  end

end

# test()

# small = Float128(1e-33)
# println(small, Float64(small))

function bench(n)
  np = n + 1
  println("========================================================")
  println("n: ", n)
  println("--------------------------------------------------------")
  println("Dense MM")
  x, w = xwlgl(np)
  d = derlgl(x)

  A = Float64.(stiff_3d_sp(w, d, 1.0, 1.0, 1.0))
  u = rand(size(A, 1))

  b = @benchmark v = $A * $u
  print_benckmark(b)
  println("--------------------------------------------------------")
  println("Kron trick")

  d64 = Float64.(derlgl(x))
  b = @benchmark v = $d64 * reshape($u, $np, $np * $np)
  print_benckmark(b)

  println("--------------------------------------------------------")
  println("Sparse CSC")
  Asp = sparse(A)
  println(nnz(Asp), " ", nnz(Asp) / length(A), " ", sum(abs.(A) .> 0.0), " ", sum(abs.(A) .> 1e-12))
  b = @benchmark v = $Asp * $u
  print_benckmark(b)

  println("--------------------------------------------------------")
  println("Sparse CSC transpose")
  b = @benchmark v = transpose($Asp) * $u
  print_benckmark(b)

  println("--------------------------------------------------------")
  println("Sparse CSR")
  a = zeros(Float64, 0)
  c = zeros(Int64, 0)
  r = ones(Int64, size(A, 1) + 1)
  for i = 1:size(A, 1)
    for j = 1:size(A, 2)
      if abs(A[i, j]) > 2 * eps(Float64)
        push!(a, A[i, j])
        push!(c, j)
      end
    end
    r[i+1] = length(a)
  end

  function mul(v::Vector{Float64}, c::Vector{Int64}, r::Vector{Int64}, u::Vector{Float64})::Vector{Float64}
    y = zeros(Float64, size(u))
    for row = 1:length(u)
      @inbounds for j = r[row]:r[row+1]
        @inbounds y[row] += v[j] * u[c[j]]
      end
    end
    return y
  end

  b = @benchmark v = $mul($a, $c, $r, $u)
  print_benckmark(b)

  println("========================================================")
end



  # for n = [1, 2, 4, 8]
  #   bench(n)
  # end