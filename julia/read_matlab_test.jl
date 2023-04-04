using MAT 
using LinearAlgebra

nel = 4
deg = 3
vars = matread("test-data/3d_stiff_deg"*repr(deg)*"_nel"*repr(nel)*".mat")

A = vars["A"]
f = vars["f"]
u = vars["un"]

u_julia = A \ f
e = u - u_julia

println(norm(e, 2))
println(norm(e, Inf))