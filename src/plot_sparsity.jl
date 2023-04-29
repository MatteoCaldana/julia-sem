using Plots

include("mesh.jl")

function spy(nel, deg)
  mesh = CartesianMesh([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], nel*[1, 1, 1])
  dof_map, _ = distribute_dof_v2(mesh, deg)

  
  default(show = true)
  for i in axes(dof_map, 2)
    println(i, " / ", size(dof_map, 2))
    dof = dof_map[:, i]
    sparsity = vec(collect(Base.Iterators.product(Base.Iterators.repeated(dof, 2)...)))
    r = rand(1:length(sparsity), Int64(round(0.02 * length(sparsity))))
    spyx = [sparsity[i][1] for i in r]
    spyy = [sparsity[i][2] for i in r]
    scatter!(spyx, spyy, mc=:red, ms=2, ma=0.5, label="")
  end

  gui()
  read(stdin, Char)
end

spy(4, 4)

