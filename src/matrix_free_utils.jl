
function materialize_linear_map(linear_map, n::Int64, m::Int64)::Matrix{Float64}
  material_map = Array{Float64}(undef, n, m)
  e = zeros(m, 1)
  for j=1:m
    e[j] = 1.0
    material_map[:, j] = linear_map(e)
    e[j] = 0.0
  end
  return material_map
end