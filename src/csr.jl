
struct CSR
  val::Vector{Float64}
  col_ind::Vector{Int64}
  row_ptr::Vector{Int64}

  function CSR(m::Matrix{Float64})
    val = zeros(Float64, 0)
    col_ind = zeros(Int64, 0)
    row_ptr = zeros(Int64, 1)
    for i in axes(m, 1)
      for j in axes(m, 2)
        if abs(m[i, j]) > eps(Float64)
          push!(val, m[i, j])
          push!(col_ind, j)            
        end
      end
      push!(row_ptr, length(val))
    end
    return new(val, col_ind, row_ptr)
  end
end