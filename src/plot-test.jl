using Plots

function spy(A)
  spyx = Vector{Int64}(undef, 0)
  spyy = Vector{Int64}(undef, 0)

  for i = 1:size(A, 1)
    for j = 1:size(A, 2)
      if abs(A[i, j]) > 1e-14
        push!(spyx, i)
        push!(spyy, j)
      end
    end
  end

  scatter(spyx, spyy)
  gui()
  read(stdin, Char)
end

scatter([0, 0], [1, 2])
gui()
read(stdin, Char)