using BenchmarkTools

function print_benckmark(b)
  io = IOBuffer()
  show(io, "text/plain", b)
  s = String(take!(io))
  println(s)
end

a = 2 .^ (rand(UInt64, 1000) .% 4)

b = @benchmark alog = 64 .- leading_zeros.(a)
print_benckmark(b)

b = @benchmark alog = Int64.(log2.(Float64.(a)))
print_benckmark(b)


