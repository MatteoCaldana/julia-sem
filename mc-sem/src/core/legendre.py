from sympy import *
from sympy.abc import x
from functools import reduce

def gcd(a, b):
    while b:      
        a, b = b, a % b
    return a

def lcm(a, b):
    return a * b // gcd(a, b)

def max_2pow(n):
  res = 0
  while 2**(res+1) < n:
    res += 1
  return res

def coeff_to_expr(coeff):
  coeff.reverse()
  pow = max_2pow(len(coeff))

  res = ""
  for i in range(pow):
    res += f"  const auto x{i+1} = x{i} * x{i};\n"
  
  den = reduce(lcm, [c.q for c in coeff])
  res += f"  return ("
  if coeff[0].p != 0:
    res += f"{coeff[0].p * den // (coeff[0].q)} "

  for i in range(1, len(coeff)):
    if coeff[i].p != 0: 
      x = ""
      for j in range(10): #max 1024
        if int(i) & int(2**j):
          x += f"x{j}*"
      x = x[:-1]
      if res[-1] != "(":
        res += "+"
      res += f" ({coeff[i].p * den // (coeff[i].q)}.) * {x} "
  res += f") / {den}.;\n"
  return res

def to_fun(name, i, content):
  fun = f"double {name}_{i:02d}(double x0) " + "{\n"
  fun += content
  fun += "}"
  return fun

for i in range(20):
  leg = Poly(legendre(i, x), x)
  leg_prime = leg.diff()
  leg_second = leg_prime.diff()
  print(to_fun('leg__', i, coeff_to_expr(leg.all_coeffs())))
  print('')
  print(to_fun('leg_p', i, coeff_to_expr(leg_prime.all_coeffs())))
  print('')
  print(to_fun('legpp', i, coeff_to_expr(leg_second.all_coeffs())))
  print('')
