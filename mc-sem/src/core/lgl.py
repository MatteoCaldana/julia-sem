import numpy as np

from scipy import special
import sympy
from sympy.abc import x

def lglnodes(N):
  N1 = N+1
  # Use the Chebyshev-Gauss-Lobatto nodes as the first guess
  x = np.cos(np.pi*(np.arange(0,N+1))/N)
  # The Legendre Vandermonde Matrix
  P = np.zeros((N1, N1))
  # Compute P_(N) using the recursion relation
  # Compute its first and second derivatives and 
  # update x using the Newton-Raphson method.
  xold = 2
  for i in range(1000): 
    if max(abs(x-xold)) < 1e-17:
      break
    xold = x
    P[:, 0] = 1    
    P[:, 1] = x
    for k in range(1, N):
      P[:, k+1] = ( (2*(k+1)-1)*x*P[:,k] - k*P[:,k-1] ) / (k+1)
    x = xold - ( x*P[:,-1] - P[:,-2] ) / ( N1*P[:,-1] )
  return x

def sym_lglnodes(n):
  legp = sympy.Poly(sympy.legendre(n, x), x).diff()
  legpp = legp.diff()
  np_legp = np.vectorize(lambda x: legp.eval(x))
  np_legpp = np.vectorize(lambda x: legpp.eval(x))
  # Use the Chebyshev-Gauss-Lobatto nodes as the first guess
  x0 = np.cos(np.pi*(np.arange(0, n+1))/n)[1:-1]
  for _ in range(100):
    x0 = x0 - np_legp(x0) / np_legpp(x0)
  return x0


def test(n, xx):
  poly = special.legendre(n)  # coefficients of n^th degree Legendre polynomial 
  polyd = poly.deriv() # coefficients of derivative of n^th degree Legendre Polynomial
  evald = np.polyval(polyd, xx) # evaluate derivative at desired coordinates(s)

  legp = sympy.Poly(sympy.legendre(n, x), x).diff()
  np_legp = np.vectorize(lambda x: legp.eval(x))
  sym_evald = np_legp(xx)
  norm = np.abs(sym_evald).max()
  return np.abs(evald).max(), norm.p / norm.q

if __name__ == "__main__":
  for i in range(3, 20):
    n1 = lglnodes(i)[1:-1]
    n2 = sym_lglnodes(i)
    print(max(abs(n1 - n2)))
    print('newt', test(i, n1))
    print('sym ', test(i, n2))
    print('/////')
    # array = ', '.join([f'{float(x):.17e}' for x in n2])
    # print(f'{{1., {array}, -1.}},')