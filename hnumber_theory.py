# Simple number theory in Python

# Ferma(n, L)			Fermat primality test of n by bases in list L
# Miller_Rabin(n, L)	M-R primality test of n by bases in list L
# Frobenius_index(n)    Frobenius index
# is_perfect_square(n)  test for perfect square
# Frobenius(n)			Frobenius primality test
# is_prime(n):          return Frobenius(n)
# Jacobi(a,b)			Jacobi symbol
# zPower(a,b,c,k0,n)    (a + b*sqrt(c))^k0 mod n, return (a1,b1)
# nextPrime(n)			prime, greater then n
# Pollard(n)     		Pollard's Rho Algorithm for Prime Factorization
# Factor(n)				return the list of prime dividers
# egcd(a,b)				return 3 numbers (d,u,v), such d=GCD(a,a) and d = u*a + v*a
# chRem2(a1,n1, a2,n2)  Chinese remainder theorem
#-iSqrt(x)              int(sqrt(x)) change to math.isqrt(x)
#-mod_inv(a,n)			1/a mod n   change for pow(a,-1,n)
#-mod_div(a,b,n)		a/b mod n   change for (a * pow(b,-1,n))%n


# Example:
# hnumer_theory.Ferma(341,[2,3])

# For impotr from another folder:
# import sys
# sys.path.insert(0, "\\w\\PythonU\\")
# import hnumber_theory as nth
# print(nth.Factor(1001))

import math

#-------------------------------------------------------------------------------
def Ferma(n, L):
  ''' Fermat primality test, L - list of bases, 
  Exampler: Ferma(341,[2,3])
  '''
  for a in L:
    if pow(a,n-1,n) != 1: return False
  return True
#------------------------------------------------------------------------------
def Miller_Rabin(n, L):
  ''' Miller Rabin primality test, L - list of bases, 
  Exampler: Miller_Rabin(341,[2,3])
  '''
  if n<2: return False
  if n<128:
    return n in set([ 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,  
                     59,61,67,71,73,79,83,89,97,101,103,107,109,113,127])
  if n%2==0 or n%3==0: return False

  def check(a, s, d, n):
    x = pow(a, d, n)
    if x == 1: return True
    for i in range(s - 1):
      if x == n - 1: return True
      x = pow(x, 2, n)
    return x == n - 1
  # -------
  s = 0; d = n - 1
  while d % 2 == 0:
    d >>= 1
    s += 1
  for a in L:
    if not check(a, s, d, n): return False
  return True

#------------------------------------------------------------------------------
'''Frobenius primality test. See details in arXiv:1807.07249 [math.NT]:
       Evaluation of the Effectiveness of the Frobenius Primality Test.
       (c)Sergei Khashin khash2@gmail.com
   This is a probabilistic method.
   The numbers on which it is wrong are called Frobenius pseudoprimes (FPP).
   Have been proved:
   a) There are no FPP < 2^64.
   b) There are many other limitations on FPP.
   c) There is no known FPP at all.
   There is a hypothesis that FPP does not exist at all.

   Example of usage see at the end of the file.
'''
#---- Jacobi symbol -----------------------------------------------------------
def Jacobi(a,b):
    if b<=1 or (b%2)==0: return 0
    res = 1
    if a<0:
      if b%4==3: res=-1
      a = -a
    if a>=b: a=a%b
    while True:
        if a==0: return 0
        # now 1 <=a < b, b odd
        while a%4 ==0: a//=4
        if a%4==2:
            b8 = b%8
            if b8==3 or b8==5: res = -res
            a = a//2
        if a==1: return res
        # now 1 < a < b, a,b odd
        # J(a,b) -> J(b,a)
        if a%4==3 and b%4==3 : res = -res
        b, a = a, b%a

# ------------------------------------------------------------------------------
def Frobenius_index(n):
    c = 0
    if n % 4 == 3: return -1
    if n % 8 == 5: return 2
    if n % 24 == 17: return 3
    c = 5
    while True:
        jcn = Jacobi(c, n)
        if jcn == 0: return -c
        if jcn == -1: break
        c = nextPrime(c)
    return c

# ------------------------------------------------------------------------------
def is_perfect_square(n):
    x = 1<<((n.bit_length()+1)//2)
    while True:
        delta = x*x-n
        if delta==0: return True
        if delta< 0: break
        x -= 1+delta//(2*x)
    return False

#------------------------------------------------------------------------------
def Frobenius(n):
    '''Frobinius primality test. See details in arXiv:1807.07249 [math.NT]:
           Evaluation of the Effectiveness of the Frobenius Primality Test.
           Sergei Khashin.
    '''
    if n<2: return False
    if n<128:
      return n in set([ 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,  
                       59,61,67,71,73,79,83,89,97,101,103,107,109,113,127])
    if not n&1: return False            # even
    sum_of_bytes = sum(n.to_bytes((7 + n.bit_length()) // 8, 'little'))  # sum of bytes of n
    for p in [3,5,17]:
        if sum_of_bytes%p==0: return False

    if is_perfect_square(n): return False
    # now n is not divisible by 2, 3, 5, 17 and is not a perfect square

    c = Frobenius_index(n)
    if c<-1 and -c != n: return False
    a, b = 1, 1
    if c<3: a=2
    #a1,b1 = Frob_pow(a,b,c,n,n)
    a1, b1 = 1, 0
    a2, b2 = a, b
    k = n
    while k>0:
        if k % 2 == 1:                                          # z1 *= z2
            a1, b1 = (a1*a2 + b1*b2*c) % n, (a1*b2 + a2*b1) % n
        b2, a2 = (2*a2*b2) % n, (a2*a2+c*b2*b2) % n             # z2 *= z2
        k>>=1                                                   # k /= 2
    return (a1==a) and (b1==n-b)


def is_prime(n): return Frobenius(n)

#------------------------------------------------------------------------------
def zPower(a,b,c,k,n):
    ''' (a + b*sqrt(c))^k mod n, retrun (a1,b1) '''
    a1, b1 = 1, 0
    a2, b2 = a, b
    while k>0:
        if k % 2 == 1:  # z1 *= z2
            a1, b1 = (a1 * a2 + b1 * b2 * c) % n, (a1 * b2 + a2 * b1) % n
        b2, a2 = (2 * a2 * b2) % n, (a2 * a2 + c * b2 * b2) % n  # z2 *= z2
        k >>= 1  # k /= 2

    return a1, b1
#------------------------------------------------------------------------------
def nextPrime(n):
    ''' prime greater then n '''
    n+= 1 + (n&1)
    while(not Frobenius(n)): n+=2
    return(n)
#------------------------------------------------------------------------------
def Pollard(n):
    ''' Pollard's Rho Algorithm for Prime Factorization
    Method is effective if there is a divisor less than ~2**50 ~10**15
    * n is integer >2 (really n>=23*23)
    * n does not divide 2, 3, 5, 7, 11, 13, 17, 19
    * n is not prime
    '''
    for a in (1, 2,3,7,11):
        if n>1000_000 and a>1: print('Pollard, n=', n, 'a=', a, '/12')
        x = 2
        for cycle in range(25):   # till 2**24 of gcd calculation
            if cycle > 21: print('           cycle=', cycle, '/25')
            y=x
            for i in range(2 ** cycle):
                x = (x * x + a) % n
                factor = math.gcd(x - y, n)
                if factor == n: y=-2; break # all factors have the same divisor, take the next a
                if factor > 1: return factor
            if y<0: break
    return y if y==-2 else -1

# -----------------------------------------------------------------------------
def _Factor1(n):
    '''
    The Factor1(n) function returns a list of divisors of the number n.
    We assume (without checking) that n>20 and has no divisors <20.
    The function is recursive. Using the function _divider(n)
    finds the divisor d of the number n and calls Factor 1(d), Factor1(n//d)
    '''
    if Frobenius(n): return [n]
    a = Pollard(n)
    if a<=1: return [a]          # не вышло, возвращаем -1
    return _Factor1(a) + _Factor1(n // a)

#------------------------------------------------------------------------------
def small_dividers(n): 	# return the list of prime dividers<20 and n1
    if n<0:
        L = [-1]
        n = -n
    else: L = []

    while not n&1:
        L.append(2)
        n >>= 1
    # now n is odd

    have_div = True
    while have_div:
        have_div = False
        sum_of_bytes = sum(n.to_bytes((7 + n.bit_length()) // 8, 'little'))  # сумма байтов числа
        for p in [3,5,17]:
            if sum_of_bytes%p==0:
                n //= p
                L.append(p)
                have_div = True
    # now n is not divisible by 2,3,5,17

    for p in (7,11,13,19):
        while n%p==0:
            L.append(p)
            n //= p
    # now n is not divisible by 2,3,5,7, 11, 13, 17, 19

    return L, n

#------------------------------------------------------------------------------
def Factor(n): 	# return the list of prime dividers
    assert type(n) == int
    if n==0: return [0]
    if n==1: return [1]
    L, n = small_dividers(n)
    # now n is not divisible by 2,3,5,7, 11, 13, 17, 19
    if n>1:
        L.extend(_Factor1(n))
    L.sort()
    return L


#------------------------------------------------------------------------------
def egcd_recurcive(a,b):
  ''' return 3 numbers (d,u,v), such d=GCD(a,a) and d = u*a + v*a '''
  if a<0: a=-a
  if b<0: b=-b
  if a==0: return(b,0,1)
  d,Y,X= egcd(b%a,a)
  return(d,X-(b//a)*Y,Y)
#------------------------------------------------------------------------------
def egcd(a,b):
  ''' return 3 numbers (d,u,v), such d=GCD(a,a) and d = u*a + v*a '''
  if a<0: a=-a
  if b<0: b=-b
  (x0, x1, y0, y1) = (1, 0, 0, 1)
  while b != 0:
    (q, a, b) = (a // b, b, a % b)
    (x0, x1) = (x1, x0 - q * x1)
    (y0, y1) = (y1, y0 - q * y1)
  return (a, x0, y0)

#------------------------------------------------------------------------------
#def mod_inv(a,n):  	# 1/a mod n  change for  pow(a,-1,n)
#  a=a%n
#  d,X,Y=egcd(a,n)
#  if d!=1:
#      raise Exception("d != 1")
#  return(X%n)
#------------------------------------------------------------------------------
#def mod_div(a,b,n): # a/b mod n  change for (a * pow(b,-1,n))%n
#  a=a%n
#  b=b%n
#  b1 = mod_inv(b,n)
#  return a*b1%n
#-------------------------------------------------------------------------------
def chRem2(a1, n1, a2, n2):			# Chinese remainder theorem
  '''
  return k:
  k = a1 mod n1
  k = a2 mod n2
  '''
  d, u1, u2 = egcd(n1,n2)     	# d = u1*n1 + u2*n2
  k, r = divmod(a2-a1, d)
  if r !=0: return None           # no solutions
  return (a1+k*u1*n1)%(n1//d*n2)

# =============================================================================
if __name__ == "__main__":

    # Fermat primality test
    assert Ferma(7*13*19, [2,3])
    # Miller_Rabin primality test
    assert not Miller_Rabin(7*13*19, [2,3])
    assert Miller_Rabin(3283981*6567961, [2,3,5,7,11])
    print('Ferma+, MR+')

    # Frobenius primality test
    # There are no counterexamples less than 2**64
    # No counterexample is known.
    assert Frobenius(17)
    assert not Frobenius(7*242863)
    assert not Frobenius(3283981*6567961)
    print('Frobenius+')

    n = 53881
    a1, b1 = zPower(a=1, b=1, c=23, k=n, n=n)
    assert a1==1 and b1==n-1
    print('zPower+')

    assert is_perfect_square(1001**2)
    assert not is_perfect_square(1001**2+77)
    assert is_perfect_square(17**8)
    assert not is_perfect_square(17**8+777)
    print('is_perfect_square+')

    # nextPrime
    assert nextPrime(6) == 7
    assert nextPrime(7) == 11

    # Factor. Only Pollard algorithm!
    assert Factor(10**30+1)== [61, 101, 3541, 9901, 27961, 4188901, 39526741]
    L = Factor(2**3 * 3**3 * 5**2 * 7**4 * 101**2 )   #
    assert L == [2, 2, 2, 3, 3, 3, 5, 5, 7, 7, 7, 7, 101, 101]
    print('Factor+')

    print('--->Factor(p*q)')
    for k in range(10):
        p = nextPrime(17**k+100000)
        q = nextPrime(7*p+10**6)
        n = p*q
        assert Factor(n) == [p, q]
        print( 'bit_length=', n.bit_length(), 'Factor+')
