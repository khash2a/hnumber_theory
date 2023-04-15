# hnumber_theory
Simple number theory in Python

Most intresting:
Frobenius primality test. See details in arXiv:1807.07249 [math.NT]:
   Evaluation of the Effectiveness of the Frobenius Primality Test.
   (c)Sergei Khashin khash2@gmail.com
   
   This is a probabilistic method.
   
   The numbers on which it is wrong are called Frobenius pseudoprimes (FPP).
   
   Have been proved:
   
   a) There are no FPP < 2^64.
   
   b) There are many other limitations on FPP.
   
   c) There is no known FPP at all.
   
   There is a hypothesis that FPP does not exist at all.

Content: 

Ferma(n, L)			    Fermat primality test of n by bases in list L

Miller_Rabin(n, L)	M-R primality test of n by bases in list L

Frobenius_index(n)  Frobenius index

is_perfect_square(n)test for perfect square

Frobenius(n)			  Frobenius primality test

is_prime(n):        return Frobenius(n)

Jacobi(a,b)			    Jacobi symbol

zPower(a,b,c,k0,n)  (a + b*sqrt(c))^k0 mod n, return (a1,b1)

nextPrime(n)			  prime, greater then n

Pollard(n)     		  Pollard's Rho Algorithm for Prime Factorization

Factor(n)				    return the list of prime dividers

egcd(a,b)				    return 3 numbers (d,u,v), such d=GCD(a,a) and d = u*a + v*a

chRem2(a1,n1, a2,n2)Chinese remainder theorem
