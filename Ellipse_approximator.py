# -*- coding: utf-8 -*-
import math as m

def Is_Even(n):
    '''RETURNS TRUE WHEN n IS EVEN'''
    half = int(n / 2)
    return abs((n/2) - half) < 0.25

def Double_Factorial(n):
    '''CALCULATES THE DOUBLE FACTORIAL OF AN INTEGER n'''
    if n < 0 or type(n) != int:
        raise(ValueError)
    elif n == 0:
        return 1
    else:
        even = Is_Even(n)
        d_fact = 1
        i = 0
        while (n-2*i)>0:
            d_fact = d_fact * (n-2*i)
            i += 1
        return d_fact

def Binomial_Coef(n, k):
    '''CALCULATES THE BINOMIAL COEFFICIENT (n choose k)'''
    if n >= k:
        float_comb = m.factorial(n) / (m.factorial(k) * m.factorial(n-k))
        return int(float_comb)
    else:
        raise(ValueError)

def Epsilon_Series(eps, n, stirling_lim = 10**9):
    '''CALCULATES THE SERIES EXPANSION OF THE COMPLETE ELLIPTIC INTEGRAL OF THE
    SECOND KIND WITH ECCENTRICITY eps UP TO n TERMS. USES STIRLING'S FORMULA
    AFTER stirling_lim TERMS TO MAKE CALCULATIONS EASIER.'''
    Circ = m.pi / 2
    series = 0
    for k in range(1,n+1):
        
        if k <= stirling_lim:
            square_term = (Double_Factorial(2*k-1) / Double_Factorial(2*k))**2
        else:
            square_term = ((2*k - 1) / (2*m.pi*k*(k-1))) * ((2*k-1)**(4*k-2) / ((2*k-2)**(2*k-2) * (2*k)**(2*k)))
                                                            
        series += square_term * (eps**(2*k) / (1 - 2*k))
    return Circ * (1 + series)

def Eta_Series(eta, inner_limit, outer_limit, stirling_lim = 10**9):
    '''CALCULATES THE PERTURBATION SERIES OF THE PERIMETER OF AN ELLIPSE, SEEN AS
    A PERTURBATION FROM A LINE. USES STIRLING'S FORMULA AFTER stirling_lim TERMS
    BOTH IN THE INNER SERIES AND IN THE OUTER SERIES.'''
    a0 = 1
    series = 0
    for n in range(1, outer_limit + 1):
        inner_series = 0
        for k in range(n, n + inner_limit + 1):
            if k <= stirling_lim:
                square_term = (Double_Factorial(2*k-1) / Double_Factorial(2*k))**2
            else:
                square_term = ((2*k - 1) / (2*m.pi*k*(k-1))) * ((2*k-1)**(4*k-2) / ((2*k-2)**(2*k-2) * (2*k)**(2*k)))
            frac = 1 / (1 - 2*k)        
            comb = Binomial_Coef(k, n)
            inner_series += square_term * comb * frac

        series += (eta**(2*n)) * (-1)**(n) * inner_series

    return a0 + (m.pi / 2)*series


print('The circumference of a circle (eps = 0, eta = 1) of radius 1 is')
print(2*m.pi)
print('The Epsilon_Series * 4 with eps = 0.9 for 10000 terms is')
print(4.686788211126458)
print('The Eta_Series * 4 with eta = 0.1 and 800 inner terms and 50 outer terms is')
print(4.062724810623535)
print('The Eta_Series * 4 with eta = 0.1 and 50 inner terms and 800 outer terms is')
print(4.050649403358031)




















