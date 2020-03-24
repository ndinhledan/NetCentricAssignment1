import matplotlib.pyplot as plt
import numpy as np
from math import e
import math

a = 0.01 # interval
p = 0.5
intervals = [(0, 10)]

def plot_slotted_aloha(r1, r2):
    g = np.linspace(r1, r2, 400)

    slotted_aloha = g * (e**(-g))
    plt.plot(g, slotted_aloha, 'r', label="Slotted Aloha")

def plot_pure_aloha(r1, r2):
    g = np.linspace(r1, r2, 400)

    pure_aloha = g * (e**(-g*2))
    plt.plot(g, pure_aloha, 'g', label="Pure Aloha")

def plot_non_persistent_csma(r1, r2):
    g = np.linspace(r1, r2, 400)

    non_persistent_csma = (g*(e**(-a*g))) / (g*(1 + 2*a) + e**(-a*g)) 
    plt.plot(g, non_persistent_csma, 'b', label="Nonpersistent CSMA")

def plot_one_persistent_csma(r1, r2):
    g = np.linspace(r1, r2, 400)

    one_persistent_csma = (g * (1 + g + a*g*(1 + g + a*g/2))*e**(-g * (1 + 2*a))) / (g*(1 + 2*a) - (1 - e**(-a*g)) + (1 + a*g)*e**(-g*(1+a)))
    plt.plot(g, one_persistent_csma, 'c', label="1-Persistent CSMA")

# def sigma(start, end):


def pi_n(g, n):
    return ((((1 + a)*g)**n)/math.factorial(n)) * e**(-g*(1 + a))

def plot_p_persistent_csma(r1, r2):
    g = np.linspace(r1, r2, 400)

    _ps = 0
    pi_0 = pi_n(g, 0)
    ps = 0
    _t = 0
    t = 0

    # p_persistent_csma = ((1 - e**(-a*g))*(_ps*pi_0 + ps * (1 - pi_0))) / ((1 - e**(-a*g))*(a*_t*pi_0 + a*t*(1 - pi_0) + 1 + a) + a*pi_0)
    
    p_persistent_csma = ((a + p)*g*(e**(-g*(a+p))) - p*g*(e**(-g*(2*a + p)))) / ((1 + a)*(1 - e**(-a*g)) + a*e**(-g*(a+p)))
    plt.plot(g, p_persistent_csma, 'y', label="p-Persistent CSMA")
    

def main():

    for start, end in intervals:
        plot_slotted_aloha(start, end)
        plot_pure_aloha(start, end)
        plot_non_persistent_csma(start, end)
        plot_one_persistent_csma(start, end)
        plot_p_persistent_csma(start, end)

    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0.)

    plt.show()

if __name__ == "__main__":
    main()