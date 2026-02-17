"""
Created on Sat May 18 10:12:59 2019

@author: bv20
 
  Wavmat -- Wavelet Transformation Matrix  
   Usage
     W = WavMat(N, h, k0, shift)
   Inputs
     n     size of matrix/length of data. Should be power of 2.
     h      low-pass filter corresponding to orthogonal WT'
            Default is Haar.
     k0     depth of transformation. Ranges from 1 to J=log2(N).
            Default is J. 
    shift  the matrix is not unique an any integer shift gives
           a valid transformation. Default is 2.
  Outputs
    W      n x n transformation matrix 

  Description
    For a quadrature mirror filter h (low pass) the wavelet
    matrix is formed. The algorithm is described in 
    [BV99] Vidakovic, B. (1999). Statistical Modeling By Wavelets, Wiley,
     on pages 115-116.
     Any shift is valid.  Size n=2048 is still managable on a standard PC.
 
   Usage
     We will mimic the example 4.3.1 from [BV99] page 112.
import numpy as np
from Wavmat import Wavmat
haar= [1/np.sqrt(2), 1/np.sqrt(2)] 
dat=np.array([1, 0, -3, 2, 1, 0, 1, 2])
W = Wavmat(8, haar, 3, 2) 
wt = np.dot(W,dat) 
print(wt)
[ 1.41421356 -1.41421356  1.         -1.          0.70710678 -3.53553391
  0.70710678 -0.70710678]





%
"""


def Wavmat(n, h=None, k0=None, shift=None): 
    import numpy as np    
    if h is None:
          h=np.array([1/np.sqrt(2), 1/np.sqrt(2)])
    if k0 is None:
          k0=int(np.log2(n))
    if shift is None:
          shift=2 
    h = np.array(h)
    flag = float 
    if sum(np.iscomplex(h) > 0):
        flag = np.dtype(complex)
    g = np.conj(h[::-1])*[(-1)**i for i in range(len(h))]
    h=np.concatenate((h, np.zeros(n))) #extend filter H by 0's to sample by modulus
    g=np.concatenate((g, np.zeros(n))) #extend filter G by 0's to sample by modulus
    
    J=int(np.log2(n))
      
    oldmat = np.identity(2**(J-k0)) 
    print(oldmat)
    for k in np.arange(k0,0,-1):   
                ubJk =  2**(J-k) 
                ubJk1 = 2**(J-k+1) 
                hmat = np.zeros((ubJk1, ubJk), dtype = flag)
                gmat = np.zeros((ubJk1, ubJk), dtype = flag)
                for  jj in np.arange(ubJk):
                    for ii in np.arange(ubJk1):
                          modulus =  (n+ii-2*jj+shift-1) % ubJk1 
                          modulus = modulus + int(modulus == 0)*ubJk1 
                          hmat[ii][jj] =  h[ modulus-1 ] 
                          gmat[ii][jj] =  g[ modulus-1 ] 
                W = np.vstack((np.dot(oldmat, hmat.T ), gmat.T))
                oldmat = W
    return(W)
