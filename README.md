# A MATLAB algorithm for thin-film free boundary problems in one spatial dimension
MATLAB code for 1D thin film equation with contact lines as discussed in the corresponding paper "" in *Journal of Computational Physics*.

## Explanation of the main file *thinfilm.m* ##

This example solves the thin-film solution with mobility exponent *n=2* and initial data h_0(x)=1/2-|x-1/2| for 0<x<1. The equilibrium contact angles at the left and right side are both |h'|=sqrt(2). The following parameters can be modified by the user:

```
L     = 1.0;  % initial domain size (0,L)
T     = 0.2;  % final time
SL    = 1.0;  % negative spreading coefficient at x=x-
SR    = 1.0;  % negative spreading coefficient at x=x+
c1    = 0.0;  % normal gravity
c2    = 0.0;  % tangential gravity
nt    = 100;  % number of time steps
npoint= 100;  % number of vertices
```

Some simple experiments that can be done with the algorithm is

1. Change the number of vertices to 
  * a small number, e.g. 10
  * a big number, e.g. 1000. This requires increasing the number of time steps as well
  
  **Result:** This shows the robustness of the algorithm.
  
2. Switch on normal gravity: Set c1=100.
  
  **Result:** The original droplet becomes flat.
  
3. Switch on tangential gravity:
  * first set setting c2=10 and set T=1, nt=1000
  * drastically increase c2=100 and set T=1, nt=4000
  
  **Result:** The first setting produces a symmetric traveling wave h(t,x)=f(x-tv), whereas the second setting produces a very asymmetric traveling wave.
4. Gradient in surface energy: Set SR=0.1, nt=800, npoint=200
  
  **Result:** Also generates a traveling wave with contact angle left and right being different. This reflects a gradient in surface energy, by which the droplet starts to move towards the smaller angle.
