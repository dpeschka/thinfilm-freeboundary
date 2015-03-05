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

## Some simple experiments with the algorithm

Each experiment here assumes that you start with the other parameters in the main file *thinfilm.m* being as stated above. The intent of these examples is to show the versatility of the method and give some more intuition for the physics/mathematics of thin-film contact line motion.
**Note:** If one experiences problems with stability, then this is usually due to an (expected) restriction in the time-step size. Solution: Increase nt or decrease T!

1. Change the number of vertices to 
  * a small number, e.g. npoint=10
  * a big number, e.g. npoint=1000. This requires increasing the number of time steps as well, e.g. nt=1000
  
  **Result:** This shows the robustness of the algorithm.
  
2. Switch on normal gravity: Set c1=100.
  
  **Result:** The original droplet becomes flat.
  
3. Switch on tangential gravity:
  * first set setting c2=10 and set T=1, nt=1000
  * drastically increase c2=100 and set T=1, nt=4000
  
  **Result:** The first setting produces a symmetric traveling wave h(t,x)=f(x-tv), whereas the second setting produces a very asymmetric traveling wave.
4. Gradient in surface energy: Set SR=0.1, nt=800, npoint=200
  
  **Result:** Also generates a traveling wave with contact angle left and right being different. This reflects a gradient in surface energy, by which the droplet starts to move towards the smaller angle.
  
## Slightly advanced experiments

1. Modification of parameters and initial data to have *dewetting like* behavior:
  Set L=50, T=200, nt=500 and after line 25 insert *h(h>1)=1* so that the corresponding lines look like
  ```
  % * create & remember initial data
  h  = L/2-abs(L/2-x); 
  h(h>1)=1;
  ```
  **Result:** Typical dewetting front with volume collected in a rim with a droplet being the final state. 
  
2. Modification of parameters for wetting/droplet spreading case with zero contact angle
  Set SL=SR=0 and increase nt=1000.
  
  **Result:** After running the program also check the first derivative of *h* using finite differences via
  ```
  xh = (x(1:end-1)+x(2:end))/2;
  dh = diff(h)./diff(x);
  plot(xh,dh);
  ```
  to investigate the smoothness (and the actual contact angle) of *h*.
