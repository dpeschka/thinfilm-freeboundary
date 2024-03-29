# A MATLAB algorithms for thin-film free boundary problems in one spatial dimension

This repository contains various evolutions of MATLAB codes that solve the thin-film problem

```math
\partial_t h - \partial_x (m(h)\partial_x \pi)=0, \qquad \pi=\frac{\delta F(h)}{\delta h}
```

with a general free energy of the form

```math
F(q)=\int_{\{x:h>0\}} \frac{1}{2}(\partial_x h)^2 + (-S(x)) + g_2 h^2 - g_1 hx\,{\rm d}x
```

encoding surface energies, spreading coefficient, normal gravity, tangential gravity in the integrand, respectively. Different codes vary in the methodology that is used for the treatment of the contact angle. We have

* `thinfilm.m`: MATLAB code for 1D thin film equation with contact lines (static contact angle) as discussed in the corresponding paper in *Thin-film free boundary problems for partial wetting* published in *Journal of Computational Physics* available under [http://dx.doi.org/10.1016/j.jcp.2015.04.041](http://dx.doi.org/10.1016/j.jcp.2015.04.041). The corresponding code is explained in detail below and `thinfilm_clm.m` and `thinfilm_clm_dual.m` mainly follow this approach.
<br>

* `thinfilm_clm.m`: MATLAB code for the 1D thin film equation with contact lines (dynamic contact angle) as discussed in the corresponding paper *Variational approach to dynamic contact angles for thin films* in *Physics of Fluids* available under [https://doi.org/10.1063/1.5040985](https://doi.org/10.1063/1.5040985)
<br>

* `thinfilm_clm_dual.m`: MATLAB code for the 1D thin film equation with dynamic contact lines, general mobility exponent in a dual variational formulation as used in the preprint 

  1. *Droplet motion with contact-line friction: long-time asymptotics in complete wetting* by Lorenzo Giacomelli, Manuel V. Gnann, Dirk Peschka
  [https://arxiv.org/abs/2302.03005](https://arxiv.org/abs/2302.03005)


## Explanation of the main file *thinfilm.m* ##
#### I. Introduction
<img src="https://github.com/dpeschka/thinfilm-freeboundary/blob/master/pics/example.png" width="40%">

This example solves the thin-film equation with mobility exponent *n=2* and initial data *h0(x)=1/2-|x-1/2|* for *0<x<1*. The equilibrium contact angles at the left and right side are both *|h'|=sqrt(2)* as we have *SL=SR=1*. The following parameters can be modified by the user:

#### II. Parameters

```matlab
L     = 1.0;  % initial domain size (0,L)
T     = 0.2;  % final time
SL    = 1.0;  % negative spreading coefficient at x=x-
SR    = 1.0;  % negative spreading coefficient at x=x+
g1    = 0.0;  % tangential gravity
g2    = 0.0;  % normal gravity
nt    = 100;  % number of time steps
npoint= 100;  % number of vertices
```

The initial domain is *(0,L)* as set by the parameter `L` and also incorporated in the initial data *h0*. One needs to have *h0(0)=h0(L)=0* and *h0(x)>0* for *0<x<L*. Then, the algorithm attempts to solve then thin-film equation for times *0<t<T*, which might fail due to

  * topological changes (unavoidable),
  * numerical instability (decrease time-step).

The user can change the contact angles at *x_+/-* by setting the values of `SL,SR` so that *|h'(x-)|=sqrt(2SL)* and *|h'(x_+)|=sqrt(2SR)*. The parameter `g1,g2` encode normal and tangential gravity. The number of time-steps is `nt`, so that `dt=T/nt`. The initial spatial resolution is `L/npoint`, which, however, will change during the evolution. Since the domain deformation is linear, the spacing/decomposition will always stay uniform.

#### III. FEM specifics

This part constructs the standard finite element infrastructure.
  * decomposition of interval (0,L) into `nelement` intervals using `npoint` vertices `x`
  * infrastructure `nd` stores the 2 vertices, attached to an element (easy in 1D since `x` ordererd)
  * `local_mass_p1` stores the mass matrix `M_ij=\int \phi_i(x)\phi_j(x) dx` for phi1(x)=1-x, phi2(x)=x for 0<x<1

and creates the initial data *h0(x)*.


```matlab
% * create element decomposition for FE method
x               =linspace(0,L,npoint)';% vertices
nelement        =npoint-1;             % no elements
nd(1:nelement,1)=1:npoint-1;           % id left point of an element
nd(1:nelement,2)=2:npoint;             % id right point of an element
local_mass_p1   =[1/3 1/6;1/6 1/3];    % mass matrix for reference [0,1]

% * create & remember initial data
h  = L/2-abs(L/2-x); 
```

#### IV. Main PDE part

```matlab
for it=1:nt       
    % construct system matrices
    build_FE_matrices % script: matrices A,S,M,Dx for FEM
    build_ALE_matrix  % script: matrices for ALE decomposition
    
    % FE problem: build right-hand-side rhs & solve
    rhs=[zeros(npoint,1);S*h+M*(2*g2*h-g1*x)];
    rhs(ndof+1)=rhs(ndof+1)+(SL+(dh(  1)^2)/2)/dh(  1);
    rhs(2*ndof)=rhs(2*ndof)-(SR+(dh(end)^2)/2)/dh(end);
    u = A\rhs; % solve for u=(hdot,pi)^t
    
    % perform ALE decomposition & update solution
    U = (I-P)\u(1:ndof); % select only u, forget p
    h = h + dt*I*U;      % update h
    x = x + dt*X*U;      % update x
    
    % plot numerical solution
    if mod(it,10)==1
        plot(x,h,'b-','LineWidth',2);
        drawnow
    end
end
```

Steps in the main part:

  * `build_FE_matrices.m`: Builds the matrices `A,I,X,P,S,M` as in the paper in 
  Since this is standard finite element calculus, just a remark the the final matrix `A = [M  Sw;-dt*S M];` is    computed as described in the paper.
  
  * `build_ALE_matrix.m`: Construction of matrices related to the ALE decomposition (straight-forward)
  **Remark:** The computation of the space derivative of *h* in the weak form `dh = M\(Dx*h);` uses the matrices `M,Dx` contructed in the finite element part.
  * Create right-hand-side *b* of the problem *Ax=b* as described with surface tension `S*h` and gravity `M*(2*g2*h-g1*x)`. The modifications at the end-points account for the contact angle due to a nonzero spreading coefficient.
  * Solve: `u = A\rhs`
  * ALE decomposition: `U = (I-P)\u(1:ndof);`
  * Update solution: `h = h + dt*I*U;` and `x = x + dt*X*U;`
 
## Some simple experiments with the algorithm

Each experiment 1.,2.,3.,4. assumes that one starts with the standard parameters in the main file *thinfilm.m* being as stated above. The intent of these examples is to show the versatility of the method and give some more intuition for the physics/mathematics of thin-film contact line motion.
**Note:** If one experiences problems with stability, then this is usually due to an (expected) restriction in the time-step size. Solution: Increase nt or decrease T!

1. Change the number of vertices to: 
  * a small number, e.g. `npoint=10`,
  * a big number, e.g. `npoint=1000`. This requires increasing the number of time steps as well, e.g. `nt=1000`.
  
  **Result:** This shows the robustness of the algorithm.
  
2. Switch on normal gravity: Set `g2=50`.
  
  **Result:** The original droplet becomes flat.
  
3. Switch on tangential gravity:
  * first set setting `g1=10` and set `T=1`, `nt=1000`
  * drastically increase `g1=100` and set `T=1`, `nt=4000`
  
  **Result:** The first setting produces a symmetric traveling wave h(t,x)=f(x-tv), whereas the second setting produces a very asymmetric traveling wave.
4. Gradient in surface energy: Set `SR=0.1`, `nt=800`, `npoint=200`.
  
  **Result:** Also generates a traveling wave with contact angle left and right being different. This reflects a gradient in surface energy, by which the droplet starts to move towards the smaller angle.
  
## Slightly advanced experiments

1. Modification of parameters and initial data to have *dewetting like* behavior:
  Set `L=50`, `T=200`, `nt=500` and after line 25 insert `h(h>1)=1` so that the corresponding lines look like
  ```matlab
  % * create & remember initial data
  h  = L/2-abs(L/2-x); 
  h(h>1)=1;
  ```
  **Result:** Typical dewetting front with volume collected in a rim with a droplet being the final state. 
  
2. Modification of parameters for wetting/droplet spreading case with zero contact angle
  Set `SL=0`, `SR=0` and increase `nt=1000`.
  
  **Result:** After running the program also check the first derivative of *h* using finite differences via
  ```matlab
  xh = (x(1:end-1)+x(2:end))/2;
  dh = diff(h)./diff(x);
  plot(xh,dh);
  ```
  to investigate the smoothness (and the actual contact angle) of *h*.
