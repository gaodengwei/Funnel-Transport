# Funnel-Transport

code of Funnel Transport: Approximating Finite-Horizon Reachable Sets

======================I. Toolboxes Required:================================= 


1. spotless:	 

URL: https://github.com/spot-toolbox/spotless

Tips: A fork of the Systems Polynomial Optimization Toolbox. Polynomial computation in the FT is supported by this toolbox. 

2. mosek: 

URL: http://www.mosek.com

Tips: MOSEK is a highly efficient commercial solver of LPs, QPs, SOCPs, SDPs, and MIPs. This is a basic element of the following toolboxes for solving optimization problems. The toolbox "sedumi" also works well for SDPs but is too slower than mosek.

3. Yalmip:	 

URL: https://github.com/yalmip/YALMIP

Tips: Yalmip is a MATLAB toolbox for optimization modeling.

4. Gloptipoly 3（optional）: 

URL: https://homepages.laas.fr/henrion/software/gloptipoly3/

Tips: Gloptipoly 3 is intended to solve, or at least approximate, the Generalized Problem of Moments (GPM), an infinite-dimensional optimization problem that can be viewed as an extension of the classical problem of moments.
 
5. Drake（optional）:	 

URL: https://drake.mit.edu/

Tips: Drake is a model-based design and verification for the robotics toolbox. LQR-trees algorithm [1] and Funnel libraries [2] algorithm are open sources in this toolbox.

6. GPOPS-II（optional）:	 

URL: http://www.gpops2.com/

Tips: GPOPS-II is the next generation of general-purpose optimal control software which provide the nominal trajectory for funnel approximation if required.

======================II. Code Introductions:================================= 

The main code of reference <Korda, M. (2020). Computing controlled invariant sets from data using convex optimization. SIAM Journal on Control and Optimization, 58(5), 2871-2899.> and the main code of reference <V. Magron(2019). Semidefinite approximations of reachable sets for discrete-time polynomial systems. SIAM Journal on Control and Optimization, 57(4).> are provided here. Funnel approximation by using FT method and SOS-based method are also provided here applying for examples "Pendulum Swing Up"(FT is faster) and "ROA approximation"(SOS is faster).

1. A moment matrix computation for inner approximation/outer approximation/approximation for a Semi-algebraic Sets constraint(a unit square).
This is the code in the second time response. We have shown the different moment matrix approximations only in state-space X. It is simple and mature. 
 
2. Code of reference-Computing controlled invariant sets from data using convex optimization
This is the code of paper <Korda, M. (2020). Computing controlled invariant sets from data using convex optimization. SIAM Journal on Control and Optimization, 58(5), 2871-2899.> which tells the truth that an "MCI" is different from a "funnel".

3. Code of reference-Semidefinite approximations of reachable sets for discrete-time polynomial systems
This is the code of paper <V. Magron(2019). Semidefinite approximations of reachable sets for discrete-time polynomial systems. SIAM Journal on Control and Optimization, 57(4).> which tells the truth that an "infinite-horizon RS" is different from the "union of finite-horizon RS(funnel)". As a comparison, we also provide a code demo of FT with Dirac measure approximation.

4. typical code of FT

1) Solving example A by FT and SOS respectively;
2) Solving example B by FT and SOS respectively; 
3) ROA approximation examples C (Including a ROA approximation by different SOS-based methods and a main code of FT.)
4) Funnel of Pendulum Swing Up example D(Including a funnel approximation by SOS-based method forked from drake toolbox and a main code of FT. )


======================III. Help:================================= 

Please add util folder and install required toolboxes above-mentioned!!!!!









