# Spherical Harmonic Transform Library

#### A collection of MATLAB routines for the Spherical Harmonic Transform and related manipulations in the spherical harmonic spectrum.

---
>
>    Archontis Politis, 2015  
>
>    Department of Signal Processing and Acoustics, Aalto University, Finland  
>
>    archontis.politis@aalto.fi
>
---

This Matlab/Octave library was developed during my doctoral research in the [Communication Acoustics Research Group] (http://spa.aalto.fi/en/research/research_groups/communication_acoustics/), Aalto University, Finland. If you would like to reference the code, you can refer to my dissertation published [here](https://aaltodoc.aalto.fi/handle/123456789/22499):

    Archontis Politis, Microphone array processing for parametric spatial audio techniques, 2016
    Doctoral Dissertation, Department of Signal Processing and Acoustics, Aalto University, Finland
    
## Description

Both real and complex SH are supported. The orthonormalised versions of SH
are used. More specifically, the complex SHs are given by:

    Y_{nm}(\theta,\phi) =
    (-1)^m \sqrt{\frac{2n+1}{4\pi}\frac{(n-m)!}{(n+m)!}} P_l^m(\cos\theta) e^{im\phi}

and the real ones as in:

    R_{nm}(\theta,\phi) = 
    \sqrt{\frac{2n+1}{4\pi}\frac{(n-|m|)!}{(n+|m|)!}} P_l^{|m|}(\cos\theta) N_m(\phi)
  
where

    N_m(\phi) = \sqrt{2} cos(m\phi},    m>0
    N_m(\phi) = 1,    m>0
    N_m(\phi) = \sqrt{2} sin(|m|\phi},  m<0

Note that the Condon-Shortley phase of (-1)^m is not introduced in the code for
the complex SH since it is included in the definition of the associated 
Legendre functions in Matlab (and it is canceled out in the code of the real SH).

The functionality of the library is demonstrated in detail in [http://research.spa.aalto.fi/projects/sht-lib/sht.html]
or in the included script TEST_SCRIPTS_SHT.m.

The SHT transform can be done by:

a) direct summation, for appropriate sampling schemes along with their
integration weights, such as the uniform spherical t-Designs, the Fliege-Maier
sets, Gauss-Legendre quadrature grids, Lebedev grids and others.

b) least-squares, weighted or not, for arbitrary sampling schemes. In this
case weights can be provided externally, or use generic weights based on the
areas of the spherical polygons around each evaluation point determined by
the Voronoi diagram of the points on the unit sphere, using the included
functions.

MAT-files containing t-Designs and Fliege-Maier sets are also included.
For more information on t-designs, see [http://neilsloane.com/sphdesigns/](http://neilsloane.com/sphdesigns/) and

    McLaren's Improved Snub Cube and Other New Spherical Designs in Three
    Dimensions, R. H. Hardin and N. J. A. Sloane, Discrete and Computational
    Geometry, 15 (1996), pp. 429-441.

while for the Fliege-Maier sets see [http://www.personal.soton.ac.uk/jf1w07/nodes/nodes.html](http://www.personal.soton.ac.uk/jf1w07/nodes/nodes.html) and

    The distribution of points on the sphere and corresponding cubature
    formulae, J. Fliege and U. Maier, IMA Journal of Numerical Analysis (1999),
    19 (2): 317-334

Some routines in the library evaluate Gaunt coefficients, which express the
integral of the three spherical harmonics. These can be evaluated either
through Clebsch-Gordan coefficients, or from the related Wigner-3j symbols.
Here they are evaluated through the Wigner-3j symbols through the formula
introduced in

    Translational addition theorems for spherical vector wave functions,
    O. R. Cruzan, Quart. Appl. Math. 20, 33:40 (1962)

which can be also found in [http://mathworld.wolfram.com/Wigner3j-Symbol.html, Eq.17.](http://mathworld.wolfram.com/Wigner3j-Symbol.html)

Finally, a few routines are included that compute coefficients of 
rotated functions, either for the simple case of an axisymmetric kernel 
rotated to some direction (\theta_0, \phi_0), or the more complex case of 
arbitrary functions were full rotation matrices are constructed from Euler 
angles. The algorithm used is according to the recursive method of Ivanic and 
Ruedenberg, as can be found in

    Ivanic, J., Ruedenberg, K. (1996). Rotation Matrices for Real 
    Spherical Harmonics. Direct Determination by Recursion. The Journal 
    of Physical Chemistry, 100(15), 6342?6347.

and with the corrections of

    Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real 
    Spherical Harmonics. Direct Determination by Recursion Page: Additions 
    and Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100.

Rotation matrices for both real and complex SH can be obtained.

For any questions, comments, corrections, or general feedback, please
contact archontis.politis@aalto.fi

---

For more details on the functions, check their help output in Matlab.

### List of MATLAB files:

* getSH.m                   -   Get SHs up to order N
* legendre2.m               -   Same as matlab Legendre extended for negative orders m<0
* leastSquaresSHT.m         -   SHT of function using least-squares
* weightedLeastSquaresSHT.m -   SHT using weighted least-squares
* directSHT.m               -   Perform SHT directly by summation for special arrangements
* leastSquaresSHT.m         -   Perform SHT by least-squares, weighted or unweighted
* inverseSHT.m              -   Perform the inverse SHT
* getTdesign                -   Returns the spherical coordinates of T-designs up to t=21
* getFliegeNodes            -   Returns FliegeMaier point set, up to order N=29 SHT

* t_designs_1_21.mat          -   Tables of t-designs, up to t=21
* fliegeMaierNodes_1_30.mat   -   Tables of Fliege-MAier sets, up to 300 points

* complex2realCoeffs.m      -   Convert SH coeffs from the complex to real basis
* complex2realCoeffsMtx.m    -   Returns complex to real SH basis transformation matrix
* real2complexCoeffs.m      -   Convert SH coeffs from the real to complex basis
* rotateCoeffs.m            -   Get SH coefficients for a rotated axisymmetric 
                                pattern
* conjCoeffs.m              -   Get the complex SH coefficients of a conjugate 
                                function
* sphConvolution.m          -   Perform spherical convolution between a function 
                                and a filter, in the SH domain
* sphMultiplication.m       -   Computes coefficients of product of two spherical 
                                functions, through Gaunt Coefficients
* getSHrotMtx               -   Obtain rotation matrix for full rotation of the 
                                coordinate system, that when applied to the SH 
                                coefficients, returns the coefficients of the 
                                rotated function

* gaunt_mtx.m   -	Construct a 3D matrix of Gaunt coefficients up to three 
                    orders, each one a separate dimension of the matrix
* w3j           -   Evaluate the Wigner-3j symbol through the Racah formula
* w3j_stirling  -   Evaluate the Wigner-3j symbol through the Racah formula, 
                    using Stirling's large factorial approximation
* sym_w3j       -   Returns a Wigner-3j in symbolic form
* wignerD       -   Returns the Wigner-D and wigner-d matrices for rotation of 
                    complex spherical harmonics

* plotSphFunctionGrid     - Plot easily spherical function defined on a regular 
                            grid
* plotSphFunctionTriangle - Plot easily spherical function defined on an 
                            irregular grid
* plotSphFunctionCoeffs   - Plot spherical function with known SH coefficients

* Fdirs2grid.m          -	Helper function for plotting, used with grid2dirs
* grid2dirs.m           -	Construct a vector of regular grid points
* sphDelaunay.m         -	Computes the Delaunay triangulation on the unit sphere
* sphVoronoi.m          -	Computes the a Voronoi diagram on the unit sphere
* sphVoronoiAreas.m     -   Computes the areas of a voronoi diagram on the
                            unit sphere
* getVoronoiWeights.m   -   Conveniently get voronoi weights for a sampling scheme
* checkCondNumberSHT.m  -   Computes the condition number of an sampling scheme
                            for a least-squares SHT
* euler2rotationMatrix  -   Euler angles to rotation matrix
* unitCart2sph          -   Get directly azimuth and elevation from unit vectors in matrix form
* unitSph2cart          -   Get directly unit vectors from azimuth and elevation in matrix form
* replicatePerOrder     -   Replicate l^th element 2*l+1 times across specified dimension
