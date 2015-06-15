# Spherical-Harmonic-Transform
A collection of MATLAB routines for the Spherical Harmonic Transform and related manipulations in the spherical harmonic spectrum.

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

The associated Legendre functions P_l^m DO NOT include the Condon-Shortley phase factor (-1)^m in their definition
(note that MATLAB by default includes this factor only in the unnormalized form of the Legendre function, while it omits it in the normalized forms)

The SHT transform can be done by
a) direct summation, for appropriate sampling schemes along with their
integration weights, such as the uniform spherical t-Designs, the Fliege-Maier
sets, Gauss-Legendre quadrature grids, Lebedev grids and others.
b) least-squares, weighted or not, for arbitrary sampling schemes. In this
case weights can be provided externally, or use generic weights based on the
areas of the spherical polygons around each evaluation point determined by
the Voronoi diagram of the points on the unit sphere, using the included
functions.
MAT-files containing t-Designs and Fliege-Maier sets are also included.
For more information on t-designs, see:

    http://neilsloane.com/sphdesigns/
and
    McLaren's Improved Snub Cube and Other New Spherical Designs in Three
    Dimensions, R. H. Hardin and N. J. A. Sloane, Discrete and Computational
    Geometry, 15 (1996), pp. 429-441.

while for the Fliege-Maier sets see

    http://www.personal.soton.ac.uk/jf1w07/nodes/nodes.html
and
    The distribution of points on the sphere and corresponding cubature
    formulae, J. Fliege and U. Maier, IMA Journal of Numerical Analysis (1999),
    19 (2): 317-334

Most of the functionality of the library is displayed in the included test
scripts.

Some routines in the library evaluate Gaunt coefficients, which express the
integral of the three spherical harmonics. These can be evaluated either
through Clebsch-Gordan coefficients, or from the related Wigner-3j symbols.
Here they are evaluated through the Wigner-3j symbols through the formula
introduced in

    Translational addition theorems for spherical vector wave functions,
    O. R. Cruzan, Quart. Appl. Math. 20, 33?40 (1962)

which can be also found in

    http://mathworld.wolfram.com/Wigner3j-Symbol.html, Eq.17.

Finally, in the library functions are included that compute coefficients of 
rotated functions, either for the simple case of an axisymmetric kernel 
rotated to some direction (\theta_0, phi_0), or the more complex case of 
arbitrary functions were wull rotation matrices are constructed from Euler 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

For a list of functions and scripts check the README.txt file.
For more details on the functions, check their help output in Matlab.
