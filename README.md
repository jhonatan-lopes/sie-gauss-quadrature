# SIE Gauss Quadratures

*MATLAB implementation for Gauss-Chebyshev and Gauss-Laguerre quadratures to solve Singular Integral Equations.*

This project presents a simple implementation for the Gauss-Chebyshev and Gauss-Laguerre quadratures applied to solving Singular Integral Equations. 

## Gauss-Chebyshev Quadrature

The Gauss-Chebyshev quadrature is applied to integrals posed over a finite range. The end-points are not included in the quadrature (opposed to Gauss-Jacobi, for instance). For a general singular integral equation of the first kind with a simple Cauchy kernel can be written in standard form as (problem posed over a finite region, e.g. a finite contact lying over [-1,1]):

<img src="https://latex.codecogs.com/svg.latex?F(t)&space;=&space;\int_{-1}^{1}\frac{w(s)p(s)}{t-s}\mathrm{d}s" title="F(t) = \int_{-1}^{1}\frac{w(s)p(s)}{t-s}\mathrm{d}s" />

where *B(s) = w(s)p(s)* is the function that we wish to find, decomposed into *w(s)*, a weight function determined by the quadrature, and *p(s)*, a smooth continuous function. After decomposing the *B(s)*, the problem becomes finding the continous function *p(s)*.

The quadrature discretises the problem as:

<img src="https://latex.codecogs.com/svg.latex?F(t_{k})&space;=&space;\sum_{i=1}^{N}\frac{W_{i}\,p(s_i)}{t_k-s_i},&space;\qquad&space;k=1,...,N&plus;n" title="F(t_{k}) = \sum_{i=1}^{N}\frac{W_{i}\,p(s_i)}{t_k-s_i}, \qquad k=1,...,N+n" />

where *t<sub>k</sub>* are the integration points, *s<sub>i</sub>* are the collocation points, and *W<sub>i<sub>* the weights for each collocation point &mdash; all determined by the quadrature.

The function `[s,t,W] = GCHEB_POINTS(N,cas)` returns the vectors of collocation points `s`, integration points `t` and weights `W` for a given number of quadrature points `N` and the expected behaviour of the sought function, represented by the case integer `cas`.

This quadrature considers four cases (see [[1]], [[2]]) according to the desired behaviour of the sought function at the ends of the integration interval:

Case  | Behaviour at -1 | Behaviour at +1
-------|-----------------|----------------
1 (I)  |     Singular    |     Singular
2 (II) |     Singular    |     Bounded
3 (III)|     Bounded     |     Singular
4 (IV) |     Bounded     |     Bounded

### Krenk's interpolation

The solution of the problem is given at specific integration points *t<sub>k</sub>*. If we wish to find the function value at different points inside the Cauchy region, it is necessary to use a suitable interpolation method. Krenk's interpolation is provided for this.

`[phi1,phim1] = KRENK1(phi,cas)` calculates the values of the smooth function `phi` (B = w*phi) at both ends of the interval (`phi1` for t=1, `phim1` for t=-1). 

`phik = KRENK2(phi,sint,cas)` calculates the values of the smooth function `phi` (B = w*phi) at the points `sint` inside the Cauchy region. 

## Gauss-Laguerre Quadrature

The Gauss-Laguerre quadrature is applied to integrals posed over a semi-infinite range [[3]]. For a general singular integral equation of the first kind with a simple Cauchy kernel can be written in standard form as (problem posed over a semi-infinite region, e.g. a finite contact lying over [0,&infin;)):

<img src="https://latex.codecogs.com/svg.latex?F(t)&space;=&space;\int_{0}^{\inf}\frac{w(s)p(s)}{t-s}\mathrm{d}s" title="F(t) = \int_{-1}^{\inf}\frac{w(s)p(s)}{t-s}\mathrm{d}s" />

where *B(s) = w(s)p(s)* is the function that we wish to find, decomposed into *w(s)*, a weight function determined by the quadrature, and *p(s)*, a smooth continuous function. After decomposing the *B(s)*, the problem becomes finding the continous function *p(s)*.

The quadrature discretises the problem as:

<img src="https://latex.codecogs.com/svg.latex?F(t_{k})&space;=&space;\sum_{i=1}^{N}\frac{W_{i}\,p(s_i)}{t_k-s_i},&space;\qquad&space;k=1,...,N&plus;n" title="F(t_{k}) = \sum_{i=1}^{N}\frac{W_{i}\,p(s_i)}{t_k-s_i}, \qquad k=1,...,N+n" />

where *t<sub>k</sub>* are the integration points, *s<sub>i</sub>* are the collocation points, and *W<sub>i<sub>* the weights for each collocation point &mdash; all determined by the quadrature.

The function `[s,t,W] = GLAG_POINTS_CALC(N,cas)` returns the vectors of collocation points `s`, integration points `t` and weights `W` for a given number of quadrature points `N` and the expected behaviour of the sought function, represented by the case integer `cas`.

The weight function is taken as *w(s)=s<sup>a</sup>e<sup>-s</sup>*. This quadrature considers two cases (see [1]) according to the desired behaviour of the sought function at zero (always bounded at infinity):

Case  | Behaviour at 0 | Weight function
-------|-----------------|----------------
1 (I)  |     Bounded    |     w(s)=s<sup>+1/2</sup>e<sup>-s</sup>
2 (II) |     Singular    |     w(s)=s<sup>-1/2</sup>e<sup>-s</sup>

While the Gauss-Chebyshev quadrature presents explicit equations for the calcuation of points and weights, the same is not true for Gauss-Laguerre. The calculation of the points is computationally heavy, as it requires variable precision. For N > 40, the range of weights is around 10<sup>-60</sup>, such that standard double precision would not be able to calculate the collocation and integration points accurately.

For such reason, some common quadrature points were precalculated and stored in the folder **DATA_GLAG**. The function `[s,t,W] = GLAG_POINTS(N,cas)` can be used to simply read from the calculated quadrature points from the file.



[1]: https://ora.ox.ac.uk/objects/uuid:c4f497db-0a6d-4f8f-b9c0-558653fc97f8
[2]: https://www.springer.com/gp/book/9780792338482
[3]: https://www.sciencedirect.com/science/article/pii/0045794981900845
