# DG1D
An object-oriented, 1-D discontinuous Galerkin finite element model for hyperbolic PDEs in MATLAB.

The DG discretization object is initialized as `DG(p, Mesh, ProbDef)`, where `p` is the DG polynomial degree, `Mesh` is a `MeshClass` object, and `ProbDef` is a structure containing the problem definition. See examples in the `\examples` directory.

Explicit Runge--Kutta, linear multistep, and multistep-multistage (multistep Runge--Kutta) methods are included with the code, along with appropriate CFL coefficients for use with the DG spatial discretization.

Requires [PolyTools](https://www.mathworks.com/matlabcentral/fileexchange/72256-polytools) and [QuadTools](https://www.mathworks.com/matlabcentral/fileexchange/72378-quadtools) toolboxes for MATLAB.
