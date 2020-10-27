## A Software Package for LQG Dynamic Rational Inattention Problems (DRIPs)
**Major updates:** 
* Oct 07, 2019: first version of **Matlab solver** for the **steady state distribution** was added.
* Feb 22, 2020: first version of **Julia solver** for the **steady state distribution** was added.
* Apr 10, 2020: first version of **Julia solver** for the **transition dynamics** was added.
### Resources
* [Afrouzi and Yang (2019): Dynamic Rational Inattention and the Phillips Curve](http://www.afrouzi.com/dynamic_inattention.pdf)
### Codes
* [Solve_RI_Dynamics.m](Matlab/Solve_RI_Dynamics.m) - A Matlab solver for **a dynamic multidimensional rational inattention** problem
* [DRIP.jl](Julia/DRIP.jl) - (Note: this is a **Deprecated** version. The Julia solver is now provided as a package. Please visit https://github.com/afrouzi/DRIPs.jl.)
  A Julia solver for **a dynamic multidimensional rational inattention** problem. Contains:
  1.  `solve_drip`: a function for solving the **steady state distribution** of actions and shocks.
  2.  `solve_trip`: a function for solving the **transition dynamics** to the steady state distribution from an aribitrary initial distribution.
  3. `dripirfs`: a function for generating the impulse response functions of actions and beliefs under the **steady state rational inattention distribution**.
### Jupyter notebooks
An interactive set of examples using Jupyter notebooks and mybinder are available at [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/choongryulyang/dynamic_multivariate_RI/master) (no software is needed on the local machine).

