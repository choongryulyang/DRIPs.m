## A Matlab Package for LQG Dynamic Rational Inattention Problems (DRIPs)

This package provides a fast and robust method for solving LQG Dynamic Ratinoal Inattention Problems based on methods developed in [Afrouzi and Yang (2019): Dynamic Rational Inattention and the Phillips Curve](http://www.afrouzi.com/dynamic_inattention.pdf). 

### Major updates
* Oct 07, 2019: first version of **Matlab solver** for the **steady state distribution** was added.
* Feb 22, 2020: first version of **Julia solver** for the **steady state distribution** was added.
* Apr 10, 2020: first version of **Julia solver** for the **transition dynamics** was added.
* Nov 18, 2020: first version **Matlab solver** for the **transition dynamics** was added. This corresponds to the v0.2.0 of the [DRIPs.jl](https://afrouzi.com/DRIPs.jl/dev/) package.

### Resources
* [Afrouzi and Yang (2019): Dynamic Rational Inattention and the Phillips Curve](http://www.afrouzi.com/dynamic_inattention.pdf)

### Codes
The following source files are included in the `src` folder:
* [`Drip.m`](src/Drip.m)
    This function defines a DRIP given inputs from the user. The function returns a structure that stores the primitives of the problem as well as the steady state information structure. Every DRIP should be defined through this function before transition dynamics can be solved. See the documentation within the `Drip.m` file for details on inputs and output objects.

* [`Trip.m`](src/Trip.m)
    This function characterizes the transition dynamics of a DRIP. Its inputs are (1) a `p` structure that is the output of `Drip.m`, and an initial prior covariance matrix that specifies the initial condition. See the documentation within the `Trip.m` file for details on inputs and output objects.

* [`irfs.m`](src/irfs.m)
    Given the information structure from a Drip or a Trip, this function returns the impulse response functions of actions and beliefs to structural shocks. See the documentation within the `irfs.m` file for details on inputs and output objects.

### Examples and Replications
* [`ex1_pricing_wo_feedback.m`](examples/ex1_pricing_wo_feedback.m)
    Solves a pricing problem with no strategic complementarity

* [`ex2_sims2010.m`](examples/ex2_sims2010.m)
    Solves the example from [Sims (2010)](http://sims.princeton.edu/yftp/RIMP/handbookChapterRI2.pdf)  and provides solution and irfs under transition dynamics based the extension in Afrouzi and Yang (2020).

### Deprecated Source Files
*  ***(Note: this is a deprecated version. This file was updated to Drip.m in the src folder.)*** 
    [Solve_RI_Dynamics.m](deprecated/Matlab/Solve_RI_Dynamics.m)
* ***(Note: this is a deprecated version. The Julia solver is now provided as a package. Please visit https://github.com/afrouzi/DRIPs.jl.)***
    [DRIP.jl](deprecated/Julia/DRIP.jl) 
