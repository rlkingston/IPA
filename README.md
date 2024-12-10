# IPA 
***I**terative **P**rojection **A**lgorithms for protein crystallography*

A program for **direct** (*ab initio*) phase determination.


# Overview

+ Can perform **direct** phase determination for protein crystals with high solvent content (generally solvent content > 70%)

+ Can be used to generate plausible molecular envelopes for crystals with much lower solvent content.

+ A useful platform for exploring the use of Iterative Projection Algorithms for crystallographic phase retrival.

Fundamental to the approach is the treatment of *ab initio* phase determination as a constraint satisfaction problem. An image is sought that is consistent with both the diffraction data and generic constraints on the density distribution in the crystal (like the solvent flatness constraint). The problem is solved using iterative projection algorithms which have good global convergence properties, and can locate the solution without any initial phase information, so long as the solution is uniquely specified by the constraints. This approach can be viewed as an evolution of the density modification procedures traditionally used in protein crystallography to improve experimentally determined phases. 

The only required inputs are the measured diffraction data (in a CCP4 mtz file), an estimate for the solvent fraction of the crystal, and some information about the atomic composition of the asymmetric unit.


# Compiling from source

The source code can be found in the zip archive associated with the latest release.  It should be readily compilable under Mac OS and Linux.  The only major dependency is the CCP4 package (or more specifically the Clipper crystallographic library). See the detailed instructions in the src directory (Install.txt). Data and a parameter file to allow solution of a test case (PDB ID 4fzn) are included in sub-directory test_case. Some bare bones documentation about the paramter file syntax is included in sub-directory doc


# Running the program

The behaviour of program IPA is controlled by a single experimental parameter file. Running the program without a parameter file will cause the program to generate a default parameter file, which the user can edit. The default parameter file will execute the two-stage envelope and phase determination procedure described in Kingston and Millane (2022). 

The compulsory edits to the default parameter file are:

(1) Adding details of the mtz file that contains the diffraction data.

The program expects estimates of the Structure Factor Amplitudes |F| and their standard deviations. Standard Bayesean treatment of weak data  is anticipated. If the data have been corrected for significant anisotropy (e.g. as is optionally done via the STARANISO server), this will cause the estimation of the overall anisotropic B-factor within program IPA to perform unpredictably, with some bad downstream consequences ...

(2) Providing an estimate of the solvent fraction of the crystal. 

The solvent fraction is generally estimated via analysis of crystal packing density (e.g. with the CCP4 program MATTHEWS_COEFF). If there is uncertainty about the solvent fraction, it is better to **underestimate** this quantity. If the solvent fraction is overestimated, *ab initio* phase determination will certainly fail (setting the solvent fraction too high is inconsistent with the solution, while setting it too low remains consistent with the solution).

(3) Providing some detail of the atomic content of the asymmetric unit (by specifying the protein sequence)

This information is used to put the Structure Factor amplitudes on an approximately absolute scale via analysis of the Patterson origin peak (using a procedure devised by Donald Rogers and further developed by Bob Blessing and David Langs). If your data is already on an absolute scale you can actually omit specifying this information, but we don't recommend it. 


For computational efficiency, the problem of **direct** phase determination is broken into two stages; initial approximation of the molecular envelope at low resolution, followed by subsequent phase determination using all of the data. At both stages, the algorithm is initiated with many different and random phase sets, which are evolved subject the constraints. Averaging across the final part of the algorithm trajectory is now performed by default, which improves the image if the algorithm has converged to the solution.  A clustering procedure is used to identify consistent results across multiple runs, which are themselves averaged to generate consensus envelopes or phase sets. **The emergence of highly consistent phase sets in the second step of the procedure is diagnostic of success**.

You can find the log and summary files for the individual runs of the program, as well as the logs for the comparison, clustering and averaging operations in:

./job_name/envelope_determination/logs

./job_name/phase_retrieval/logs
 
The outputs of the  comparison, clustering and averaging operations can be found in:

./job_name/envelope_determination/consensus

./job_name/phase_retrieval/logs/consensus


In particular, the program will output the proposed solution (both Fourier coefficients and maps) to the directory phase_retrieval/consensus, for inspection and interpretation. As with all *ab initio* phase determination methods, the procedure may generate the true solution or its inverse. The correct hand of the reconstructed image must be determined by inspection. 


The above describes what will happen when executing the program using the default parameter file. However the program IPA allows the flexible use of many different iterative projection algorithm, in combination with different real space constraints on the crystallographic image. The user is not restricted to the supplied protocol, you can invent your own. For details on how the program is controlled via the parameter file, see the bare bones documentation. 


# References

Version 1.1 of the program effectively automates the procedures described in:

Kingston, R.L. and Millane R.P. (2022) **A general method for directly phasing diffraction data from high–solvent content protein crystals.** IUCrJ, 9, https://doi.org/10.1107/S2052252522006996

Version 1.2 of the program  has many additional features, some of which are described in

Barnett, M.J., Millane R.P. and Kingston, R.L. **Analysis of crystallographic phase retrieval using iterative projection algorithms** Acta Cryst D, 80, https://doi.org/10.1107/S2059798324009902



# Contact

Address correspondence to: r.l.kingston@auckland.ac.nz

Report issues to r.l.kingston@auckland.ac.nz
