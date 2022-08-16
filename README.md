# IPA 
***I**terative **P**rojection **A**lgorithms for protein crystallography*

A program for *ab initio* phase determination.


# Overview

+ The program is capable of *ab initio* phase determination for protein crystals with high solvent content (solvent fraction > 70%) 

+ It can also be used to generate plausible molecular envelopes for crystals with much lower solvent content.

*Ab initio* phase determination is treated as a constraint satisfaction problem, in which an image is sought that is consistent with both the diffraction data and generic constraints on the density distribution in the crystal. The problem is solved using  iterative projection algorithms which have good global convergence properties, and can locate the correct solution without any initial phase information. 

The only required inputs are the measured diffraction data (in a CCP4 mtz file) and an estimate for the solvent fraction of the crystal.


# The Current Release

The current release (v1.0.0) contains just the source code, and program outputs for a successful *ab intio* phase determination experiment on a test case (PDBID 2W4M, solvent fraction 0.77). The source code is for purposes of inspection and verification, and can't be directly compiled as distributed. 

Mac OS binary executables will be shared very shortly - in the mean time please contact us if you want to make use of the program. Contact details are below. 

# Running the program

The behavior of program IPA is controlled by a single experimental parameter file. Running the program without a parameter file will cause the program to generate a default parameter file, which the user can edit. 

The only compulsory edits are:

(1) Adding details of the mtz file that contains the diffraction data.

(2) Providing an estimate of the solvent fraction of the crystal. 

The solvent fraction is generally estimated via analysis of crystal packing density (e.g. with the CCP4 program MATTHEWS_COEFF). If there is uncertainty about the solvent fraction, it is better to **underestimate** this quantity. If the solvent fraction is overestimated, *ab initio* phase determination will certainly fail (setting the solvent fraction too high is inconsistent with the solution, while setting it too low remains consistent with the solution).

For computational efficiency, the program breaks the problem of *ab initio* phase determination into two stages; initial approximation of the molecular envelope at low resolution, followed by subsequent phase determination using all of the data. At both stages, the algorithm is initiated with many different and random phase sets, which are evolved subject the constraints. A clustering procedure is used to identify consistent results across multiple runs, which are then averaged to generate consensus envelopes or phase sets. **The emergence of highly consistent phase sets is diagnostic of success**.

The results are stored in 4 subdirectories:
+ envelope_determination
+ envelope_consensus
+ phase_determination
+ phase_consensus

If successful, the program will output the proposed solution (both Fourier coefficients and maps) in directory phase_consensus, for inspection and interpretation. As with all *ab initio* phase determination methods, the procedure may generate the true solution or its inverse. The correct hand of the reconstructed image must be determined by inspection. 

# Reference

Version 1.0 of the program effectively automates the procedures described in:

Kingston, R.L. and Millane R.P. (2022) **A general method for directly phasing diffraction data from high–solvent content protein crystals.** IUCrJ, In Press 

# Contact

Address correspondence to: r.l.kingston@auckland.ac.nz

Report issues to michael.barnett@auckland.ac.nz or r.l.kingston@auckland.ac.nz
