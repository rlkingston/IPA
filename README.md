# IPA 
***I**terative **P**rojection **A**lgorithms for protein crystallography*

A program for *ab initio* phase determination and phase improvement.


# Overview

+ The program is capable of *ab initio* phase determination for protein crystals with high solvent content (solvent fraction > 70%) 

+ It can also be used to generate plausible molecular envelopes for crystals with much lower solvent content.

The only required input are the measured diffraction data (in a CCP4 mtz file) and an estimate for the solvent fraction of the crystal.


# The Current Distribution

# Running the program

The behavior of program IPA is controlled by a single experimental parameter file. Running the program without a parameter file will generate a default file. The user can then edit the default file. The only compulsory inputs are (1) Some details of the mtz file that contains the diffraction data, and (2) an estimate of the solvent fraction of the crystal. The latter is generally obtained through analysis of crystal packing density (e.g. with the CCP4 program MATTHEWS_COEFF). If there is uncertainty about the solvent fraction, it is better to **underestimate** this quantity. If the solvent fraction is overestimated, *ab initio* phasing will certainly fail.


# Contact

Address correspondence to: r.l.kingston@auckland.ac.nz
Report issues to michael.barnett@auckland.ac.nz or r.l.kingston@auckland.ac.nz
