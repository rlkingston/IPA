# Version  1.0.0

+ Effective automation of the procedures described in Kingston and Millane (2022), IUCrJ, 9, 648-665.
+ Envelope and phase determination steps are parallelized.

# Version  1.1.0

Major changes

+ Envelope comparison and clustering now allows for inversion symmetry.
+ The way the apodization function is specified has changed.
+ Checks on the magnitude of missing terms are now based entirely on the Wilson model. 
+ Agreement of the data with the Wilson model is reported.
+ New procedures introduced for controlling the amplitudes of missing data during Fourier space projection onto the constraints  
+ A much more flexible method of defining the update rules and constraints was introduced. 
+ During DBSCAN clustering of envelopes, "k-distances" are now evaluated. Optionally use the k-distance distribution to set threshold value ε. 
+ Internal handling of Fourier data was reworked. Fourier coefficients beyond the current "effective resolution limit" are now removed, speeding up the calculations with heavily apodized data.
+ Ouput from analysis of overall scale and B-factor (method of Rogers) is now reported in the summary log file for the job.
+ Parameterization of the envelope determination step has changed, significantly improving algorithm performance for some test cases. 

# Version  1.2.0

Major changes

+ Code restructuring to allow user to define the major steps to be perfomed and how the inputs and outputs of those steps are handled.
+ Updated parameter file syntax, allowing better control of what happens at each iterate.
+ Better procedures for outlier detection and removal in the imput data set. 
+ Addition of RRR and RAAR algorithms and their reversed variants.
+ Implementation of procedures to add error to an existing  phase set using appropriate circular probability distributions.  
+ Changes to the way the apodzation function is applied.
+ Better procedures for handling random number generation, facilitating algorithm comparisons.
+ Algorithms optimized to increase speed of execution.
+ Added ability to average phases over the algorithm trajectory, following convergence to the solution, and calculate the phase retrival transfer function.
