// IPA: Iterative Projection Algorithms for protein crystallography

// Contributors: Richard Kingston, Michael Barnett

// Verson History
// 1.0.0 June 2022. First public release - Automates the procedures described in Kingston & Millane. IUCrJ
// 1.1.0 Sept 2022

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <random>
#include <thread>
#include <atomic>
#include <chrono>
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

#include "ipa-lib.h"
#include "util-lib.h"
#include "ipa_worker.h"
#include "consensus_worker.h"
#include "ipa_manager.h"
#include "exp_manager.h"
#include "svn-version.h"

#ifndef IPA_PROGRAM
#define IPA_PROGRAM

int main( int argc, char** argv)
{
  // Start counting.
  clock_t iterate_start_time = clock();

  // Print version to the log.
  util::version_fingerprint(std::cout);

  // Echo the command line input when we execute CCP4CommandInput ... it will get sent to std:cout
  std::cout << "\nCommand Line Input:\n";
  bool echo = true;
  CCP4CommandInput args( argc, argv, echo );

  // Start off by checking if a custom experiment parameter filename was provided...
  // (use the default name otherwise.)
  std::string experimentparams = "experiment.params"; // default name
  if (argc == 2)
  {
    // A custom name for the experiment parameter file has been provided...
    experimentparams = std::string(args[1]);
    std::cout << "Will use experiment parameter file: \"" << experimentparams << "\"" << std::endl;
  }
  
  // Run the experiment using the named experiment parameter file.
  exp_manager my_experiment;
  bool success = my_experiment.run(experimentparams);
  
  if (success)
  { // On completion, report some information about the run time to std::cout.
    int n_seconds_total = (clock() -  iterate_start_time)/CLOCKS_PER_SEC;
    int n_seconds = n_seconds_total%60;
    int n_minutes = (n_seconds_total%3600)/60;
    int n_hours = n_seconds_total/3600;
    std::cout << "Total time taken: [" << n_hours << ":" << n_minutes << ":" << n_seconds << "] Hour:Min:Sec. \n" << std::endl;
    std::cout << "Program IPA ran successfully" << std::endl;
    return 0;
  }

  return 1;
}

#endif