// IPA: Iterative Projection Algorithms for protein crystallography
//  Header for the IPA manager class  

/*
 This class takes the program inputs, and passes them onto ipa_worker instances, resulting in the execution of individual runs of the IP algorithms.
 The number of threads is specified by the user, however this is handled in a safe manner so that the number of threads cannot exceed the hardware limit of the machine.
 Furthermore, during the phase determination step, this manager repeatedly clusters and attempts to form a consensus from the phase sets produced by the IP algorithms. 
 If highly consistent phase sets are detected,  this almost certainly means phase retrival has been succcessful, and the user is alerted.
 This allows the user to manually inspect the output maps, and terminate the program if desired. 
 */

#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>

#include <thread>
#include <atomic>
#include <chrono>

#include "ipa_worker.h"
#include "util-lib.h"
#include "svn-version.h"
#include "consensus_worker.h"

#ifndef IPA_MANAGER
#define IPA_MANAGER

class ipa_manager {
  
public:
  // ------------------------------------------------------- //
  // ---------------- PUBLIC Methods ----------------------- //
  // ------------------------------------------------------- //
  
  ipa_manager();
  
  int perform_ipa_experiment(const ipa_settings RunSettings, const int n_threads);
  
  void safely_update_thread_count(const int n_threads);
  
  // ------------------------------------------------------- //
  // ---------------- PUBLIC Variables --------------------- //
  // ------------------------------------------------------- //
  
  bool quiet_mode = false; // Set this to true if we don't want to print the update screens, for debug as updating screen can hide obscure console messages.
  
  std::atomic<int> iterations_performed;
  
  std::string current_job_name = "Undefined Job"; // Might want to rename this as "jobname" is confusing in context with exp_managers jobname which dictates folder structure.
  
  std::string current_user_job_name = "Undefined Jobname"; // Might want to rename this as "jobname" is confusing in context with exp_managers jobname which dictates folder structure.
  
  consensus_phase_settings my_consensus_phase_settings;   // Set these consensus phase settings when phase determination is happening.

  std::vector<cluster_data> consensus_solutions;
  bool consensus_solutions_is_valid = false;

  std::atomic<int> n_solutions_found;
  std::atomic<bool> solutions_are_found;

  bool infinity_mode = false; // Turn to true during phase determination to endlessly perform runs until a consensus is found.
  
private:
  // ------------------------------------------------------- //
  // ----------------------- Methods ----------------------- //
  // ------------------------------------------------------- //
  
  // IPA threaded job:
  int run_ipa_job(std::atomic<bool> &FinishTheThread,
                  std::atomic<float> &MyProgress,
                  std::atomic<float>& CUSUM,
                  int MyThreadID,
                  int run_id,
                  ipa_settings JobSettings);
  
  // Phase Consensus Worker:
  void initialise_consensus_worker();
  bool run_phase_consensus_threaded(std::atomic<bool>& FinishTheThread,
                                    std::atomic<float>& MyProgress);
  
  // Update screen and associated helper functions:
  void reprint_update_screen(std::atomic<float> ThreadProgress[],
                             std::atomic<bool> ThreadStatus[], 
                             std::atomic<float> CUSUMS[]);
  
  // ------------------------------------------------------- //
  // ----------------------- Variables --------------------- //
  // ------------------------------------------------------- //
  
  float total_progress; // float between 0 and 1, indicating the total progress of this step.
  
  int update_rate = 1000; // How often the console will update during a run in milliseconds, default 1000 ms.
  
  int total_runs = -1; // The total number of runs, is updated by settings, but can be updated during infinity mode.
  
  int max_threads; // Safely updated via settings.
  
  int iterations_since_last_consensus = 0;
  
  int iterations_total; // Total iterations to be completed.
  
  std::atomic<int> runs_complete; // Number of IPA iterations completed.
  
  bool performing_consensus_phasing = false;
    
  time_t start_time;
  
  std::vector<std::string> threading_comms;
  
  consensus_worker my_consensus_phase_worker;
};

#endif
