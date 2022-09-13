// IPA: Iterative Projection Algorithms for protein crystallography
// IPA manager class  

#include "ipa_manager.h"

ipa_manager::ipa_manager()
{
  /*
   IPA Manager Constructor
   */
  
  // Initialise the threading communication ids, encoded in the Integer of ThreadProgress float.
  //(whereas the total progress is encoded in the fractional component)
  threading_comms = {
    "Undefined",                              // 0
    "Idle",                                   // 1
    "Initialising",                           // 2
    "Apodization",                            // 3
    "IPA: Difference Map",                    // 4
    "IPA: Relax Reflect Relax",               // 5
    "IPA: Error Reduction",                   // 6
    "IPA: Reverse Relax Reflect Relax",       // 7
    "Calculating Statistics",                 // 8
    "Checking for Consistent Phase Sets"      // 9
  };
}

int ipa_manager::perform_ipa_experiment(const ipa_settings run_settings,const int n_threads)
{
  /*
   This function initialises and executes the operations for an entire experiment. 
   It is the primary point of control for threading, as the atomic variables which control the threads can not be member variables, and will go out of scope without careful parsing.
   The function takes the settings for an IPA experiment, and will pass the settings onto the instanced worker objects.
   */
  
  // record the start time...
  time(&start_time);
  iterations_performed = 0;
  
  solutions_are_found = false;
  n_solutions_found = 0;
  
  total_runs = run_settings.n_runs;
  runs_complete = 0; // Reset this.
  
  if (run_settings.job_is_phase_determination && run_settings.check_for_consensus_phases)     // if we are phasing, we should init the consensus worker here upon request
  {
    initialise_consensus_worker();
    std::cout << "Phase Determination: Consensus worker initialised. \nA check for consistent solutions will be made at the end of each threaded run cycle." << std::endl;
    performing_consensus_phasing = true; // member variable to store this info.
  }
  
  //std::cout << "Total runs requested: " << total_runs << std::endl;
  iterations_total = run_settings.n_total_iterates * total_runs; // n_total_iterates is calculated at the parameter parsing steps. Assumes all runs are identical.
  
  //std::cout << "Total iterates per run: " << run_settings.n_total_iterates << std::endl;
  
  //std::cout << "Splitting job across " << n_threads << " threads" << std::endl;
  
  /*
  // Find out how many iterates per run are required across all runs.
  // This currently assumes they are all identical runs...
  ipa_worker calcworker; // Calc worker is a little hacky but oh well.
  calcworker.settings = run_settings; // Give the calc worker the settings.
  calcworker.calculate_total_number_of_iterates();
  
  if (run_settings.n_total_iterates == 0)
  { // only do old way if its 0 i.e. not set...
    // This will become redundant once new way is fully tested.
    iterations_total = calcworker.n_total_iterates * total_runs;
  }
  */
  
  // confirm the number of threads being used for this job.
  safely_update_thread_count(n_threads);
  
  // Declare (and fill) important threading objects;
  std::atomic<bool> finished_thread_flags[max_threads];
  std::atomic<float> ThreadProgress[max_threads];
  std::atomic<float> CUSUM_Values[max_threads];
  
  for (int i = 0; i < max_threads; i++)
  {
    finished_thread_flags[i] = true;
    ThreadProgress[i] = 0.0;
    CUSUM_Values[i] = 0;
  }
  
  std::thread Threads[max_threads];
  
  std::cout << "Initialised " << max_threads << " threads for " << total_runs << " jobs.\n" << std::endl;
  std::cout << std::endl;
  
  // Clear some space for the update screen, so previous logs are not lost:
  // 75 lines is 'about right', but this can be refined.
  for (int i = 0; i < 75; i ++)
  {
    std::cout << " " << std::endl;
  }
  
  // No std::couts within this loop!
  for (int run = 0; run < total_runs; run++)
  {
    /*
     We start by identifying if there is an available free thread, during this time we are also updating the summary screen for the user's convenience to monitor progress.
     */
    
    // Increment total runs if we haven't found a solution, and infinity mode is turned on.
    // safety check for phasing is also in place - this shouldnt ever happen under envelope consensus.
    if ((run == total_runs-1) && infinity_mode && !solutions_are_found && performing_consensus_phasing)
    {
      total_runs++;
    }
    
    int next_thread_id = -1;
    while (next_thread_id < 0)
    {
      for (int id = 0; id < max_threads; id++)
      {
        if (finished_thread_flags[id] == true)
        {
          next_thread_id = id;
          break;
        }
        
      }
      
      // Update the screen, and don't check this more often than 1000 ms.
      total_progress = float(iterations_performed) / float(iterations_total);
      reprint_update_screen(ThreadProgress,finished_thread_flags, CUSUM_Values);
      
      std::this_thread::sleep_for(std::chrono::milliseconds(update_rate)); // Sleep before updating console.
    }
    
    // Ensure we join() this thread before replacing with new one.
    if (Threads[next_thread_id].joinable()) Threads[next_thread_id].join();
    
    // Check whether we need to run a phase consensus step.
    // Note using the modulus of max threads for the run number allows for the phasing consensus to wait for the nth thread to finish before running a consensus (as all threads are likely to finish within seconds of each other.)
    
    // We only try to find a consensus after an entire threading run is done (as most threads finish within seconds of each other)
    // Thus to do a run, the consensus worker must be idle, it there must be at least 2 runs completed already, and the last finished run must be divisible by the number of threads.
    if (my_consensus_phase_worker.is_idle && (runs_complete%max_threads == 0) && (run > 1) && (iterations_since_last_consensus < runs_complete) && performing_consensus_phasing)
    {
      // A Batch of threaded runs just finished, lets process them!
      // Add all the files from the list.
      std::ifstream input_file(my_consensus_phase_settings.phase_list_file);
      std::vector<clipper::String> file_list;
      clipper::String line;
      
      while(std::getline(input_file, line))
      {
        my_consensus_phase_worker.add_phase_input_to_list(line);
      }
      
      input_file.close();
      
      iterations_since_last_consensus = runs_complete;
      finished_thread_flags[next_thread_id] = false;
      Threads[next_thread_id] = std::thread(&ipa_manager::run_phase_consensus_threaded,
                                            this,
                                            std::ref(finished_thread_flags[next_thread_id]),
                                            std::ref(ThreadProgress[next_thread_id]));
      // We need to decrement run so that we don't lose a run every time consensus happens...
      run--;
    }
    else // No consensus to perform, start the next IPA run thread.
    {
      finished_thread_flags[next_thread_id] = false;
      Threads[next_thread_id] = std::thread(&ipa_manager::run_ipa_job,
                                            this,
                                            std::ref(finished_thread_flags[next_thread_id]),
                                            std::ref(ThreadProgress[next_thread_id]),
                                            std::ref(CUSUM_Values[next_thread_id]),
                                            next_thread_id,
                                            run,
                                            run_settings);
    }
  }
  
  // All runs have been assigned ... Wait for threads to all finished....
  bool check = true;
  while (check)
  {
    check = false;
    for (int id = 0; id < max_threads; id++)
    {
      if (finished_thread_flags[id] == false)
      {
        check = true;
      }
    }
    
    // Everything has now finished ... in the case of phase determination, run the consensus worker a final time if requested 
    if (check == false && (iterations_since_last_consensus < runs_complete) && performing_consensus_phasing)
    {
      
      // Ensure we join() this thread before using it.
      if (Threads[0].joinable()) {Threads[0].join();}
      
      // Do one last phase consensus now that everything is finished.
      //std::cout << "Final check for consistent solutions during the phase determination step " << std::endl;
      // Add all the files from the list.
      //std::cout << "Opening file... " << std::endl;
      std::ifstream input_file(my_consensus_phase_settings.phase_list_file);
      std::string line;
      while(std::getline(input_file, line)) {
        my_consensus_phase_worker.add_phase_input_to_list(line);
        //file_list.push_back(line);
      }
      //std::cout << "Added files to worker list." << std::endl;

      
      iterations_since_last_consensus = runs_complete;
      finished_thread_flags[0] = false;
      Threads[0] = std::thread(&ipa_manager::run_phase_consensus_threaded,
                                            this,
                                            std::ref(finished_thread_flags[0]),
                                            std::ref(ThreadProgress[0]));
      check = true;
    }
    
    // Update the screen, and don't check this more often than 100 ms.
    total_progress = float(iterations_performed) / float(iterations_total);
    reprint_update_screen(ThreadProgress,finished_thread_flags, CUSUM_Values);
    std::this_thread::sleep_for(std::chrono::milliseconds(update_rate));
  }
  
  total_progress = 1;
  
  std::cout << "Rejoining all threads." << std::endl;
  // Join all threads back together...
  for (int i = 0; i < max_threads ; i++) {
    if (Threads[i].joinable()) Threads[i].join();
  }
  
  std::cout << "Threading is complete!\n" << std::endl;
  
  return 1;
}

int ipa_manager::run_ipa_job(std::atomic<bool> &FinishTheThread, std::atomic<float> &MyProgress, std::atomic<float>& CUSUM, int MyThreadID, int run_id, ipa_settings JobSettings)
{
  /*
   This function is responsible for generating a worker instance, which will perform a single run (The function is called once per thread).
   Importantly, only variables that are flagged atomic are parsed to the worker instances, and the settings are copied in memory for thread safety.
   
   The worker updates the ThreadProgress with a number dictating its current state
   When the run concludes the bool FinishTheThread is set to true
   */
  
  int output;
  
  ipa_worker MyIpaWorker(JobSettings);
  //MyIpaWorker.settings = JobSettings; // Here we set the settings to what we want.
  //MyIpaWorker.look_for_own_settings = false;

  int run_id_offset = run_id + JobSettings.first_run_id_offset;
  MyIpaWorker.run_experiment(MyProgress, std::ref(iterations_performed), CUSUM, run_id_offset);
  
  FinishTheThread = true;
  MyProgress = 1;
  runs_complete++;
  return output;
}

void ipa_manager::initialise_consensus_worker()
{
  // Helper function to initialise the consensus worker for phase determination in desired state.
  my_consensus_phase_worker.settings = my_consensus_phase_settings;
  my_consensus_phase_worker.is_idle = true;
}

bool ipa_manager::run_phase_consensus_threaded(std::atomic<bool>& FinishTheThread, std::atomic<float>& MyProgress)
{
  /*
   Here we run the consensus_worker to check if the current list of solutions contains any which are very similar, in which case we average the results to produce a consensus
   This is typically run in a threaded manner, importantly additional solutions cannot be added whilst it is running.
   Finished phase files are added to the phase_list file.
   I think we will want to run this after every nth run, where n is the number of threads, as the runs will tend to finish in batches of how many threads are being used.
   */
  
  MyProgress = 9.0;
  bool success = my_consensus_phase_worker.check_for_consensus_phases(MyProgress);
  consensus_solutions.clear();
  
  // Maybe not with the outputs to std::cout here, seeming as we're all threaded.
  FinishTheThread = true;

  my_consensus_phase_worker.return_list_of_consensus_solutions(consensus_solutions, std::cout);
  consensus_solutions_is_valid = true;
  
  n_solutions_found = consensus_solutions.size();
  FinishTheThread = true;
  MyProgress = 1;
  return success;
}

void ipa_manager::safely_update_thread_count(int n_threads)
{
  /*
   Simple function to control the private max_threads variable. Although setting the number of available threads to be greater than the hardware number of threads/cores would not be catastrophic, it would slow down each thread as multiply threads would run on same cores.
   */
  
  // parse the thread count and give warnings if necessary...
  int n_hardware_thread_count = std::thread::hardware_concurrency();
  if (n_threads > n_hardware_thread_count) {
    std::cerr << "You have requested more threads (" << n_threads << ") than are supported by hardware (" << n_hardware_thread_count << ")." << std::endl;
    n_threads = n_hardware_thread_count;
  }
  
  if (n_threads < 1) n_threads = 1; // Sanity check for bad values.
  
  max_threads = n_threads; // Set the number of threads we're using.
}

void ipa_manager::reprint_update_screen(std::atomic<float> ThreadProgress[], std::atomic<bool> ThreadStatus[], std::atomic<float> CUSUMS[])
{
  if (quiet_mode)
  {
    return;
  }
  /*
   This function is called whenever a visual update of the progress of the program is desired. Intention is to pass floats that represent the total progress and individual threaded processes, so these can be visually supplied to the user as loading bars.
   
   Other feedback individual to each thread (like current outputs) may also be passed in the future to fancy it up a bit. i.e. what the job is currently doing etc...
   
   Only works on OS or consoles which allow for ANSI escape characters... which is probably most of them.
   
   Importantly, whilst this function will over-write itself, other outputs to std::cout will slowly accumulate, taking up memory. It is best to use custom log files for all outputs now.
   */
  
  std::cout << "\033c"; // Clear the console using an ANSI escape hack.
  std::cout << "\033[H"; // Different kind of clear...
  
  int cl = 78; // Custom length (width of window) of output. Currently hard coded to 80, as that is a default console size across many platforms. However there are ways to dynamically read this and adjust on the fly which are OS specific. Or we can change it to whatever we like best, and people can change the size of their console to match.
  
  // Get the timing formation
  time_t current_time;
  time(&current_time);
  int seconds_total, seconds_estimated;
  seconds_total = int(difftime(current_time,start_time));
  
  // calculate the estimated time left by wildly extrapolating.
  seconds_estimated = int(float(seconds_total) / total_progress) - seconds_total;
  
  int s = seconds_total%60;
  int m = (seconds_total/60)%60;
  int h = (seconds_total/3600)%24;
  int d = (seconds_total/86400);
  
  int se = seconds_estimated%60;
  int me = (seconds_estimated/60)%60;
  int he = (seconds_estimated/3600)%24;
  int de = (seconds_estimated/86400);
  
  //current_time = localtime(&timeobj);
  //std::chrono::time_point<std::chrono::> duration = current_time - start_time;
  
  // Header Information
  util::print_section_line(cl, std::cout);
  util::print_empty_line(cl, std::cout);
  std::cout << util::parse_line(cl,"Program IPA: Iterative Projection Algorithms for protein crystallography v " + std::string(xmake_string(release_version)));
  util::print_empty_line(cl, std::cout);
  std::cout << util::parse_line(cl,"Job Directory: /" + current_user_job_name + "/");
  std::cout << util::parse_line(cl,current_job_name);
  util::print_empty_line(cl, std::cout);
  util::print_section_line(cl, std::cout);
  std::cout << "\n";
  
  // Phase determination information: Simple check to alert user that a consensus was found...
  if (performing_consensus_phasing)
    {
      char buf[256];
      util::print_section_line(cl, std::cout);
      std::cout << util::parse_line(cl,"Identification of highly similar phase sets using clustering: ");
      if (n_solutions_found > 1)
      {
        int n = n_solutions_found; // Not atomic.
        std::cout << util::parse_line(cl,std::to_string(n) + " clusters have been identified!");
        std::cout << util::parse_line(cl,"Check /phase_consensus/" + current_user_job_name + "/ folder for output .mtz and .ccp4 files!");
      }
      else if (n_solutions_found == 1)
      {
          std::cout << util::parse_line(cl,"A cluster has been identified!");
          std::cout << util::parse_line(cl,"Check /phase_consensus/" + current_user_job_name + "/ folder for output .mtz and .ccp4 files!");
      }
      else
      {
        std::cout << util::parse_line(cl,"No clusters have been identified. Phase determination has not yet succeeded");
        util::print_empty_line(cl, std::cout);
      }
      
      sprintf(buf,"Clustering was last performed after run %i: ", iterations_since_last_consensus);
      std::cout << util::parse_line(cl,std::string(buf));

      util::print_section_line(cl, std::cout);
      std::cout << "\n";
    }
  // Total Progress Information
  util::print_section_line(cl, std::cout);
  util::print_empty_line(cl, std::cout);
  char buf[256];
  
  int iii = iterations_performed;
  int runs_completed = runs_complete; // So that Atomic<int> plays nicely.
  sprintf(buf,"Total Progress:   %i/%i     Runs Complete: %i/%i",iii, iterations_total, runs_completed, total_runs);
  std::cout << util::parse_line(cl,std::string(buf));
  // Loading bar
  util::print_progress_bar(cl, total_progress, std::cout);
  
  // Time remaining estimates
  if (de > 0 || d > 0) {
    sprintf(buf,"Time taken (remaining):       %02i days, %02i:%02i:%02i (%02i days, %02i:%02i:%02i)",d ,h, m, s, de, he, me, se);
  } else {
    sprintf(buf,"Time taken (remaining):       %02i:%02i:%02i (%02i:%02i:%02i)",h, m, s, he, me, se);
  } // Only print days when necessary.
  
  std::cout << util::parse_line(cl,std::string(buf));
  util::print_empty_line(cl, std::cout);
  util::print_section_line(cl, std::cout);
  std::cout << "\n";
  
  // Threading Status Information
  util::print_section_line(cl, std::cout);
  util::print_empty_line(cl, std::cout);
  
  for (int th = 0; th < max_threads; th ++)
  {
    // Decode the fraction and non fractional parts of the progress float
    float nonfrac;
    float frac = std::modf(ThreadProgress[th],&nonfrac);
    int myCUSUM = CUSUMS[th];
    
    if (nonfrac > 9) nonfrac = 0; // Out of bounds safety check.
    
    // Prepare status string
    std::string status;
    if (ThreadStatus[th] == true)
    {status = "Idle";} else {status = threading_comms[nonfrac];}
    
    // Prepare final string
    char buffer[512];
    sprintf(buffer,"Thread[%i] CUSUM = %i :",th, myCUSUM);
    
    // Parse into line
    std::cout << util::parse_line(cl,std::string(buffer) + status);
    
    // Include progress bar
    util::print_progress_bar(cl, frac, std::cout);
  }
  
  // Tail
  util::print_empty_line(cl, std::cout);
  util::print_section_line(cl, std::cout);
  
  std::cout << std::endl;
}




