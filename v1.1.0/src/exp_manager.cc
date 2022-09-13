// IPA: Iterative Projection Algorithms for protein crystallography
// Experiment manager class

#include "exp_manager.h"

// TODO ... develop a sinple log file  for the experiment manager, which gives a condensed summary of the entire experiment and its results
// TODO ... ensure consistent handling of all Boolean settings in the experiment parameter file - in some cases true/false input is allowed and in others, the appearance of a keyword simply defaults the setting to true. 
 
// ------------------------------------------------------- //
// ----------------------- FORTRAN ----------------------- //
// ------------------------------------------------------- //

// Fortran subroutines required for run_initialization, from module atom_count

extern "C"  { void  count_protein_atoms(char *protein_sequence_filename, int *protein_sequence_nchains, int *atom_counts, char *output_filename); }
extern "C"  { void  count_na_atoms(char *dna_sequence_filename, bool * is_dna, int *dna_sequence_nchains, int *atom_counts, char *output_filename); }

// Fortran subroutines required for run_initialization, from module rogers_analysis

extern "C"  {void FTABLE(int *IZ, double *RM, double *A1, double *B1, double *A2, double *B2, double *A3, double *B3, double *A4, double *B4, double *C); }
extern "C"  {void METRIC(double *CELL, double *VCELL, double *G, double *GINV); }
extern "C"  {void PRELIM(int *n_data, int *index_h, int *index_k, int *index_l, double *normalized_intensity, double *CELL, double *VCELL, double *GINV, double *P_axial, double *PMIN, double *P0, double *PP, double *R); }
extern "C"  {void PLOT1(double *P_axial, char *output_filename); } 
extern "C"  {void GRID(double *RADIUS, double *GINV, double *U, double *V, double *W); }
extern "C"  {void PSUM(int *n_data, int *index_h, int *index_k, int *index_l, double *normalized_intensity, double *VCELL, double *U, double *V, double *W, double *Patterson); } 
extern "C"  {void PADJ(double *Patterson, double *PMIN, double *P0); }
extern "C"  {void PFIT(double *U, double *V, double *W, double *Patterson, double *PMIN, double *P0, double *PP, double *VCELL, double *F000, double *SUMZSQ, double *rogers_fraction, int *NCYCLE, double *R); }
extern "C"  {void PLOT2(double *U, double *V, double *W, double *Patterson, double *PMIN, double *P0, double *PP, double *RADIUS, double *G, char *output_filename); }
extern "C"  {void SOLVE(double *P0, double *PP, double *VCELL, double *G, double *GINV, char *crystal_family, double *SCALEK, double *B, double *U_CIF, double *U_STAR,  double *BISO, double *UISO, char *output_filename); }
extern "C"  {void PRINT_SOLUTION(double *RADIUS, int *NCYCLE, double *R, double *SCALEK,  double *B, double *U, double *BISO, double *UISO, char *output_filename); }

// ------------------------------------------------------- //
// ----------------------- METHODS ----------------------- //
// ------------------------------------------------------- //

bool exp_manager::run(std::string experimentparams)
{

  // Look for the experiment.params file
  // If no file is found, then generate a template and exit
  my_param_name = experimentparams;
  int settings_check = import_settings_from_file(envelope_settings, phasing_settings, experimentparams, std::cout);
  if (settings_check == 1)
  { // 1 = error code for could not find or open the parameter file:
    generate_default_experiment_params_file(experimentparams); // Generate a default file.
    std::cout << "Terminating run, please rerun with valid experiment parameter file." << std::endl;
    return false;
  }
  
  // Declare the working directories for each major folder based on the job name provided by user:
  dir_envelope_determination = "envelope_determination/" + job_name + "/";
  dir_envelope_consensus = "envelope_consensus/" + job_name + "/";
  dir_phase_determination = "phase_determination/" + job_name + "/";
  dir_phase_consensus = "phase_consensus/" + job_name + "/";
  dir_final_outputs = "final_outputs/" + job_name + "/";

  // We can open up a summary log file now too:
  summary_log_filename = job_name + "_summary_log.txt";
  summary_log.open(summary_log_filename, std::ostream::trunc); // We wipe previous log files of this filename.
  util::version_fingerprint(summary_log);

  // The validation of data and settings (TODO .. needs to be expanded, and program halted on egregious input violations.)
  found_input_file = validate_input_mtz();
  
  std::cout << "\u001b[31m"; // Following is coloured red for errors:
  validate_settings(envelope_settings);
  validate_settings(phasing_settings);
  std::cout << "\u001b[0m" << std::endl; // Reset colour code.
  
  // Print an output of the settings - typically for debugging - but also so the user can have one last check before hitting go.
  output_settings_summary(std::cout);
  std::cout << "Review the above settings, and press ENTER to continue..." << std::endl;
  std::cin.get();
  
  // Now we start and run each JOB in turn...
  output_settings_summary(summary_log); // Send the settings to the summary log now we have agreed to play.

  if (JOB_envelope_determination || JOB_phase_determination)
  {
    JOB_initialisation = true; // This is for now hard coded as INIT has become necessary. Except for a lone Phase or Envelope Consensus run ...
  }
  /*
  // TODO ... We don't need the anistropic B-factor analysis for the consensus jobs, when they are performed alone ... in which case the code below could control whether initialization is run or not  
  // However this would require that establishing the folder structure was split off from the rest of the initialization step
  if ( (JOB_phase_consensus ||  JOB_envelope_consensus) && !JOB_envelope_determination && !JOB_phase_determination)
    {JOB_initialisation = false;}
  else
    {JOB_initialisation = true;}
  */
    
  if (JOB_initialisation)
  {
    std::cout << "RUNNING JOB: INITIALISATION\n" << std::endl;
    run_initialisation();

    // Now we have the scale and overall B-Factor evaluated, and the Wilson model <I> calculated, we print out a comparison with the original data.
    util::print_intensity_statistics_data_vs_model(measured_f_all, mean_i_wilson, 50, summary_log); // hmmm why is this not working...

    util::print_intensity_statistics_data_vs_model(measured_f_all, mean_i_wilson, 50, std::cout);
  }
  
  if (JOB_envelope_determination)
  {
    // DO ENVE DETERMINATION
    std::cout << "RUNNING JOB: ENVELOPE DETERMINATION IN..." << std::endl;
    for (int i = 3; i > 0; i--)
    {
      std::cout << i << "..." << std::endl;
      std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }
    std::cout << std::endl;

    util::print_section_line(cl,summary_log);
    util::print_empty_line(cl, summary_log);
    summary_log << util::parse_line(cl, " ~ ENVELOPE DETERMINATION ~ ");
    util::print_empty_line(cl, summary_log);
    util::print_section_line(cl,summary_log);
    unpack_settings_as_log(summary_log, envelope_settings);
    
    // Perform Threaded Runs
    ipa_manager new_manager;
    new_manager.current_job_name = " ~ Envelope Determination ~ ";
    new_manager.current_user_job_name = job_name;
    new_manager.quiet_mode = no_display;
    new_manager.perform_ipa_experiment(envelope_settings, n_threads);
  }
  
  if (JOB_envelope_consensus)
  {

    util::print_section_line(cl,summary_log);
    util::print_empty_line(cl, summary_log);
    summary_log << util::parse_line(cl, " ~ ENVELOPE CONSENSUS ~ ");
    util::print_empty_line(cl, summary_log);
    util::print_section_line(cl,summary_log);

    // DO ENVELOPE MASK CONSENSUS
    std::cout << "RUNNING JOB: ENVELOPE CONCENSUS IN..." << std::endl;
    for (int i = 3; i > 0; i--)
    {
      std::cout << i << "..." << std::endl;
      std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }
    std::cout << std::endl;

    bool consensus_found = false;
    if (envelope_consensus_epsilon.size() == 0)
    {
      std::cout << " No epsilon values were provided for envelope consensus!\n The default value, " << envelope_consensus_mask_settings.epsilon << ", will be used." << std::endl;
      consensus_found = run_envelope_consensus();
    }

    for (int i = 0 ; i < envelope_consensus_epsilon.size(); i ++ )
    {
      std::cout << " Running consensus using epsilon value: " << envelope_consensus_epsilon[i] << std::endl;
      envelope_consensus_mask_settings.epsilon = envelope_consensus_epsilon[i]; // Set the current epsilon to use.
      consensus_found = run_envelope_consensus();
    }
    
    if (!consensus_found && JOB_phase_determination)
    {
      std::cout << "\nPhase determination was requested, however no consensus envelope emerged from the envelope determination step. If no initial envelopes were provided by the user, phase determination will be initiated without a starting envelope. The chance of success is small" << std::endl;
    }
  }
  
  // We now iterate over all candidate envelopes for the phase determination step, these are identified by the list_of_working_masks array.
  // This list is populated by a) the parameter file setting ... "-phasing mskin-target-working file.ccp4" 
  //                           b) the output envelopes from the envelope_consensus step.

  if (JOB_phase_determination)
  {

    // TODO ... Work out what to do if we are provided with a starting phase set (so scenario is "phase improvement" rather than ab initio phase determination)
    // If we initiate with the same starting phase set (what's enabled at present) the algorithm will behave deterministically
    // No matter how many times we run it, the outcome will always be the same. 
    // Best would probably be to add some random noise to the starting phase set, thus initiating slightly differently each time.
    // This could be a useful way to remove model bias if the starting phase set were derived from molecular replacement ? To be tested.
    // In this case we would probably bypass the envelope determination step 
    
    util::print_section_line(cl,summary_log);
    util::print_empty_line(cl, summary_log);
    summary_log << util::parse_line(cl, " ~ PHASE DETERMINATION ~ ");
    util::print_empty_line(cl, summary_log);
    util::print_section_line(cl,summary_log);
    unpack_settings_as_log(summary_log, phasing_settings);

    if (user_input_mask_count > 0)
      {
        std::cout << "\n" << user_input_mask_count << " envelopes were input via settings." << std::endl;
        summary_log << util::parse_line(cl, std::to_string(user_input_mask_count) + " envelopes were input via parameter file.");
      }
    
    if ( (list_of_working_masks.size() - user_input_mask_count) > 0)
      {
        std::cout << "\n" << list_of_working_masks.size() - user_input_mask_count << " consensus envelopes were generated by clustering." << std::endl;
        summary_log << util::parse_line(cl, std::to_string(list_of_working_masks.size() - user_input_mask_count) + " consensus envelopes were generated by clustering.");

        }
    
    if (list_of_working_masks.size() > 0)
    {
      std::cout << "\nEach envelope will now be used sequentially for Phase Determination." << std::endl;
      if (list_of_working_masks.size() > 1) infinity_mode = false; // We force infinity mode off here, otherwise we will never work our way through multiple envelopes.
    }
  
    //std::cout << "\n" << list_of_working_masks.size() << " consensus envelopes have been identified, or input via settings, each envelope will now be used sequentially for Phase Determination." << std::endl;

    // Deal with the special case, where there are no input working masks to start phase determination from...
    if (list_of_working_masks.size() == 0)
    {
      // There are no masks, however we can execute the phase determination runs without a starting mask...
      std::cout << " No starting envelopes are available, the Phase Determination step will be initiated without an envelope." << std::endl;
      summary_log << util::parse_line(cl, " No starting envelopes are available, the Phase Determination step was initiated without an envelope.");
      phasing_settings.n_iterations_fix_mask = 0; // Make sure this is 0, else problems.
      phasing_settings.input_working_mask = false; // Just in case.
      list_of_working_masks.push_back("NONE"); // Ensures the loop below runs exactly once, with a NONE mask.
    }
    
    // Initiate phase determination using all working masks (i.e. all inputs from the user, plus all outputs from the envelope concensus process)
    for (int mask_id = 0; mask_id < list_of_working_masks.size() ; mask_id ++)
    {
      // Before we start the phase determination run, we alter the offset values, so results get saved under different runs IDs. Thus no over-write of previous solutions.
      phasing_settings.input_filename_target_mask_working = list_of_working_masks[0]; // Set the phase settings input mask from list.
      int new_offset = (phasing_settings.n_runs * mask_id) + phasing_settings.first_run_id_offset ; // We create an offset value dependent on the mask_id
      std::cout << "Phase determination for envelope [" << mask_id <<  "]: " << list_of_working_masks[mask_id] << " will be associated with runs " << new_offset << " to " << new_offset + phasing_settings.n_runs << std::endl;
      phasing_settings.first_run_id_offset = new_offset; //
      
      // DO PHASE DETERMINATION
      std::cout << "RUNNING JOB: PHASE DETERMINATION IN..." << std::endl;
      for (int i = 3; i > 0; i--)
      {
        std::cout << i << "..." << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
      }
      std::cout << std::endl;
      
      // Perform a threaded run
      ipa_manager new_manager;
      new_manager.current_job_name = " ~ Phase Determination ~ "; // Pass in a nice job name.
      new_manager.current_user_job_name = job_name; // also give it the user defined job name for easy reference.
      new_manager.infinity_mode = infinity_mode; // Pass in if we really want to commit to this.
      new_manager.quiet_mode = no_display;
      new_manager.my_consensus_phase_settings = phase_consensus_phase_settings; // Pass in IPA settings.
      new_manager.perform_ipa_experiment(phasing_settings, n_threads);
      
      // all done now with phase determination
      
      if (new_manager.consensus_solutions_is_valid)
      {
        int iterates = new_manager.consensus_solutions.size();
        for (int i = 0 ; i < iterates ; i ++ )
        {
          list_of_solutions.push_back(new_manager.consensus_solutions[i]); // Copy the array into our list of solutions... we push back here, and slowly add to the list 
        }
        // It's worth printing some summary information for the runs associated with each envelope, as there may be a consensus phase set generated. 
        // However we don't do this for the last envelope on the list - as we will be printing the summary again at the end.
        if ((mask_id - 1) == list_of_working_masks.size())
        {
          print_summary_information_page(std::cout);
          print_summary_information_page(summary_log);
        }
      }
    }
  }
  
  // Generally run automatically during JOB PHASE DETERMINATION, but we also have the capability to run it in a standalone fashion
  
  if (JOB_phase_consensus)
  {
    // DO PHASE CONCENSUS
    std::cout << "RUNNING JOB: PHASE CONCENSUS IN..." << std::endl;
    for (int i = 3; i > 0; i--)
    {
      std::cout << i << "..." << std::endl;
      std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }
    std::cout << std::endl;
    
    bool consensus_found = run_phase_consensus();
    if (!consensus_found)
    {
      // Indicate failure of the procedure.
      std::cout << "Phase consensus did not identify any highly consistent phase sets, and no consensus was generated." << std::endl;
      summary_log << util::parse_line(cl, "Phase consensus did not identify any highly consistent phase sets, and no consensus was generated.");
    } else 
    {
      summary_log << util::parse_line(cl, std::to_string(list_of_solutions.size()) + " highly consistent phase sets were identified");
    }
  }
  
  /*
   From here we have succeeded and found a consensus phase set, and should print some sort of comprehensive summary information.
   */
  
  print_summary_information_page(std::cout);
  print_summary_information_page(summary_log);
  
  return true;
}

bool exp_manager::run_initialisation()
{
  //This function generates the required folder structures, and handles the one time estimation of the overall scale and B-factor.
   
  //We start by creating a guaranteed folder structure, which other parts of the program use. The folder structure is based on the allowed jobs, and within those folders, subdirectories labelled based on the user input job_name.
    
  // We currently do a cheeky "mkdir" system call to make the required folder structures, this has the benefit of working on most OS's. 
  
  // Standard libraries for filesystem aren't included in MacOS 10.14, We need MacOS 10.15 or later to have a more 'standard' way of solving this problem.
  // TODO: implement a more standard way of creating the folder structure, don't worry about versions of Mac OS preceding 10.15 ..
  
  /*
   
   ->root /envelope_determination
              /job_name
          /envelope_consensus
              /job_name
          /phasing_determination
              /job_name
          /phase_consensus
              /job_name
          /final_outputs/
              /job_name
   */
  

  // Create the main output folders:
  system("mkdir envelope_determination");
  system("mkdir envelope_consensus");
  system("mkdir phase_determination");
  system("mkdir phase_consensus");
  
  //Create sub-folders within the main ones, based on the requested job_name.
  {
    char input[512];
    std::strcpy(input,"mkdir envelope_determination/");
    std::strcat(input, job_name.c_str());
    system(input);
  }
  
  {
    char input[512];
    std::strcpy(input,"mkdir envelope_consensus/");
    std::strcat(input, job_name.c_str());
    system(input);
  }
  
  {
    char input[512];
    std::strcpy(input,"mkdir phase_determination/");
    std::strcat(input, job_name.c_str());
    system(input);
  }
  
  {
    char input[512];
    std::strcpy(input,"mkdir phase_consensus/");
    std::strcat(input, job_name.c_str());
    system(input);
  }
  
//  {
//    char input[100];
//    std::strcpy(input,"mkdir final_output/");
//    std::strcat(input, job_name.c_str());
//    system(input);
//  }
  
  // Calculate scale, overall isotropic and anisotropic B-factor utilising C++ driver for Fortran module Rogers
  // The method involves fitting a Gaussian to the Patterson origin peak
  // This analysis needs to be done before executing envelope or phase determination runs.

  util::print_section_line(cl,summary_log);
  util::print_empty_line(cl, summary_log);
  summary_log << util::parse_line(cl, " ~ INITIALISATION ~ ");
  util::print_empty_line(cl, summary_log);
  util::print_section_line(cl,summary_log);

  if (found_input_file)
  {
    
    std::cout << "Estimating scale, the overall anisotropic B-factor, and the equivalent overall isotropic B-factor by the method of Rogers" << std::endl;
        
    // It also creates model-based estimates for  <I> and puts these in a shared data object. 
    target_b_factor_data = estimate_B_from_diffraction_data(measured_f_all, summary_log); 
    
    std::cout << "The overall isotropic B-factor for the dataset is estimated as: " << target_b_factor_data.B_ISO << " Å^2" << std::endl;  
    // TODO: Report scale and overall anisotropic B-factor here as well, at the moment they are restricted to the Rogers output and easily missed ...
    
    // Add the information about scale and B into the settings files for envelope and phase determination
    envelope_settings.target_b_factor_data = b_factor_data(target_b_factor_data); // Copy the B-factor data into the settings objects.
    phasing_settings.target_b_factor_data = b_factor_data(target_b_factor_data); // Copy the B-factor data into the settings objects.
  }
  else
  {
    std::cout << "No input data were found (the mtz file could not be located) so no B-factor can be estimated Abort." << std::endl;
    return false;
  }
  
  // Time to ensure our settings structures have the predicted <I> and the measured Structure factor amplitudes embedded in them
  // This is required before we pass on copies of the settings to the workers, via threading.
   
  initialise_hkl_info_and_data_objects_for_settings(envelope_settings);
  initialise_hkl_info_and_data_objects_for_settings(phasing_settings);
  
  // TODO: Generate an envelope from known phases if provided... Set it in the settings. This could be useful for automated debugging and program development.

  return true;
}

bool exp_manager::run_envelope_consensus()
{
  
  // Create custom log file for this envelope_consensus run.
  std::ofstream envelope_log_stream;
  char log_filename[512];
  snprintf(log_filename, sizeof(log_filename), "%s%s_%0.2f_env_cons_log.txt", dir_envelope_consensus.c_str(), job_name.c_str(), envelope_consensus_mask_settings.epsilon);
  envelope_log_stream.open (log_filename);
  if (!envelope_log_stream.good())
  {
    std::cerr << "Bad envelope consensus log file generation!";
    return false; // Unsuccesful run
  }
  util::version_fingerprint(envelope_log_stream);
  envelope_log_stream << "Start of Envelope Consensus Job." << std::endl;
  
  
  // Retrieve a list of files from the envelope_determination/outputs/ directory ...
 
  // TODO If an envelope determination job is killed and relaunched with the same job name, without file ID offset, the masks and log files get overwritten, but mask list file is not
  // This will result in duplication of the filenames within mask list, and failure of the subsequent attempt to determine an envelope consensus.
  // Need to devise soem catch for this issue
  
  
  if (envelope_consensus_mask_settings.mask_list_file == "DEFAULT")
  {
    // The name of the mask list file has not been set - try finding the default for this job name.
    envelope_consensus_mask_settings.mask_list_file =  dir_envelope_determination + "mask_list.txt";
  }
  
  std::ifstream input_file(envelope_consensus_mask_settings.mask_list_file);
  std::vector<clipper::String> file_list;
  clipper::String line;
  
  if(!input_file) {
    envelope_log_stream << "File list could not be opened: " << envelope_consensus_mask_settings.mask_list_file << "\n" << std::endl;
    return 1;
  }
  while(std::getline(input_file, line)) {
    file_list.push_back(line);
  }
  
  if (file_list.size() < 1) {
    envelope_log_stream << " No envelopes were specified in the envelope_determination/outputs/mask_list.txt file!" << std::endl;
  }
  
  envelope_consensus_mask_settings.file_list = file_list; // wont be out of scope memory-wise so this should be fine rather than copying.
  
  std::vector<std::string> list_of_consensus_files; // Create somewhere to store a list of the consensus envelopes
  
  
  // Perform envelope clustering and "averaging"
  ipa_functions::compute_consensus_mask(envelope_consensus_mask_settings,
                                        envelope_log_stream,
                                        dir_envelope_consensus, // important for job name.
                                        list_of_consensus_files);
 
                                        
  bool found_consensus_envelopes = false;
  if (list_of_consensus_files.size() > 0) // If any output files were generated...
  { 
    util::print_empty_line(cl,summary_log);
    summary_log << util::parse_line(cl, "A total of " + std::to_string(list_of_consensus_files.size()) + " consensus envelopes were generated.");
    found_consensus_envelopes  = true;
    phasing_settings.input_working_mask = true; // Must set this true now, as irrespective of whether any envelopes were input by the user, we have some now
		
    for (int i = 0 ; i < list_of_consensus_files.size() ; i ++ ) // For all new files
    {
      bool do_add = true; // Determines if we are adding this file, defaults to true.
      
      for (int j = 0; j < list_of_working_masks.size() ; j ++ ) // Check against all files already added.
      {
        if (list_of_consensus_files[i] == list_of_working_masks[j])
        {
          do_add = false; // found a duplicate, don't add
          break;
        }
      }
      
      if (do_add) // Add if unique.
      {
        list_of_working_masks.push_back(list_of_consensus_files[i]);
        summary_log << util::parse_line(cl, list_of_consensus_files[i]);

      }
    }
  // TODO Currently we just work through each mask generated in turn, and execute the phase determination step
  // In the future use some heuristic criteria (connectivity, conformity with the solvent fraction etc etc) to put the masks in some rank order, with the "best" masks used first

  
  util::print_empty_line(cl,summary_log);
    
  } 
  else 
  {
    summary_log << util::parse_line(cl, "Envelope Consensus step failed to generate any consensus envelopes due to a lack of agreement.");

    return false; // No consensus masks were generated.
  }
  
  envelope_log_stream.close();
  
  return true;
  //Successful closure.
}

bool exp_manager::run_phase_consensus()
{
  // Retrieve a list of files from the phase_determination/outputs/ directory 
  
  if (phase_consensus_phase_settings.phase_list_file == "DEFAULT")
  {
    // The name of the phase list file has not been set - try finding the default for this job name.
    phase_consensus_phase_settings.phase_list_file =  dir_phase_determination + "phase_list.txt";
  }
  
  std::ifstream input_file(phase_consensus_phase_settings.phase_list_file);
  std::vector<clipper::String> file_list;
  clipper::String line;
  
  while(std::getline(input_file, line)) {
    file_list.push_back(line);
  }

  util::print_section_line(cl,summary_log);
  util::print_empty_line(cl, summary_log);
  summary_log << util::parse_line(cl, " ~ PHASE CONSENSUS SUMMARY ~ ");
  util::print_empty_line(cl, summary_log);
  util::print_section_line(cl,summary_log);
  util::print_empty_line(cl, summary_log);

  if (file_list.size() < 1)
  {
    std::cout << " No phase sets were found in the phase_determination/outputs/phase_list.txt file!" << std::endl;
    summary_log << util::parse_line(cl, " No phase sets were found in the phase_determination/outputs/phase_list.txt file!");
    util::print_empty_line(cl, summary_log);
    util::print_section_line(cl,summary_log);
    summary_log << std::endl;
    return false;
  }
  
  
  //Create the consensus phasing object.
  consensus_worker my_consensus_worker(phase_consensus_phase_settings);
  my_consensus_worker.is_idle = true;
  
  
  std::cout << " Total size of input file list = " << file_list.size() << std::endl;
  summary_log << util::parse_line(cl," Total number of input phase sets to be analyzed: " + std::to_string(file_list.size()));

  for (int i = 0 ; i < file_list.size() ; i ++ )
  {
    //std::cout << "Adding file \"" << file_list[i] << "\" to consensus worker..." << std::endl;
    bool success = my_consensus_worker.add_phase_input_to_list(file_list[i]);
    if (!success)
    {
      std::cout << "Failed to add file to consensus worker: " << file_list[i] << " " << std::endl;
    } else 
    {
      summary_log << util::parse_line(cl," File " + std::to_string(i) + ": " +file_list[i]);
    }
  }

  util::print_empty_line(cl, summary_log);
  
  std::atomic<float> myfloat; // gotta pass something.
  
  // Do a consensus.
  my_consensus_worker.check_for_consensus_phases(myfloat);
  
  my_consensus_worker.return_list_of_consensus_solutions(list_of_solutions, std::cout);
  my_consensus_worker.return_list_of_consensus_solutions(list_of_solutions, summary_log);

  return (list_of_solutions.size() > 0); // return true if solutions were found.
  //return true; // Can't ever get here because of if /else, but stops compile warning.
}

void exp_manager::print_summary_information_page(std::ostream& custom_outlog)
{
  /*
   Output the summary information from the most recent phase consensus run.
   This can be called at any time, but mainly for calling after a phase consensus or phase determination job.
   */
  
  //  Start with  reminder of the clustering parameters (epsilon and MinPoints)
  //  Then probably report (sequentially, for each envelope in list_of_working_masks)

  //  Where do I locate the envelope used to initiate phasing?
  //  The number of clusters generated with this envelope (hopefully 1 !!! but technically we need to loop over all the clusters here)
  //  The Cluster membership (in terms of run identifiers or filenames)
  //  The sample circular variance for the cluster
  //  Where do I find the consenus mtz / map files for the cluster ?

  util::print_section_line(cl, custom_outlog);
  util::print_empty_line(cl, custom_outlog);
  custom_outlog << util::parse_line(cl, "               SUMMARY OF SOLUTIONS:   ");
  util::print_empty_line(cl, custom_outlog);

  
  for (int i = 0 ; i < list_of_solutions.size(); i ++ )
  {
    summarise_cluster_data(list_of_solutions[i], custom_outlog);
  }
  
  //util::print_empty_line(cl, std::cout);
  util::print_section_line(cl, custom_outlog);

  custom_outlog << std::endl;
}

std::ostream& exp_manager::summarise_cluster_data(cluster_data& data, std::ostream& outlog)
{
  /*
   Helper function to summarise a single cluster data object.
   
   struct cluster_data
   {
     float  epsilon;
     int cluster_id;
     int n_members; // Redundant for a call to list_of_members.size()?
     std::vector<std::string> list_of_members; // We remember which files went into this cluster, now.
     std::string consensus_filename; // The output average of the members.
     std::string consensus_filename_inverse = "NONE"; // if space group allows, this will be populated. else its None.
     
     float scv;
     float mean_phase_diff_all_with_known;
     float mean_phase_diff_all_with_known_inverse;
   };
   
   */
  
  util::print_section_line(cl, outlog);
  outlog << util::parse_line(cl, "Consensus Phase Summary: ");
  outlog << util::parse_line(cl, "Solution File:");
  //  For the achiral space groups we write out the solution and its inverse, in the original spacegroup
  //  For the chiral spacegroups  we write out the solution and its inverse, switching to the enantimorphic space group when we invert
  outlog << util::parse_line(cl, data.consensus_filename);
  outlog << util::parse_line(cl, "Inverse of Solution:");
  outlog << util::parse_line(cl, data.consensus_filename_inverse);
  util::print_empty_line(cl, outlog);
  outlog << util::parse_line(cl, "Both the solution and its inverse are written ... the correct hand will need to be determined by inspection");
  outlog << util::parse_line(cl, "If the space group is chiral, inversion will be associated with a change of space group");
    
  int number_of_members = data.list_of_members.size(); // Grab this from here for redundancy from n_members.
  
  util::print_empty_line(cl, outlog);
  outlog << util::parse_line(cl, "Number of members: " + std::to_string(number_of_members));
  outlog << util::parse_line(cl, "Sample Circular Variance: " + std::to_string(data.scv));
  outlog << util::parse_line(cl, "Mean Phase Difference across cluster: " + std::to_string(data.mean_phase_difference_across_cluster));
  if (data.mean_phase_diff_all_with_known != 0) // Why not check if a known phase set was input in settings ? Then the comparison must have been made. Seems safer than this numerical check 
  {
    outlog << util::parse_line(cl, "Comparison with the Known Phase Set: ");
    outlog << util::parse_line(cl, "Mean Phase Difference: " + std::to_string(data.mean_phase_diff_all_with_known));
    bool space_group_is_chiral = util::space_group_is_chiral(phasing_settings.hkl_target.spacegroup().spacegroup_number());
    if (!space_group_is_chiral) {outlog << util::parse_line(cl, "Mean Inverse Phase Difference: " + std::to_string(data.mean_phase_diff_all_with_known_inverse));}
  }
  util::print_empty_line(cl, outlog);
  outlog << util::parse_line(cl, "List of Members: ");
  for (int i = 0 ; i < number_of_members ; i ++ )
  {
    outlog << util::parse_line(cl, " ~ " + data.list_of_members[i]);
  }
  
  util::print_empty_line(cl, outlog);
  
  return outlog;
}

int exp_manager::import_settings_from_file(ipa_settings& env_settings, ipa_settings& pha_settings, std::string filename, std::ostream& outlog)
{
  /*
   This function will parse a experiment.params file, and populate the entirety of an ipa_settings struct with default values, or those present within the file.
   
   Some checks are performed at the end to ensure that settings are sensible and no human error has been made.
   
   An integer is returned,and is indicative of error found:
   0 = No Error
   1 = Could not find or open the experiment.params file.
   2 = Bad Algorithm Instructions.
   3 = Unrecognised Job.
   4 = Undefined Settings input.
   */
  
  // Set some Boolean varaibles so we can easily track some essential input 
  JOB_initialisation = false;
  JOB_envelope_determination = false;
  JOB_envelope_consensus = false;
  JOB_phase_determination = false;
  JOB_phase_consensus = false;
  
  n_bad_inputs_count = 0;
  found_input_file = false;
  found_input_file_columns = false;
  custom_job_id = false;
  found_job_list = false;
  found_beta_params = false;
  found_filter_params = false;
  found_solvent_radius_params = false;
  core_b_factor_provided = false;
  core_solvent_fraction_provided = false;
  
  list_of_working_masks.clear(); // ensure the list of working masks is empty.
  
  n_threads = std::thread::hardware_concurrency(); // default number of threads is the total number of possible threads for this machine.
  
  // Open file
  std::fstream fparams(filename);
  if (!fparams.good())
  {
    outlog << "Could not find or open the experiment parameter file, a default file will be generated" << std::endl;
    return 1; // Could not find or open the params file.
  }
  
  outlog << "Found experiment parameter file, parsing..." << std::endl;
  
  std::string line;
  int current_line = 0; // Debug Counter indicating the line we are on.
                        // Read the file line by line...
  
  while (std::getline(fparams, line))
  {
    current_line++;
    
    remove_comments_from_string(line);
    
    //outlog << "Line " << current_line << ":" << line << std::endl;
    
    if (line.length() > 0) // Only safe lines here...
    {
      char delim = '-';
      if (line.front() == delim)
      {
        std::vector<std::string> line_list;
        parse_string_into_string_array(line_list, line, deliminators);
        
        //  DEBUGGING EACH WORD IF NEEDED.
        //outlog << "The list of words are: ";
        //for (int i = 0 ; i < line_list.size(); i ++)
        //{
        //    outlog << " [" << i << "]: " << line_list[i] << "   ";
        //}
        //outlog << std::endl;
        
        if (line_list[0] == "-JOBS") {
          // Read the JOBS setting here.
          read_job_list_from_file(line_list, outlog);
        }
        else if (line_list[0] == "-core") {
          // Read the core setting here.
          read_string_with_core_label(line_list,outlog);
        }
        else if (line_list[0] == "-envelope")
        {
          read_string_setting_into_struct(env_settings,line_list,outlog);
          
        }
        else if (line_list[0] == "-phasing")
        {
          read_string_setting_into_struct(pha_settings,line_list,outlog);
          
        }
        else if (line_list[0] == "-mskconsensus")
        {
          // doesn't need to pass the settings as there's only ever 1.
          read_mask_consensus_setting_into_struct(line_list,outlog);
        }
        else if (line_list[0] == "-phsconsensus")
        {
          // doesn't need to pass the settings as there's only ever 1.
          read_phase_consensus_setting_into_struct(line_list,outlog);
        }
        else
        {
          // abort, unrecognised setting type...
          outlog << "[Line #" << current_line << "] Unrecognised setting type: " << line_list[0] << std::endl;
          n_bad_inputs_count++;
          return 4;
        }
      }
    }
  } // end of While Loop;
  
  // Set some job specific constants.
  envelope_settings.working_dir = "envelope_determination/" + job_name + "/";
  phasing_settings.working_dir = "phase_determination/" + job_name + "/";
  
  phase_consensus_phase_settings.working_dir = "phase_consensus/" + job_name + "/";
  
  envelope_settings.job_is_envelope_determination = true;
  phasing_settings.job_is_phase_determination = true;
  fparams.close();
  
  return 0;
}

// TODO ... rename this - is "validate" the best name for a function that also "imports" ... maybe import_FsigF_data would be better
bool exp_manager::validate_input_mtz()
{
  /*
   Imports the Structure Factor amplitudes from a mtz file, using the user supplied column labels, this should work, and if not, is a first step for the user to correct their inputs.
   At the same , we import all the associated crystallographic information into hkl_target, which is passed on via struct to all works
   This function should only be called AFTER the settings files have been initialised from the parameter file.
   The filename and column labels are core settings, which are pushed to both phase and envelope settings by default. Hence we can access through either
   */
  
  std::cout << "\nReading data from the user specified file: " << envelope_settings.input_filename_target_data << std::endl;

  clipper::CCP4MTZfile mtzin;                                               // Declare new mtz object
  
  clipper::HKL_info hkl_input;
  clipper::HKL_data<clipper::data64::F_sigF> data_input;

  mtzin.open_read( envelope_settings.input_filename_target_data );          // open the mtz file for reading
  mtzin.import_hkl_info( hkl_input );                                      // Reads cell etc and initialises hkl_target.
  data_input.init(hkl_input, hkl_input.cell());                       // Initialises the F_SigF data to the hkl_target info.
  mtzin.import_hkl_data( data_input, envelope_settings.target_col_fo);  // Copies the F_sigF data into measured_f_all
  mtzin.close_read();                                                       // close the mtz file after reading.
  
  // Check if we need to truncate the data...
  double resolution;
  if (high_resolution_truncation > 0)
  {
    std::cout << "Truncating data to " << std::setprecision(2) << high_resolution_truncation << " Å resolution ..." << std::endl;
    resolution = high_resolution_truncation; // INSERT THE SETTING HERE.
  }
  else
  {
    resolution = hkl_input.resolution().limit(); // No change.
  }

  hkl_target.init(hkl_input.spacegroup(), hkl_input.cell(), clipper::Resolution(resolution));
  hkl_target.generate_hkl_list(); // Make the hkl reflections.
  measured_f_all.init(hkl_target, hkl_target.cell()); // Initialise an data object

  // Copy the data into the truncated dataset. Missing values are also copied.
  clipper::HKL_info::HKL_reference_index ih;
  for (ih = hkl_target.first(); !ih.last(); ih.next() )
  {
    measured_f_all[ih].f() = data_input[ih.hkl()].f();
    measured_f_all[ih].sigf() = data_input[ih.hkl()].sigf();
  }

  std::cout << "Reading of data is complete" << std::endl; 
  return true;
}

bool exp_manager::validate_settings(ipa_settings& settings)
{
  // Perform some validation of the envelope and phasing settings and compute some essential information - e.g. the total number of iterates
    
  settings.n_total_iterates = 0;
  int new_total_iterates = 0;
  int n_rules = settings.UpdateRuleRegime.size();
  for (int i = 0; i < n_rules; i++)
  {
    new_total_iterates += settings.UpdateRuleRegime[i].n_iterates;
  }
  settings.n_total_iterates = new_total_iterates;
  
  if (settings.n_apodization_last_iterate > settings.n_total_iterates)
  {
    std::cout << "Final Apodization step was set beyond the final iterate, Final apodization step will be set to be the final iterate" << std::endl;
    settings.n_apodization_last_iterate = settings.n_total_iterates;
    n_bad_inputs_count++;
  }
  
  if (phase_consensus_phase_settings.phase_list_file == "DEFAULT")
  {
    // The phase list file has not been set - try find one thats default for this job name.
    phase_consensus_phase_settings.phase_list_file =  dir_phase_determination + "phase_list.txt";
    // else don't worry about it, the user can get an error if its set wrong when its failed to be opened.
  }
  
  // Some basic checks on the data and input
  
  if ( settings.target_col_fo == "NONE" )
  {
    std::cerr << "No F's provided for the target structure." << std::endl;
    n_bad_inputs_count++;
    return false;
  }
  
  /*  //TODO:  Is this still desired behaviour...
  if ( settings.input_starting_phases && (settings.n_runs > 1 ))
  {
    std:cerr << " Phase set is input - only a *single* run will be performed "; settings.n_runs = 1;
    n_bad_inputs_count++;
  }
   */
    
    
  if ( ( JOB_phase_determination && (settings.n_iterations_fix_mask > 0) ) &&  !(settings.input_working_mask  ||  JOB_envelope_consensus) ) // We have to input a mask or generate one, if we are going to fix the mask during phase determination.
  {
    std::cerr << "To fix the mask during phase determination, a mask needs to be input, or the mask needs to be generated via an envelope consensus job" << std::endl;
    n_bad_inputs_count++;
    return false;
  }
  
  if (settings.compute_mask_interval < 1)
  {
    std::cerr << " Output mask interval cannot be less than 1 - resetting to 1";
    settings.compute_mask_interval = 1;
    n_bad_inputs_count++;
  }
    
  if ( ( settings.mask_weighting_factor_mean < 0.0) || (  settings.mask_weighting_factor_mean > 1.0) )
  {
    std::cerr << "mask weighting factor should be  > 0 and < 1\n";
    n_bad_inputs_count++;
    return false;
  }
  
  if ( (settings.erase_islands_threshold < 0.0) || ( settings.erase_islands_threshold > 1.0) )
  {
    std::cerr << "Threshold value for island erasure should be > 0 and < 1\n";
    n_bad_inputs_count++;
    return false;
  }
  
  if ( (settings.fill_voids_threshold < 0.0) || ( settings.fill_voids_threshold > 1.0) )
  {
    std::cerr << "Threshold value for void filling should be > 0 and < 1\n";
    n_bad_inputs_count++;
    return false;
  }
  
  
  // If phase_output_interval or mask_output_interval = 0 then we only output the phases/mask on the final cycle
  // Otherwise we output the phases/mask at the specified intervals
  // If we don't output the phases at every cycle then turn off phase reversion, since this requires a complete record of the progression through the iterates.
  
  //  if ( (output_phases == false) || ( (output_phases == true) && (phase_output_interval != 1) ) ) {revert_to_best_phase_set = false;}
  
  // if (calculate_mask_agreement && !input_known_mask)
  //   { std::cerr << "Sorry, I can't calculate mask agreement, if no known mask is input" << std::endl; return 1; }
  
  //  should also check here for failure to assign solvent content, B-factor etc, etc, etc
  
  return true;
}

bool exp_manager::read_job_list_from_file(std::vector<std::string>& line, std::ostream& outlog)
{
  /*
   Determines the jobs to be performed from the experiment.params line, assumes a -JOBS flag has been found at position 0 of the line.
   */
  for (int i = 1; i < line.size(); i++)
  {
    // for each job, turn on if present.
    if (line[i] == "INIT")
    {
      JOB_initialisation = true; // Remove this, as INIT has become obligate? Or continue to allow it, so user could perform Rogers analysis in stand-alone fashion?
      continue;
    }
    
    if (line[i] == "ENVELOPE_DETERMINATION")
    {
      JOB_envelope_determination = true;
      continue;
    }
    if (line[i] == "ENVELOPE_CONSENSUS")
    {
      JOB_envelope_consensus = true;
      continue;
    }
    if (line[i] == "PHASE_DETERMINATION")
    {
      JOB_phase_determination = true;
      continue;
    }
    if (line[i] == "PHASE_CONSENSUS")
    {
      JOB_phase_consensus = true;
      continue;
    }
    
    // If we get to this point in the loop, the job is unrecognised.
    outlog << std::endl << "BAD JOB INPUT: " << line[i] << std::endl;
    return false;
  }
  
  found_job_list = true;
  return true;
}

bool exp_manager::read_string_with_core_label(std::vector<std::string>& lines, std::ostream& outlog)
{
  /*
   parse a -core labelled setting. Typically these are reserved for critical input. Bfactor, Solvent content, and directory locations to input files.
   */
  
  if (lines.size() == 1)
  {
    /*
     We can't get here if lines size == 0, as it needed the core label, however if there is simply only a core label, throw an error.
     */
    std::cout << "ERROR: -core labelled parameter line contains no further input" << std::endl;
    n_bad_inputs_count++;
    return false;
  }
  
  if (lines[1] == "no-display")
  {
    no_display = true;
    if (lines.size() > 2) no_display = convert_string_into_bool(lines[2]);
    return true;
  }
  
  if (lines[1] == "input-data-mtz-filename")
  {
    mtz_input_file = lines[2]; // Storing the filename this way seems redundant, but is useful later for analysis of the overall B-factor
    envelope_settings.input_filename_target_data = mtz_input_file;
    phasing_settings.input_filename_target_data = mtz_input_file;
    
    // TODO: eliminate this clumsy method of passing the data to the envelope consensus function - just a throwback to its origns as a standalone proram 
    envelope_consensus_mask_settings.input_filename_target_data = mtz_input_file;
    envelope_consensus_mask_settings.input_Fourier_amplitudes = true;
    
    found_input_file = true;
  }
  
  if (lines[1] == "input-data-fsigf-column-labels") // Labels for the measured structure factor amplitudes and their standard deviations
  {
    if (lines.size() > 2)
    {
      envelope_settings.target_col_fo = lines[2];
      phasing_settings.target_col_fo = lines[2];
      
    // TODO: eliminate this clumsy method of passing the data to the envelope consensus function - just a throwback to its origns as a standalone proram 
      envelope_consensus_mask_settings.target_col_fo = lines[2];
      
      found_input_file_columns = true;
    }
    else
    {
      std::cout << "The column labels for the input MTZ were not defined at -core input-data-fsigf-column-labels " << std::endl;
      n_bad_inputs_count++;
      return false;
    }
  }
  
  if (lines[1] == "input-data-phiw-column-labels") // Labels for the starting phases and weights, if we are not doing ab initio phase determination.
  {
    if (lines.size() > 2)
    {
      envelope_settings.target_col_pw = lines[2];
      phasing_settings.target_col_pw = lines[2];

      envelope_settings.input_starting_phases = true;
      phasing_settings.input_starting_phases = true;
    }
  }
  
  if (lines[1] == "known-solution-mask-filename") // A file containing the known mask (i.e. a mask calculated from the solution). For testing algorithm performance
  {
    // If the user specifies a known mask in the -core section, add it everywhere 
    envelope_settings.input_filename_known_mask = lines[2];
    phasing_settings.input_filename_known_mask = lines[2];
    envelope_consensus_mask_settings.input_filename_known_mask = lines[2];
    
    envelope_settings.input_known_mask = true;
    phasing_settings.input_known_mask = true;
    envelope_consensus_mask_settings.input_known_mask = true;
  }
  
  if (lines[1] == "known-solution-phase-filename") // A file containing the known phases (i.e. phases calculated from the solution). For testing algorithm performance
  {
    // If the user specifies a known phase set in the -core section, add them everywhere
    envelope_settings.input_filename_known_phase_set = lines[2];
    phasing_settings.input_filename_known_phase_set = lines[2];
    phase_consensus_phase_settings.input_filename_known_phase_set = lines[2];
    
    envelope_settings.input_known_phases = true;
    phasing_settings.input_known_phases = true;
    phase_consensus_phase_settings.input_known_phases = true;
  }
  
  if (lines[1] == "known-solution-fp-column-labels")  // Labels for the known amplitudes and phases
  {
    if (lines.size() > 2)
    {
      envelope_settings.known_phase_set_fp_column_labels = lines[2];
      phasing_settings.known_phase_set_fp_column_labels = lines[2];
      envelope_settings.input_known_phases = true;
      phasing_settings.input_known_phases = true;
      
      phase_consensus_phase_settings.known_phase_set_fp_column_labels = lines[2]; 
      phase_consensus_phase_settings.input_known_phases = true;
    }
    else
    {
      std::cout << "The column labels for the known phase MTZ were not defined at -core known-solution-fp-column-labels " << std::endl;
      n_bad_inputs_count++;
      return false;
    }
  }
  
  if (lines[1] == "input-starting-envelope") // Input a starting envelope for the phase determination step. User can repeatedly add as many envelopes as they wish 
  {
    if (lines.size() > 2)
    {
      phasing_settings.input_filename_target_mask_working = lines[2];
      list_of_working_masks.push_back(lines[2]);
      user_input_mask_count++; // Increment this so we know how many user input masks there are.
      
      // TODO: CHECK DIRECTORY/FILE IS VALID - need a general solution for this.
    
      //envelope_settings.input_working_mask = true;
      phasing_settings.input_working_mask = true;
    }
  }
  
  if (lines[1] == "job-name") // A name for our job
  {
    job_name = lines[2];
    envelope_settings.job_name = job_name;
    phasing_settings.job_name = job_name;
    
    custom_job_id = true;
  }

  if (lines[1] == "truncate-data")
  {
    if (lines.size() > 2)
    {
      high_resolution_truncation = stof(lines[2]);
    }
  }
  
  
  if (lines[1] == "input-data-solvent-fraction") // The solvent fraction
  {
    envelope_settings.fraction_solvent = stof(lines[2]);
    phasing_settings.fraction_solvent = stof(lines[2]);
    
    envelope_consensus_mask_settings.fraction_solvent_target = stof(lines[2]);
    
    core_solvent_fraction_provided = true;
    return true;
  }
  
  if (lines[1] == "threads") // the number of threads 
  {
    if (lines.size() > 2)
    {
    n_threads = stoi(lines[2]);
    }
    else
    {
      std::cout << "The number of threads was not defined" << std::endl;
      n_bad_inputs_count++;
      return false;
    }
    return true;
  }
  
  if (lines[1] == "verbose")
  {
    bool verbose_setting = true; // defaults to true, if there is no explicit true/false assignment 
    if (lines.size() > 2) {
      verbose_setting = convert_string_into_bool(lines[2]);
    }
    envelope_settings.verbose = verbose_setting;
    phasing_settings.verbose = verbose_setting;
    envelope_consensus_mask_settings.verbose = verbose_setting;
    phase_consensus_phase_settings.verbose =verbose_setting;
  }
  
  if (lines[1] == "really-random")
  {
    bool random_setting = true; // defaults to true, if there is no explicit true/false assignment 
    if (lines.size() > 2) {
      random_setting = convert_string_into_bool(lines[2]);
    }
    envelope_settings.really_random = random_setting;
    phasing_settings.really_random = random_setting;
  }
  
  if (lines[1] == "statistics-interval")
  {
    if (lines.size() == 3)
    {
      envelope_settings.report_statistics_interval = stoi(lines[2]);
      phasing_settings.report_statistics_interval = stoi(lines[2]);
    }
    else
    {
      outlog << "Issue with core input for statistics interval" << std::endl;
      n_bad_inputs_count++;
      return false;
    }
  }
  
  if (lines[1] == "protein-sequence")
  {
    if (lines.size() == 3)
    {
      snprintf(protein_sequence_filename, sizeof(protein_sequence_filename), "%s", lines[2].c_str());
      protein_sequence_is_provided = true;
    }
  }
  if (lines[1] == "rna-sequence")
  {
    if (lines.size() == 3)
    {
      snprintf(rna_sequence_filename, sizeof(rna_sequence_filename),"%s", lines[2].c_str());
      rna_sequence_is_provided = true;
    }
  }
  if (lines[1] == "dna-sequence")
  {
    if (lines.size() == 3)
    {
      snprintf(dna_sequence_filename, sizeof(dna_sequence_filename),"%s", lines[2].c_str());
      dna_sequence_is_provided = true;
    }
  }
  
  if (lines[1] == "protein-n-chains")
  {
    if (lines.size() == 3)
    {
      protein_sequence_nchains = stoi(lines[2]);
    }
  }
  if (lines[1] == "dna-n-chains")
  {
    if (lines.size() == 3)
    {
      dna_sequence_nchains = stoi(lines[2]);
    }
  }
  if (lines[1] == "rna-n-chains")
  {
    if (lines.size() == 3)
    {
      dna_sequence_nchains = stoi(lines[2]);
    }
  }
  
  if (lines[1] == "protein-copies-in-asu")
  {
    if (lines.size() == 3)
    {
      protein_sequence_ncopies_in_asu = stoi(lines[2]);
    }
  }
  if (lines[1] == "dna-copies-in-asu")
  {
    if (lines.size() == 3)
    {
      dna_sequence_ncopies_in_asu = stoi(lines[2]);
    }
  }
  if (lines[1] == "rna-copies-in-asu")
  {
    if (lines.size() == 3)
    {
      rna_sequence_ncopies_in_asu = stoi(lines[2]);
    }
  }
  if (lines[1] == "infinity-mode")
  {
    infinity_mode = true;  /// Defaults to true if there is no explicit true/false assignment 
    if (lines.size() == 3) 
    {
      infinity_mode = convert_string_into_bool(lines[2]);
    }
  }
  return true;
}

bool exp_manager::read_string_setting_into_struct(ipa_settings& SettingsToUpdate, std::vector<std::string>& info, std::ostream& outlog)
{
  /*
   This function is designed to evaluate the string array for a line of an experiment.params file, and copy the appropriate values into the associated envelope or phasing struct. The passed std::vector of string information has already been checked to contain an appropriate header, and will be devoid of comments.
   */
  
  int max_size = info.size(); // Useful for safety checks and avoiding undefined behaviour.
  
  if (max_size < 3)
  {
    outlog << "Unexpected size for line: " << info[0] << " " << info[1] << std::endl;
  }
  
  // Initial "complex" line inputs for >1 setting inputs.
  if (info[1] == "algorithm")
  {
    bool success = read_string_algorithm_info_into_struct(SettingsToUpdate, info, outlog);
    return success;
  }
  else if (info[1] == "beta-param")
  {
    bool success = read_string_beta_params_info_into_struct(SettingsToUpdate.beta_params_look_up_table, info, outlog);
    SettingsToUpdate.found_beta_params = success;
    return success;
  }
  else if (info[1] == "filter-param")
  {
    bool success =  read_string_filter_params_info_into_struct(SettingsToUpdate.filter_radius_params_look_up_table, info, outlog);
    SettingsToUpdate.found_filter_params = success;
    return success;
  }
  else if (info[1] == "solvent-param")
  {
    bool success = read_string_solvent_params_info_into_struct(SettingsToUpdate.solvent_fraction_multiplier_params_look_up_table, info, outlog);
    SettingsToUpdate.found_solvent_fraction_multiplier_params = success;
    return success;
  }

  
  
  // Iterate through all other 'simple' settings via if/else. Ensure all have at least 1 name and 1 entry (i.e. size() >= 3)
  else if ( info[1] == "mtzin-target" && info.size()>=3) {
    SettingsToUpdate.input_filename_target_data = info[2]; // Name of the input mtz file carrying the data for the target structure
  } else if ( info[1] == "target-fo" ) {
    SettingsToUpdate.target_col_fo = info[2]; // Column path for observed F and sigma in the mtz
  } else if ( info[1] == "n-runs" ) {
    SettingsToUpdate.n_runs = clipper::String(info[2]).i(); // Total runs to perform.
  } else if ( info[1] == "target-phifom" ) {
    SettingsToUpdate.target_col_pw = info[2]; // Column path for Phase and Weight (Figure of Merit) in the mtz- relevant where we are improving an existing phase set - not performning ab initio phase determination
    SettingsToUpdate.input_starting_phases = true;
  } else if ( info[1] == "target-phifom-calc" ) {
    SettingsToUpdate.target_col_pw_calc = info[2]; // Column path for calculated Phase and Weight (Figure of Merit) in the mtz - this is useful for checking algorithm performance
    SettingsToUpdate.input_known_phases = true;
  } else if ( info[1] == "mskin-target-working" || info[1] == "input-starting-envelope") { // (mskin is included for backwards compatibility with old param files for now: 12.09.22)
    SettingsToUpdate.input_working_mask = true;
    SettingsToUpdate.input_filename_target_mask_working = info[2]; // Name of the input working mask file for the target structure.
    list_of_working_masks.push_back(info[2]); // Lets append to the list of working masks to go through.
    user_input_mask_count++; // Increment this so we know how many user input masks there are.
  } else if ( info[1] == "mskin-target-known" ) {
    SettingsToUpdate.input_known_mask = true;
    SettingsToUpdate.input_filename_known_mask = info[2]; // Name of the input known mask file for the target structure. Useful for checking algorithm performance
                                                                 // Can also now be set in -core settings! Leaving here for experiment.params legacy option.
  } else if ( info[1] == "n-iterations-fix-mask" ) {
    SettingsToUpdate.n_iterations_fix_mask = clipper::String(info[2]).i(); // if non-zero then the input working mask will be fixed for the specified number of iterations (otherwise it is continually updated)
  } else if ( info[1] == "n-sphere-phase-generation" ) {
    SettingsToUpdate.n_sphere_phase_generation = clipper::String(info[2]).i(); // Number of spheres to employ in real space method for generating starting phases
    SettingsToUpdate.fourier_space_phase_generation = false;
  } else if ( info[1] == "rad-sphere-phase-generation" ) {
    SettingsToUpdate.rad_sphere_phase_generation = clipper::String(info[2]).f64(); // Radius of spheres to employ in real space method for generating starting phases
    SettingsToUpdate.fourier_space_phase_generation = false;
  } else if ( info[1] == "max-attempts" ) {
    SettingsToUpdate.max_attempts = clipper::String(info[2]).i(); // Maximum number of attempts at random sphere positioning before we give up
  } else if ( info[1] == "prot-density" ) {
    SettingsToUpdate.expected_mean_density_protein = clipper::String(info[2]).f64(); // Expected mean electron density in protein region (electrons/Angstrom^3)
  } else if ( info[1] == "solvent-density" ) {
    SettingsToUpdate.expected_mean_density_solvent = clipper::String(info[2]).f64(); // Expected mean density in solvent region (electrons/Angstrom^3)
  } else if ( info[1] == "probability-threshold" ) {
    SettingsToUpdate.p_threshold = clipper::String(info[2]).f64(); // Probability cutoff used when deciding whether to accept the estimated values for the missing data
  } else if ( info[1] == "replace-largei-with-expected-value" ) {  
    SettingsToUpdate.intensity_multiplier = clipper::String(info[2]).f64(); // scale factor for the expected value
    SettingsToUpdate.replace_largeI_with_square_of_expectedF = true; 
  } else if ( info[1] == "replace-largei-with-attenuatedi" ) {
    SettingsToUpdate.intensity_multiplier = clipper::String(info[2]).f64(); // scale factor for the intensity estimate
    SettingsToUpdate.replace_largeI_with_square_of_expectedF = false;
  } else if ( info[1] == "beta-function" ) {
    SettingsToUpdate.beta_function = clipper::String(info[2]).i(); // Integer indicating which function is use for calculation of beta
  } else if ( info[1] == "linear-freq-change" ) {
    SettingsToUpdate.linear_freq_change = true;
  } else if ( info[1] == "linear-period-change" ) {
    SettingsToUpdate.linear_freq_change = false; // not strictly needed - this is the default setting
  } else if ( info[1] == "remove-low-res-data" ) {
    SettingsToUpdate.remove_low_res_data = true;
    SettingsToUpdate.low_res_cutoff = clipper::String(info[2]).f64(); // Cutoff for removal of low resolution data
  } else if ( info[1] == "really-random" ) {
    SettingsToUpdate.really_random = true;
    if (info.size()>2) {
      // check explicitness
      SettingsToUpdate.really_random = (info[2] == "true" || info[2] == "True" || info[2] == "TRUE");
    }
  } else if ( info[1] == "use-test-set" ) {
    SettingsToUpdate.use_test_set = true;
  } else if ( info[1] == "estimate-missing-amplitudes-on-first-iterate" ) {
    SettingsToUpdate.estimate_missing_amplitudes_on_first_iterate = true;
  } else if ( info[1] == "n-test-sets" ) {
    SettingsToUpdate.use_test_set = true;
    SettingsToUpdate.n_test_sets = clipper::String(info[2]).i();
  } else if ( info[1] == "mask-redefinition-interval" ) {
    SettingsToUpdate.compute_mask_interval = clipper::String(info[2]).i(); // Defines how many iterations elapse before redefining the mask
    
    
  } else if ( info[1] == "mask-from" ) { // Slightly more complicated.
    if (info[2] == "variance")
    {SettingsToUpdate.compute_mask_from_local_variance = true;} else
      if (info[2] == "mean")
      {SettingsToUpdate.compute_mask_from_local_mean = true;} else
        
        if (info.size() > 3 && (info[2] == "mean-and-variance" || info[2] == "variance-and-mean"))
        {
          SettingsToUpdate.compute_mask_from_local_mean = false;
          SettingsToUpdate.compute_mask_from_local_variance = false;
          SettingsToUpdate.compute_mask_from_local_mean_and_variance = true;
          SettingsToUpdate.mask_weighting_factor_mean = clipper::String(info[4]).f64();
        } else {
          outlog << "WARNING ... Unrecognized argument:\t";
          for (int i = 0; i<info.size(); i++)
          {
            outlog << info[i] << " ";
          }
          outlog << std::endl;
          
          n_bad_inputs_count++; // We count these now.
        }
  }
  else if ( info[1] == "check-for-consensus-during-phasing") {
    // Here we can toggle whether we perform ongoing consensus phase checks throughout phase determination jobs.
    SettingsToUpdate.check_for_consensus_phases = true; // Default is true anyway.
    if (info.size()>2){SettingsToUpdate.check_for_consensus_phases = convert_string_into_bool(info[2]);}
  }

  else if ( info[1] == "keep-centric-phases") {
    // This one takes the format keep-centric-phases bool float bool, where the first bool toggles it on, the float dictates the % of centric reflections to use, and the second bool activates a control, where acentrics are used instead.
    SettingsToUpdate.keep_known_centric_phases = true;
    if (info.size()>2)
    {
      SettingsToUpdate.keep_known_centric_phases = convert_string_into_bool(info[2]);
    }
    if (info.size() > 3)
    {
      // Alter the percent of kept reflections.
      SettingsToUpdate.centric_random_percent = clipper::String(info[3]).f64();
    }
    if (info.size() > 4 )
    {
      SettingsToUpdate.control_for_phase_copying = convert_string_into_bool(info[4]);
    }
    if (info.size() > 5)
    {
      SettingsToUpdate.n_fix_centric_phases = clipper::String(info[5]).i();
    }
  }
  
  else if ( info[1] == "fix-trusted-phases")
  {
    if (info.size() > 2)
    {
      SettingsToUpdate.fix_trusted_phases = clipper::String(info[2]).i();
    } else {
      outlog << "Warning, fix_trusted_phases expects an integer value, representing the number of iterates for which to fix the known phases." << std::endl;
    }
  }
    
  else if ( info[1] == "free-lunch" ) {
    SettingsToUpdate.free_lunch = true; // Whether or not we desire a free lunch.
    free_lunch = true;
    
    if (info.size() > 2 )
    {
      SettingsToUpdate.free_lunch_scale_factor = clipper::String(info[2]).f64();// The Fourier space volume scale factor for the Free lunch algoriuthm - used to calculate the extent of resolution extension. Defaults to 1.0 (no resolution extension)
    }
    
  }  else if (info[1] == "id-offset") {
    SettingsToUpdate.first_run_id_offset = clipper::String(info[2]).i(); //Applies an offset to the run_id, so that we don't have to overwrite runs within the same job if desired.
  }  else if (info[1] == "success-check-rate") {
    SettingsToUpdate.phase_consistency_interval = clipper::String(info[2]).i(); // How often we check for self-consistency of the phase set.
  }  else if (info[1] == "resolution-gaussian-width-cutoff") {
    SettingsToUpdate.resolution_gaussian_width_cutoff = clipper::String(info[2]).i(); // Specifies the the 1/nth maximum height of the Gaussian apodization function, used to calculate the effective resolution. Defaults to 100.
  } else if ( info[1] == "shannon-rate" ) {
    SettingsToUpdate.shannon_rate = clipper::String(info[2]).f64(); // Controls the gridding of the maps at any effective resolution. 
  } else if ( info[1] == "erase-islands" ) {
    SettingsToUpdate.erase_islands=true;
    SettingsToUpdate.erase_islands_interval = clipper::String(info[2]).i(); // Defines how many iterations elapse between enforcement of mask connectivity (the erasure of islands)
  } else if ( info[1] == "erase-islands-threshold" ) {
    SettingsToUpdate.erase_islands=true;
    SettingsToUpdate.erase_islands_threshold = clipper::String(info[2]).f64(); // Defines the threshold (as a fraction of largest set connected set size) for erasure of islands during mask connectivity enforcement
  } else if ( info[1] == "fill-voids" ) {
    SettingsToUpdate.fill_voids=true;
    SettingsToUpdate.fill_voids_interval = clipper::String(info[2]).i(); // Defines how many iterations elapse between enforcement of mask connectivity (the filling of voids)
  } else if ( info[1] == "fill-voids-threshold" ) {
    SettingsToUpdate.fill_voids=true;
    SettingsToUpdate.fill_voids_threshold = clipper::String(info[2]).f64(); // Defines the threshold (as a fraction of largest set connected set size) for removal of voids during mask connectivity enforcement
  } else if ( info[1] == "envelope-editing-limit" ) {
    SettingsToUpdate.envelope_editing_limit = clipper::String(info[2]).f64(); // Defines an effective resolution limit for envelope editing. At resolutions below this limit, envelope editing will be performed. At resolutions above this limit it will not. By default, there is no limit, and envelope editing will be performed irrespective of resolution
  } else if ( info[1] == "output-mask" ) {
    SettingsToUpdate.output_mask=true; // Output mask at intermediate steps ... currently inactivated. TODO: fix or remove
    SettingsToUpdate.mask_output_interval = clipper::String(info[2]).i(); // Defines how many iterations elapse between output of the mask
  } else if ( info[1] == "output-phases" ) {
    SettingsToUpdate.output_phases=true; // Output phases at intermediate steps ... currently inactivated. TODO: fix or remove
    SettingsToUpdate.phase_output_interval = clipper::String(info[2]).i(); // Defines how many iterations elapse between output of the phases
  } else if ( info[1] == "output-envelopes-and-phases" ) {
    SettingsToUpdate.output_envelopes_and_phases=true; // Output both envelopes and phase sets at the end of a run. By default only envelopes are output during envelope determination, and phase sets are output during phase determination
  } else if ( info[1] == "n-apodization-steps" ) {
    SettingsToUpdate.apodize=true;
    SettingsToUpdate.n_apodization_steps = clipper::String(info[2]).i(); // Number of apodization steps. Apodization is systematically changed at each step, smoothly transitioning between an initial and final state. If there is no apodization specified at the last step, this is run without any apodization.
  } else if ( info[1] == "n-first-apodization-offset" ) {
    SettingsToUpdate.apodize=true;
    SettingsToUpdate.n_apodization_first_iterate = clipper::String(info[2]).i(); // The iterate at which apodization starts changing
  } else if ( info[1] == "n-apodization-last-iterate" ) {
    SettingsToUpdate.apodize=true;
    SettingsToUpdate.n_apodization_last_iterate = clipper::String(info[2]).i();// The iterate at which apodization stops changing
  } else if ( info[1] == "initial-apodization-state" ) {
    SettingsToUpdate.apodize=true;
    SettingsToUpdate.initial_apodization_sigma = clipper::String(info[2]).f64(); // Specifies the Initial apodization state, in terms of the standard deviation of the Gaussian function
  } else if ( info[1] == "final-apodization-state" ) {
    SettingsToUpdate.upper_apodization_limit = true;
    SettingsToUpdate.final_apodization_sigma = clipper::String(info[2]).f64(); // Specifies the final apodization state, in terms of the standard deviation of the Gaussian function
  } else if ( info[1] == "gradient-matching" ) {
    SettingsToUpdate.gradient_matching = true; // Perform gradient matching
    SettingsToUpdate.gradient_matching_interval = clipper::String(info[2]).i(); // Defines how often we perform gradient matching
  } else if ( info[1] == "histogram-matching" ) {
    SettingsToUpdate.histogram_matching = true; // Perform histogram matching
    SettingsToUpdate.histogram_matching_interval = clipper::String(info[2]).i(); // Defines how often we perform histogram matching
  } else if ( info[1] == "no-apodize-reference" ) {
    SettingsToUpdate.apodize_reference_data = false; // When peforming histogram matching, the reference data will not be apodized (otherwise the apodization scheme applied to the target is applied to the reference)
  } else if ( info[1] == "n-density-bins" ) {
    SettingsToUpdate.n_density_bins = clipper::String(info[2]).i(); // Number of bins used for the electron density histograms
  } else if ( info[1] == "n-gradient-bins" ) {
    SettingsToUpdate.n_gradient_bins = clipper::String(info[2]).i(); // Number of bins used for the gradient magnitude histograms
  } else if ( info[1] == "resolution-bins" ) {
    SettingsToUpdate.res_bins = clipper::String(info[2]).i(); // Number of resolution bins for the computation of various agreement statistics
  } else if ( info[1] == "verbose" ) {
    SettingsToUpdate.verbose = (info[2] == "true"); // long form output useful for detailed analysis of program function
    // Catch for people who don't spell apodize like an American ??
  } else if ( info[1] == "Apodise" ) {
    SettingsToUpdate.apodize = true;
    SettingsToUpdate.n_apodization_steps = clipper::String(info[2]).i();
  }
  else
  {
    // Display a warning message and the associated line of issue.
    outlog << "WARNING ... Unrecognized argument:\t";
    for (int i = 0; i<info.size(); i++)
    {
      outlog << info[i] << " ";
    }
    outlog << std::endl;
    
    n_bad_inputs_count++; // We count these now.
  }
  
  return true;
}

bool exp_manager::read_string_algorithm_info_into_struct(ipa_settings& settingstruct, std::vector<std::string>& info, std::ostream&outlog)
{
  /*
   This handles the complicated except to the line of an algorithm setting, which dictates the number of iterates, the update rule, and the constraints 
   */
  int max_size = info.size(); // For quick reference.
  
  lower_string(info[1]);
  if (info[1] != "algorithm" && info.size() > 2) // Size must be at least 3. Algorithm must be in position 1. 0th ignored.
  {
    outlog << "A line to interpret algorithm was not valid at the first check." << std::endl;
    return false;
  }
  
  ipa_setting_algorithm_container container;
  
  container.n_iterates = stoi(info[2]); // 2nd entry is the number of iterates for this step.
  
  std::vector<IPA_updaterule_command> temp_commands;
  std::vector<int> temp_n_commands;
  std::vector<IPA_constraint_command> temp_constraints;
  std::vector<int> temp_n_constraints;

  bool found_commands = false;     // Set to true, if we find update rules.
  bool found_constraints = false;  // Set to true, if we find constraints.

  for (int i = 3; i < max_size-1; i ++ ) // Iterate through remaining duplexes (we add one WITHIN the loop to do duplexes, incase we see some "other" input, we can ignore)
  {
    // Parse algorithm UpdateRule type:
    if (info[i+1] == "DifMap") 
    {
      temp_commands.push_back(DifMap);
      temp_n_commands.push_back(stoi(info[i])); // Number comes first.
      found_commands = true;
    }
    else if (info[i+1] == "ErrRed") 
    {
      temp_commands.push_back(ErrRed);
      temp_n_commands.push_back(stoi(info[i])); // Number comes first.
      found_commands = true;
    }
    else if (info[i+1] == "ReReRe") 
    {
      temp_commands.push_back(ReReRe);
      temp_n_commands.push_back(stoi(info[i])); // Number comes first.
      found_commands = true;
    } 
    else if (info[i+1] == "RevRRR") 
    {
      temp_commands.push_back(RevRRR);
      temp_n_commands.push_back(stoi(info[i])); // Number comes first.
      found_commands = true;
    } 
    // Parse algorithm Constraints type:
    else if (info[i+1] == "pHistogramMatch_sFlatten")
    {
      temp_constraints.push_back(pHistogramMatch_sFlatten);
      temp_n_constraints.push_back(stoi(info[i])); // Number comes first.
      found_constraints = true;
    }
    else if (info[i+1] == "pHistogramMatch_sSmooth")
    {
      temp_constraints.push_back(pHistogramMatch_sSmooth);
      temp_n_constraints.push_back(stoi(info[i])); // Number comes first.
      found_constraints = true;
    }

    else if (info[i+1] == "psSkeletonize") // This doesn't function yet, but it can be parsed (for debugging)
    {
      temp_constraints.push_back(psSkeletonize);
      temp_n_constraints.push_back(stoi(info[i])); 
      found_constraints = true;
    }
    else if (info[i+1] == "pGlobicize_pHistogramMatch_sFlatten")
    {
      temp_constraints.push_back(pGlobicize_pHistogramMatch_sFlatten);
      temp_n_constraints.push_back(stoi(info[i])); 
      found_constraints = true;
    }

    // Else we found something weird.
    else
    {
      if (info[i].length() != 0) {
        outlog << " Unrecognised Algorithm: " << info[i+1] << std::endl;
        n_bad_inputs_count++;
        break; // End interpretation, bad input.
      } else break; // This just ignores empty lines.
    }
    
    i++; // Increment an additional integer if we found something, as that's useful.
    // Parse nth times to perform aglorithm
    //temp_n_commands.push_back(stoi(info[i+1]));
  }

  // If ONLY Iterates were entered, then we can fill the rest of the input with the defaults, i.e. Difference Map algorithm, with constraints Histogram matching and solvent flatenning, every iterate ...
  if (found_commands == false)
  {
    // Add default commands for this run.
    temp_commands.push_back(DifMap);
    temp_n_commands.push_back(1);
  }

  if (found_constraints == false)
  {
    // Add default constraint for this run.
    temp_constraints.push_back(pHistogramMatch_sFlatten);
    temp_n_constraints.push_back(1);
  }

  // All done, copy the appropriate stuff to the right places.
  
  // Add the vectors to the container
  container.commands = temp_commands;
  container.n_commands = temp_n_commands;
  container.constraints = temp_constraints;
  container.n_constraints = temp_n_constraints;
  
  // Add the algorithms steps to the experiment regime.
  settingstruct.UpdateRuleRegime.push_back(container);
  
  return true;
}

bool exp_manager::read_string_beta_params_info_into_struct(std::vector<beta_param>& writeparams,std::vector<std::string>& lines,std::ostream& outlog)
{
  if (lines.size() == 9) // MUST be size 9, error otherwise.
  {
    beta_param temp_param = {stoi(lines[2]),stod(lines[3]),stod(lines[4]),stof(lines[5]),stof(lines[6]),stof(lines[7]),stof(lines[8])};
    writeparams.push_back(temp_param);
    return true; // Success
  }
  else
  {
    outlog << "There was an error reading a beta-param value, wrong number of entries." << std::endl;
    return false; // Failure
  }
}

bool exp_manager::read_string_filter_params_info_into_struct(std::vector<filter_param>& writeparams, std::vector<std::string>& lines, std::ostream& outlog)
{
  if (lines.size() == 4) // MUST be size 4, error otherwise.
  {
    filter_param temp_param = {stoi(lines[2]),stod(lines[3])};
    writeparams.push_back(temp_param);
    return true; // Success
  }
  else
  {
    outlog << "There was an error reading a filter-param value, wrong number of entries." << std::endl;
    return false; // Failure
  }
}

bool exp_manager::read_string_solvent_params_info_into_struct(std::vector<solvent_fraction_param>& writeparams, std::vector<std::string>& lines, std::ostream& outlog)
{
  if (lines.size() == 4) // MUST be size 4, error otherwise.
  {
    solvent_fraction_param temp_param = {stoi(lines[2]),stod(lines[3])};
    writeparams.push_back(temp_param);
    return true; // Success
  }
  else
  {
    outlog << "There was an error reading a solvent-fraction-param value, wrong number of entries." << std::endl;
    return false; // Failure
  }
}

void exp_manager::parse_string_into_string_array(std::vector<std::string>& stringArray, const std::string& input, const std::vector<char>& delims)
{
  //std::cout << "Parsing: \"" << input << "\"" << std::endl;
  /*
    This function takes a string and breaks off the initial segment based on delimiters and passes on the rest for processing recursively.
    Its a bit crude but it works, and is only used to parse our parameter file. It also allows us to parse by multiple delimiters now. 
   */
  //stringArray.clear(); // This should be done outside of the function. As it is now recursive.
  int first_delim = -1;
  for (int c = 0 ; c < delims.size(); c++ )
  {
    int new_delim = input.find_first_of(delims[c]);
    if ((new_delim < first_delim && new_delim>=0) || (first_delim == -1 && new_delim >= 0) )
    {
      first_delim = new_delim;
      //std::cout << "Delim " << c << " Found = " << first_delim << std::endl;
    }
  }

  if (first_delim > 0)
   // We found a delimiter! Append the first part, and send the rest to recursion land.
  {
    std::string word = input.substr(0,first_delim);
    if (word.size() > 0) {stringArray.push_back(word);} // Add the word if its valid, and exists.
    if (input.size() - (first_delim + 1) > 0)
    { // Recurse the leftovers, if they're valid and exist.
      std::string leftovers = input.substr(first_delim+1);
      parse_string_into_string_array(stringArray, leftovers, delims);
    }
  }
  else if (first_delim == 0 && input.size() > 1) // The delimiter is first, this happens if theres a pair of deliminators.
  {
    std::string word = input.substr(1);
    parse_string_into_string_array(stringArray, word, delims);
  }
  else // No more delimitres found, we can stop here, ensure whats left is something and chuck it in.
  {
    if (input.size() > 0 && first_delim != 0)
    {
      //std::cout << "Saving: \"" << input << "\"" << std::endl;
      stringArray.push_back(input);
    }
  }
}

void exp_manager::remove_comments_from_string(std::string& line)
{
  /*
   Simple function takes a string reference, and will remove anything right of a comment specifier, i.e. a // or # characters.
   */
  
  // Identify position of comment specifiers:
  int comment_t1 = line.find("//");
  int comment_t2 = line.find("#");
  int commentposition = -1;
  
  // Pick the left most specifier:
  if (comment_t1 != -1 && comment_t2 != -1)
  {
    commentposition = std::min(comment_t1, comment_t2);
  }
  else
  {
    commentposition = std::max(comment_t1, comment_t2);
  }
  // If the line is not empty, and a comment specifier was found, remove the comment from the string.
  if (line.length() > 0 && commentposition != -1)
  {
    line.erase(line.begin() + commentposition,line.end());
  }
}

void exp_manager::lower_string(std::string& string)
{
  // Enforce lowercase.
  for (int i = 0; i < string.length(); i++)
  {
    string[i] = std::tolower(string[i]);
  }
}

std::string exp_manager::return_algorithm_name_from_enum(IPA_updaterule_command myalg)
{
  /*
  THIS IS NOW REDUNDANT AND PRESENT IN THE IPA FUNCTIONS LIBRARY.

   Simple helper function which provides
   a string of the full algorithm name
   from the the codes enum for logging
   purposes.
   */
  switch(myalg)
  {
    case IPA_updaterule_command::DifMap : return "Difference Map";
    case IPA_updaterule_command::ErrRed : return "Error Reduction";
    case IPA_updaterule_command::ReReRe : return "Relax Reflect Relax";
    case IPA_updaterule_command::RevRRR : return "Reverse Relax Reflect Relax";
  }
  
  return "Undefined";
}

void exp_manager::generate_default_experiment_params_file(std::string paramname)
{
  // Generates a experiment.params text file containing almost all the necessary information to run an experiment
  // The settings here are those recommended for optimal performance of the algorithm 
  // As currently configured this version automatically performs the procedure described in Kingston and Millane
  
  // This function is executed when the program is run but cannot find a experiment parameter file.
  // The program then generates a parameter file using this function and stops. 
  // The user is prompting to edit this file and run again.
  
  // TODO ... Check that the output here is consistent with default values for all simple variables
  
  ipa_settings defaults; // Create a new settings object for automating default values
  std::ofstream savefile; // Create new file out stream.
  savefile.open (paramname);
  if (!savefile.good())
  {   // Simple sanity check.
    std::cerr << "Bad " << paramname << " parameter file generation!";
    return;
  }
  // Print some tutorial information first.
  savefile << "# *------------------------------------------------*" << std::endl;
  savefile << "# |                                                |" << std::endl;
  savefile << "# | Autogenerated experiment.params example file   |" << std::endl;
  savefile << "# | Version: " << xmake_string(release_version) <<  "                                |" << std::endl;
  savefile << "# | Comments after // or # will be ignored,        |" << std::endl;
  savefile << "# | Only jobs outlined under -JOBS are performed,  |" << std::endl;
  savefile << "# | Settings begin with the - (hyphen) character,  |" << std::endl;
  savefile << "# | Spaces or Tab characters separate inputs.      |" << std::endl;
  savefile << "# |                                                |" << std::endl;
  savefile << "# *------------------------------------------------*" << std::endl;

  savefile << std::endl;

  savefile << "# JOB LIST OPTIONS: INIT ENVELOPE_DETERMINATION ENVELOPE_CONSENSUS PHASE_DETERMINATION PHASE_CONSENSUS" << std::endl;
  savefile << "-JOBS INIT ENVELOPE_DETERMINATION ENVELOPE_CONSENSUS PHASE_DETERMINATION # Delete or add desired jobs." << std::endl;
  savefile << std::endl;
  savefile << "# *-------------- CORE SETTINGS -------------------*" << std::endl;
  savefile << "# | These are tagged by the -core keyword          |" << std::endl;
  savefile << "# *------------------------------------------------*" << std::endl;
  savefile << std::endl;
  savefile << "# *-----*        COMPULSORY SETTINGS:        *-----*" << std::endl;
  savefile << "-core input-data-mtz-filename ./inputs/XXX.mtz # Path to the mtz file holding the diffraction data." << std::endl;
  savefile << "-core input-data-fsigf-column-labels /*/*/[FP,SIGFP] # The labels for the Structure Factor Amplitudes and their standard deviations in the mtz file" << std::endl;
  savefile << "-core input-data-solvent-fraction 0.XXX # The solvent fraction - must be between 0 and 1" << std::endl;
  savefile << std::endl;
  savefile << "# *-----*       RECOMMENDED SETTINGS:        *-----*" << std::endl;
  savefile << "-core job-name jobdefault # a customized name for the job "  << std::endl;
  savefile << "-core threads 10 # The maximum number of threads/cores to use for the calculations." << std::endl;
  savefile << std::endl;
  savefile << "# *-----*       OPTIONAL SETTINGS:           *-----*" << std::endl;
//  savefile << "#-core input-data-phiw-column-labels /*/*/[PHI,WEIGHTS] # Optionally use a starting phase set." << std::endl;
  savefile << "-core verbose false # Amount of information present in the log-files. Verbose true results in *very* verbose output. " << std::endl;
  savefile << "-core statistics-interval 1 # Controls how frequently the statistics reporting the agreement with the constraints are calculated." << std::endl;
  savefile << "-core infinity-mode false # If true, keep trying to find a solution to the phase problem until the universe suffers heat death." << std::endl;
  savefile << "#-core truncate-data 3.0 # Truncate the high-resolution data on input." << std::endl;
  savefile << std::endl;
  savefile << "# *-----*   Compositional Settings:    *-----* " << std::endl;
  savefile << "#-core protein-sequence ./inputs/XXX.txt # Name of a file containing the sequences of any protein chains present in the asymmetric unit,  Use single letter code, white spaces are ignored. There are equivalent keywords for dna or rna." << std::endl;
  savefile << "#-core protein-n-chains 1 # the number of distinct chains present in the sequence file. There are equivalent keywords for dna and rna." << std::endl;
  savefile << "#-core protein-copies-in-asu 1 # the number of copies of the chains present in the asu. There are equivalent keywords for dna and rna." << std::endl;
  savefile << std::endl;
  savefile << "# *------------------------------------------------*" << std::endl;
  savefile << "# |               ADVANCED SETTINGS                |" << std::endl;
  savefile << "# |      Please onnly change these if you          |" << std::endl;
  savefile << "# |      understand what you are doing.            |" << std::endl;
  savefile << "# *------------------------------------------------*" << std::endl;
  savefile << std::endl;
  savefile << "# *------- ENVELOPE DETERMINATION SETTINGS --------*" << std::endl;
  savefile << "# |    These are tagged by the -envelope keyword   |" << std::endl;
  savefile << "# *------------------------------------------------*" << std::endl;
//  savefile << "-envelope algorithm 1500 1475 DifMap 25 ErrRed 1 pHistogramMatch_sFlatten # A block of iterations, followed by the update rules to be applied and the constraints to be used. These will be repeated in the sequence given until the block ends" << std::endl;
  savefile << "-envelope algorithm 1525 1500 DifMap 25 ErrRed 1 pHistogramMatch_sFlatten # A block of iterations, followed by the update rules to be applied and the constraints to be used. These will be repeated in the sequence given until the block ends" << std::endl;
  savefile << "-envelope mask-from variance # Controls which local map is used for envelope generation. Accepts mean, variance amd mean-and-variance. In the latter case a number between 0 and 1 gives the weighting." << std::endl;  
  savefile << "-envelope mask-redefinition-interval 1 # The number of iterations that elapse before we redefine the protein envelope" << std::endl;
  savefile << "-envelope n-density-bins 250 # The # of bins used for remapping of density value histograms" << std::endl;
  savefile << "-envelope resolution-bins 40 # The # of bins used for reporting resolution-dependent statistics" << std::endl;
  savefile << "-envelope prot-density 0.44 # The expected mean electron density in protein region (electrons/Angstrom^3)" << std::endl;
  savefile << "-envelope solvent-density 0.38 # The expected mean electron density in solvent region (electrons/Angstrom^3)" << std::endl;
//  savefile << "-envelope n-runs 50 # The total number of runs to be executed" << std::endl;
  savefile << "-envelope n-runs 60 # The total number of runs to be executed" << std::endl;
  savefile << "-envelope initial-apodization-state 0.09079 # The initial apodization state - specified as the standard deviation of a Gaussian function in Fourier space (1/Å)" << std::endl;
  savefile << "-envelope remove-low-res-data 25 # Cutoff for removal of ultra-low resolution data. Data below this limit will be treated as missing throughout" << std::endl;
  savefile << "-envelope probability-threshold 0.995 # Probability cutoff used during Fourier space projection. When extimated intensities for the missing data cross this threshold, operations are performed to reduce their magnitude" << std::endl;
  savefile << "-envelope replace-largei-with-attenuatedi 0.6 # During Fourier space projection, missing data with improbably large reconstructed intensities will be scaled using the supplied multiplier " << std::endl;
//  savefile << "-envelope beta-function 1 # Controls the parametric function via which parameter beta is systematically varied.  1 = square pulse train. 2 = Flattened sinusoid. " << std::endl;
  savefile << "-envelope beta-function 2 # Controls the parametric function via which parameter beta is systematically varied.  1 = square pulse train. 2 = Flattened sinusoid. " << std::endl;
  savefile << "-envelope linear-period-change true # When the function varying beta is chirped, period will be shift linearly with iterate if true, frequency will shift linearly with iterate if false" << std::endl;
  savefile << "-envelope resolution-gaussian-width-cutoff 120 # Specifies the fraction of maximum height that will be used to compute the effective resolution of the data, when apodizing with a Gaussian function " << std::endl;
  savefile << "-envelope shannon-rate 1.32 # Controls the gridding of the maps at any effective resolution" << std::endl;
  //savefile << "//-envelope id-offset 0 // Offset the generated logs and masks by this offset number, useful if you want to generated masks n to 20 after already generating masks 0 to n and cut the run short." << std::endl;
  //savefile << "//-envelope fix-trusted-phases -1 // Number of iterations to fix the trusted phase set" << std::endl;
  savefile << std::endl;
  savefile << "# Parameters controlling the value of beta, used in the Difference Map and RRR algorithms" << std::endl;
  savefile << "#" << std::endl;
  savefile << "#                   Iterate Midline Amplitude Phase  Period  Func2Param  Func1Param" << std::endl;
  savefile << "-envelope beta-param      0   0.75    0.120    0.50     56       4            1" << std::endl;
  savefile << "-envelope beta-param    699   0.75    0.120    0.50     80       4            1" << std::endl;
  savefile << "-envelope beta-param   1399   0.75    0.025    0.50    144       2            1" << std::endl;
  savefile << "-envelope beta-param   1499   0.75    0.025    0.50    144       2            1" << std::endl;
//  savefile << "-envelope beta-param 0       0.75    0.03      0.50   2       1           0.5" << std::endl;
//  savefile << "-envelope beta-param 1474    0.75    0.03      0.50   2       1           0.5" << std::endl;
//  savefile << "-envelope beta-param 1475    0.0     0.00      0.50   1       1           1" << std::endl;
//  savefile << "-envelope beta-param 1499    0.0     0.00      0.50   1       1           1" << std::endl;
  savefile << "#" << std::endl;
  savefile << "# Parameters controlling the the filter radius, used for generating the molecular envelope" << std::endl;
  savefile << "# Specify iterate and radius, the program will linearly interpolate between the specified states" << std::endl;
  savefile << "#" << std::endl;
  savefile << "#                     Iterate Radius" << std::endl;
  savefile << "-envelope filter-param      0 13.25 " << std::endl;     
  savefile << "-envelope filter-param    699 12.75 " << std::endl;     
  savefile << "-envelope filter-param   1399 12.25  " << std::endl;    
  savefile << "-envelope filter-param   1499 12.25 " << std::endl;     
//  savefile << "-envelope filter-param 0       10.8  " << std::endl;
//  savefile << "-envelope filter-param 999      8.0  " << std::endl;
//  savefile << "-envelope filter-param 1499     8.0  " << std::endl;
  savefile << "#" << std::endl;
  savefile << "# Parameters controlling the solvent fraction multiplier" << std::endl;
  savefile << "# Specify iterate and multiplier, the program will linearly interpolate between the specified states" << std::endl;
  savefile << "#" << std::endl;
  savefile << "#                      Iterate    Multiplier" << std::endl;
  savefile << "-envelope solvent-param 0           1.04" << std::endl;    
  savefile << "-envelope solvent-param 699         1.02" << std::endl;     
  savefile << "-envelope solvent-param 1399        1.00" << std::endl;     
  savefile << "-envelope solvent-param 1499        1.00" << std::endl;     
//  savefile << "-envelope solvent-param 0          1.0" << std::endl;
//  savefile << "-envelope solvent-param 1499       1.0" << std::endl;
  savefile << "" << std::endl;

  savefile << "# *-------- ENVELOPE CONSENSUS SETTINGS -----------*" << std::endl;
  savefile << "# | These are tagged by the -mskconsensus keyword  |" << std::endl;
  savefile << "# *------------------------------------------------*" << std::endl;
//  savefile << "-mskconsensus min-points 5 # DBSCAN parameter MinPoints" << std::endl;
  savefile << "-mskconsensus min-points 4 # DBSCAN parameter MinPoints" << std::endl;
  savefile << "-mskconsensus epsilon-as-percentile 16.67 # DBSCAN threshold value ε, specified as a percentile of the k-distance distribution" << std::endl;
  savefile << "-mskconsensus erase-islands-threshold 0.05 # defines a threshold for removal of small connected sets during mask connectivity enforcement. Specified as a fraction of the largest set connected set size " << std::endl;
  savefile << "-mskconsensus fill-voids-threshold 0.05 # defines a threshold for removal of small connected sets during mask connectivity enforcement. Specified as a fraction of the largest set connected set size " << std::endl;
  savefile << "-mskconsensus distance-measure 2 # the distance measure used for clustering 1 = 1-|cc|, 2 = √(1-cc^2) " << std::endl;
  savefile << "-mskconsensus allow-inversion # inversion of envelopes will be allowed (if permitted by the space group) when clustering and averaging " << std::endl;
  savefile << "-mskconsensus verbose true # Turn on verbose output for the mask consensus step only" << std::endl;
  savefile << "" << std::endl;
  savefile << "# *-------- PHASE DETERMINATION SETTINGS ----------*" << std::endl;
  savefile << "# |    These are tagged by the -phasing keyword    |" << std::endl;
  savefile << "# *------------------------------------------------*" << std::endl;
  savefile << "#-phasing mskin-target-working /envelope_consensus/jobname/XXXX.ccp4 # Set this to add a specific envelope to the list that will be used during phase determination." << std::endl;
  savefile << "-phasing algorithm 7200 1 DifMap  1 pHistogramMatch_sFlatten # A block of iterations, followed by the update rules to be applied and constraints to be used. These will be repeated in the sequence given until the block ends" << std::endl;
  savefile << "-phasing algorithm 900 200 DifMap 25 ErrRed 1 pHistogramMatch_sFlatten # A block of iterations, followed by the update rules to be applied and the constraints to be used. These will be repeated in the sequence given until the block ends" << std::endl;
  savefile << "-phasing mask-from variance" << std::endl;
  savefile << "-phasing mask-redefinition-interval 1" << std::endl;
  savefile << "-phasing n-iterations-fix-mask 10 # the number of iterations for which the mask will be fixed" << std::endl; 
  savefile << "-phasing n-apodization-steps 31 # the number of apodizations steps " << std::endl;
  savefile << "-phasing n-first-apodization-offset 0 # the number of iterations before we begin adjusting the apodization " << std::endl;
  savefile << "-phasing n-apodization-last-iterate 7440 # the end of the final apodization step" << std::endl; // Due to how division works, and 31 apo steps works, this is actually finishing with the 31st apo step (30th counting from 0) on iterate 7200
  savefile << "-phasing initial-apodization-state 0.157130378 # The initial apodization state - specified as the standard deviation of a Gaussian function in Fourier space (1/Å)" << std::endl;
  savefile << "-phasing remove-low-res-data 25" << std::endl;
  savefile << "-phasing n-density-bins 250" << std::endl;
  savefile << "-phasing resolution-bins 40" << std::endl;
  savefile << "-phasing prot-density 0.44" << std::endl;
  savefile << "-phasing solvent-density 0.38" << std::endl;
  savefile << "-phasing n-runs 20" << std::endl;
  savefile << "-phasing probability-threshold 0.99995 # Probability cutoff used during Fourier space projection. When estimated intensities for the missing data cross this threshold, action is taken to reduce their magnitude" << std::endl;
  savefile << "-phasing replace-largei-with-expected-value 1.0 # During Fourier space projection, missing data with improbably large intensities will be replaced with E[|F|]^2, scaled using the supplied multiplier " << std::endl;
  savefile << "-phasing beta-function 1 # Controls the parametric function via which parameter beta is systematically varied.  1 = square pulse train. 2 = Flattened sinusoid. " << std::endl;
  savefile << "-phasing linear-period-change true # When the function varying beta is chirped, period will be shift linearly with iterate if true, frequency will shift linearly with iterate if false" << std::endl;
  savefile << "-phasing resolution-gaussian-width-cutoff 120 # Specifies the fraction of maximum height that will be used to compute the effective resolution of the data, when apodizing with a Gaussian function " << std::endl;
  savefile << "-phasing shannon-rate 1.32; # Controls the gridding of the maps at any effective resolution" << std::endl;
  savefile << "#-phasing free-lunch 1.0 # Data will be extended past the measured resolution limit, based on the supplied Fourier space volume multiplier (1/Å^3)" << std::endl;
  //savefile << "//-phasing id-offset 0 // Offset the generated logs and masks by this offset number, useful if you want to generated masks n to 20 after already generating masks 0 to n and had cut the run short." << std::endl;
  //savefile << "//-phasing fix-trusted-phases -1" << std::endl;
  savefile << std::endl;
  savefile << "# Parameters controlling the value of beta, used in the Difference Map and RRR algorithms" << std::endl;
  savefile << "#" << std::endl;
  savefile << "#                   Iterate Midline Amplitude Phase  Period  Func2Param  Func1Param" << std::endl;
  // This variant gets the flipping of beta in the sequence used for the IUCrJ paper
  savefile << "-phasing beta-param  0       0.7375  0.0625    0.250    120   1           0.5" << std::endl;
  savefile << "-phasing beta-param  7199    0.7375  0.0625    0.250    120   1           0.5" << std::endl;
  savefile << "-phasing beta-param  7200    0.1     0.6500    0.750    200   1           0.5" << std::endl;
  savefile << "-phasing beta-param  7399    0.1     0.6500    0.750    200   1           0.5" << std::endl;
  savefile << "-phasing beta-param  7400    0.0     0.0000    0.750    200   1           1  " << std::endl;
  savefile << "-phasing beta-param  7424    0.0     0.0000    0.750    200   1           1  " << std::endl;
  savefile << "-phasing beta-param  7425    0.1     0.6500    0.625    200   1           0.5" << std::endl;
  savefile << "-phasing beta-param  7624    0.1     0.6500    0.625    200   1           0.5" << std::endl;
  savefile << "-phasing beta-param  7625    0.0     0.0000    0.625    200   1           1  " << std::endl;
  savefile << "-phasing beta-param  7649    0.0     0.0000    0.625    200   1           1  " << std::endl;
  savefile << "-phasing beta-param  7650    0.1     0.6500    0.500    200   1           0.5" << std::endl;
  savefile << "-phasing beta-param  7849    0.1     0.6500    0.500    200   1           0.5" << std::endl;
  savefile << "-phasing beta-param  7850    0.0     0.0000    0.500    200   1           1  " << std::endl;
  savefile << "-phasing beta-param  7874    0.0     0.0000    0.500    200   1           1  " << std::endl;
  savefile << "-phasing beta-param  7875    0.1     0.6500    0.375    200   1           0.5" << std::endl;
  savefile << "-phasing beta-param  8074    0.1     0.6500    0.375    200   1           0.5" << std::endl;
  savefile << "-phasing beta-param  8075    0.0     0.0000    0.375    200   1           1  " << std::endl;
  savefile << "-phasing beta-param  8099    0.0     0.0000    0.375    200   1           1  " << std::endl;
  savefile << "" << std::endl;
  savefile << "# Parameters controlling the the filter radius, used for generating the molecular envelope" << std::endl;
  savefile << "# Specify iterate and radius, the program will linearly interpolate between the specified states" << std::endl;
  savefile << "#" << std::endl;
  savefile << "#                     Iterate Radius" << std::endl;
  savefile << "-phasing filter-param  0       8.0" << std::endl;
  savefile << "-phasing filter-param  8099    8.0" << std::endl;
  savefile << "" << std::endl;
  savefile << "# Parameters controlling the solvent fraction multiplier" << std::endl;
  savefile << "# Specify iterate and multiplier, the program will linearly interpolate between the specified states" << std::endl;
  savefile << "#" << std::endl;
  savefile << "#                      Iterate    Multiplier" << std::endl;
  savefile << "-phasing solvent-param  0          1.0" << std::endl;
  savefile << "-phasing solvent-param  8099       1.0" << std::endl;
  savefile << "-phasing check-for-consensus-during-phasing true # Whether or not to check for consensus from finished runs in between threads." << std::endl;
  savefile << "" << std::endl;
  savefile << "# *----------- PHASE CONSENSUS SETTINGS ------------*" << std::endl;
  savefile << "# | These are tagged by the -phsconsensus keyword   |" << std::endl;
  savefile << "# *------------------------------------------------*" << std::endl;
  savefile << "-phsconsensus epsilon 40.0 # DBSCAN threshold value ε, specified in terms of mean absolute phase difference" << std::endl; // Based on test cases, ε of 50.0° may generate the very occasional false positive. ε of 40.0° generates no false positives, but occasional false negatives ...
  savefile << "-phsconsensus min-points 2 # DBSCAN parameter MinPoints" << std::endl;
  savefile << "-phsconsensus verbose true # Turn on verbose output for the phase consensus step only" << std::endl;
  savefile << std::endl;
  savefile << std::endl;

  savefile << "# *-----*        DEBUG and TEST SETTINGS:        *-----*" << std::endl;
  savefile << "#-core really-random false # for TESTING and DEBUGGING: If false, generates the same random numbers on each invocation of the program " << std::endl;
  savefile << "#-core known-solution-mask-filename ./inputs/XXXX.ccp4 # for TESTING and DEBUGGING: A file containing the known mask" << std::endl;
  savefile << "#-core known-solution-phase-filename ./inputs/XXXX.mtz # for TESTING and DEBUGGING: A file containing known phases " << std::endl;
  savefile << "#-core known-solution-fp-column-labels */*/[FC,PHIC] # for TESTING and DEBUGGING: The labels for the known Fourier amplitudes and phases in the mtz file. Amplitiudes are ignored"  << std::endl;
  savefile << "#-phsconsensus shift-to-known-origin # for TESTING and DEBUGGING: Consensus phase sets will be shifted to the same origin as a known phase set, and written out "  << std::endl;
  savefile << "#-mskconsensus shift-to-known-origin # for TESTING and DEBUGGING: Consensus envelopes will be shifted to the same origin as a known envelope, and written out"  << std::endl;

  savefile << std::endl;
  savefile.close();
  
  std::cout << "A default parameter file (Filename: " << paramname << ") has been automatically generated in the directory where program IPA was run.\n" << std::endl;
  std::cout << "Please update the following -core settings with compulsory user input before rerunning program IPA:\n" << "-core input-data-mtz-filename filename.mtz\n" << "-core input-fsigf-column-labels /*/*/[FP,SIGFP]\n" << "-core solvent-fraction 0.XX\n" << std::endl;
}

void exp_manager::output_settings_summary(std::ostream& outlog)
{
  /*
   Call this for a comprehensive output of all pertinent settings.
   */

  std::string boldme = "";
  if (&outlog == &std::cout) boldme = cBold;
  
  util::print_section_line(cl,outlog);
  outlog << util::parse_line(cl, "      SUMMARY OF INPUT SETTINGS:      ", 0, boldme);
  util::print_section_line(cl,outlog);
  util::print_empty_line(cl,outlog);

  if (custom_job_id)
  {
    outlog << util::parse_line(cl, " Settings for Job \"" + job_name + "\"");
    outlog << util::parse_line(cl, " Outputs for this job will be stored in a /" + job_name + "/ folder hierarchy.");
  }
  
  
  if (core_solvent_fraction_provided && envelope_settings.fraction_solvent != 0.0 && phasing_settings.fraction_solvent != 0.0)
  {
    outlog << util::parse_line(cl, " Solvent Fraction provided: " + std::to_string(envelope_settings.fraction_solvent));
  }
  else
  {
    // This check should really be done previously, and picked up appropriately with red text.
    outlog << util::parse_line(cl, " Solvent Fraction MUST BE PROVIDED, use -core solvent-fraction 0.X, where 0.X is the fractional solvent", 0, cRed);
    n_bad_inputs_count++;
  }
  
   if (!found_job_list)
  {
    outlog << " No -JOB list was found! Indicate which jobs are to be performed using the -JOBS keyword within the parameter file." << std::endl;
  }
  
  if (envelope_settings.really_random == false)
  {
    outlog << util::parse_line(cl, " Pseudo-random numbers will be used for debugging purposes.");
  }
  
  if (envelope_settings.verbose == true) // Note that verbose logging can now be turned on individually for the envelope and phase consensus jobs ...
  {
    outlog << util::parse_line(cl, " Verbose logging is enabled during envelope and phase determination");
  } else {
    outlog << util::parse_line(cl, " Verbose logging is disabled during envelope and phase determination");
  }
  
  util::print_empty_line(cl,outlog);

  if (JOB_initialisation)
  {
    outlog << util::parse_line(cl, "      INITIALISATION:      ",0,boldme);
    outlog << util::parse_line(cl, " Will create Folder hierarchy for this job.");
    if (!core_b_factor_provided)
    {
      outlog << util::parse_line(cl, " An overall Anisotropic B-Factor will be calculated.");
    }
    if (protein_sequence_is_provided || rna_sequence_is_provided || dna_sequence_is_provided)
    {
      outlog << util::parse_line(cl, " The Follow Compositional information will be used to help evaluate the B-Factor:");
      if (protein_sequence_is_provided)
      {
        outlog << util::parse_line(cl, " Protein composition: (" + std::to_string(protein_sequence_nchains) + " chains, " + std::to_string(protein_sequence_ncopies_in_asu) + " copies in ASU )\n Sequence from filename: " + protein_sequence_filename);
      }
      if (dna_sequence_is_provided)
      {
        outlog << util::parse_line(cl, " DNA composition: (" + std::to_string(dna_sequence_nchains) + " chains, " + std::to_string(dna_sequence_ncopies_in_asu) + " copies in ASU )\n Sequence from filename: " + dna_sequence_filename);
      }
      if (rna_sequence_is_provided)
      {
          outlog << util::parse_line(cl, " RNA composition: (" + std::to_string(rna_sequence_nchains) + " chains, " + std::to_string(rna_sequence_ncopies_in_asu) + " copies in ASU )\n Sequence from filename: " + rna_sequence_filename);
      }
    }
    else
    {
      outlog << util::parse_line(cl, " No compositional information was provided to aide B-factor estimation");
    }
    util::print_empty_line(cl,outlog);
  }
  
  if (JOB_envelope_determination)
  {
    outlog << util::parse_line(cl, "      ENVELOPE DETERMINATION:      ",0,boldme);
    unpack_settings_as_log(outlog, envelope_settings);
    util::print_empty_line(cl,outlog);
  }
  
  if (JOB_envelope_consensus)
  {
    outlog << util::parse_line(cl, "      ENVELOPE CONSENSUS:      ",0,boldme);
    unpack_settings_as_log(outlog, envelope_consensus_mask_settings);
    util::print_empty_line(cl,outlog);
  }
  
  if (JOB_phase_determination)
  {
    outlog << util::parse_line(cl, "      PHASE DETERMINATION:      ",0,boldme);
    unpack_settings_as_log(outlog, phasing_settings);
    util::print_empty_line(cl,outlog);
  }
  
  if (JOB_phase_consensus)
  {
    outlog << util::parse_line(cl, "      PHASE CONSENSUS:      ",0, boldme);
    unpack_settings_as_log(outlog, phase_consensus_phase_settings);
    util::print_empty_line(cl,outlog);
  }
  
  // Display the error count summary:
  util::print_empty_line(cl,outlog);
  std::string counts = std::to_string(n_bad_inputs_count) + " errors were detected when parsing settings.";
  int n_hacks = 0;
  if ( n_bad_inputs_count != 0 ) {
    counts = "\u001b[31m" + counts + "\u001b[0m";  // RED if bad
    n_hacks = 9;
  }
  outlog << util::parse_line(cl, counts, n_hacks);
  util::print_section_line(cl,outlog);
}

bool exp_manager::unpack_settings_as_log(std::ostream& outlog, consensus_mask_settings& settings)
{

  outlog << util::parse_line(cl, " DBSCAN parameter MinPoints is " + std::to_string(settings.min_pts));

  if (settings.epsilon_is_absolute)
  {
    outlog << util::parse_line(cl, " DBSCAN parameter ε is specified in absolute terms", 1);
  }
  else
  {
    outlog << util::parse_line(cl, " DBSCAN parameter ε is specified as a percentile of", 1);
    outlog << util::parse_line(cl, " the k-distance distribution");
  }

  if (envelope_consensus_epsilon.size() > 1) // There is more than 1 epsilon value to use.
  {
    outlog << util::parse_line(cl, std::to_string(envelope_consensus_epsilon.size()) + " epsilon values will be used: ");
    // Make a string list of the values.
    std::string list_of_values = "";
    for (int i = 0; i < envelope_consensus_epsilon.size(); i++)
    {
      list_of_values = list_of_values + std::to_string(envelope_consensus_epsilon[i]);
      if (i == envelope_consensus_epsilon.size() - 1)
      {
        list_of_values += ".";
      }
      else
      {
        list_of_values += ", ";
      }
    }
    outlog << util::parse_line(cl, list_of_values);
  }
  else if (envelope_consensus_epsilon.size() == 1)
  {
    outlog << util::parse_line(cl, " An epsilon value of " + std::to_string(envelope_consensus_epsilon[0]) + " will be used: ");
  }
  else if (envelope_consensus_epsilon.size() <= 0)
  {
    outlog << util::parse_line(cl, " Warning: No epsilon values were provided in the input settings file!", 0, cRed);
    outlog << util::parse_line(cl, " A default value of " + std::to_string(settings.epsilon) + " will be used.", 0, cRed);
  }
  return true;
}

bool exp_manager::unpack_settings_as_log(std::ostream &outlog, consensus_phase_settings &settings)
{
  outlog << util::parse_line(cl, " Consensus phases will be determined based on the files found in phase_determination/" + job_name);
  outlog << util::parse_line(cl, " DBSCAN parameter MinPoints is " + std::to_string(settings.min_pts));
  outlog << util::parse_line(cl, " An epsilon value of " + std::to_string(settings.epsilon_input) + " will be used: ");
  return true;
}

bool exp_manager::unpack_settings_as_log(std::ostream& outlog, ipa_settings& settings)
{
  /*
   Unpack any interesting settings that are derived from the ipa_settings struct to be checked and logged before a run here.
   */
  std::string name;
  if (settings.job_is_envelope_determination) name = "Envelope";
  if (settings.job_is_phase_determination) name = "Phase";

  outlog << util::parse_line(cl," Execute " + std::to_string(settings.n_runs) + " " + name + " determination runs, using " + std::to_string(n_threads) + " threads.");
  
  if (settings.first_run_id_offset != 0)
  {
    outlog << util::parse_line(cl, " Run IDs will start from " + std::to_string(settings.first_run_id_offset));
  }
    // Debuging output.
  outlog << util::parse_line(cl, " " + name + " IPA Update Regime: ");
  
  for (int i = 0; i < settings.UpdateRuleRegime.size(); i++)
  {
    ipa_setting_algorithm_container& container = settings.UpdateRuleRegime[i];
    
    std::string my_string = "  > For " + std::to_string(container.n_iterates) + " iterates, repeat ";

    for (int i = 0; i < container.n_commands.size(); i ++)
    {
      std::string alg_name = util::return_name_from_enum(container.commands[i]);
      
      if (i>0 && i != container.n_commands.size() - 1) {my_string = my_string + ", ";}
      if (i == container.n_commands.size() - 1 && i!=0) {my_string = my_string + ", and ";}
      my_string = my_string + alg_name;
      if (container.n_commands[i] == 1) {
        my_string =  my_string + " for " + std::to_string(container.n_commands[i]) + " iterate";
      }
      if (container.n_commands[i] != 1) {
        my_string =  my_string + " for " + std::to_string(container.n_commands[i]) + " iterates";
      }
    }
    if (container.n_constraints.size() == 1)
    {
      my_string = my_string + ", and repeat the constraint ";
    }

    if (container.n_constraints.size() > 1)
    {
      my_string = my_string + ", and repeat the constraints: ";
    }
    
    for (int i = 0; i < container.n_constraints.size(); i ++)
    {
      std::string alg_name = util::return_name_from_enum(container.constraints[i]);
      
      if (i>0 && i != container.n_constraints.size() - 1) {my_string = my_string + ", ";}
      if (i == container.n_constraints.size() - 1 && i!=0) {my_string = my_string + ", and ";}
      my_string = my_string + alg_name;
      if (container.n_constraints[i] == 1) {
        my_string = my_string + " for " + std::to_string(container.n_constraints[i]) + " iterate";
      }
      if (container.n_constraints[i] != 1) {
        my_string = my_string + " for " + std::to_string(container.n_constraints[i]) + " iterates";
      }
      
    }
    my_string = my_string + ".";
    outlog << util::parse_line(cl,my_string);
  }

  util::print_empty_line(cl,outlog);
  outlog << util::parse_line(cl," The total number of iterates per run is " + std::to_string(settings.n_total_iterates) + ".");
  
  if ( !settings.apodize )
  {
    outlog << util::parse_line(cl,"No apodization will be performed.") << std::endl;
  }
  else
  {
    if (settings.n_apodization_steps == 0)
    {
      outlog << util::parse_line(cl," An unchanging apodizaton function will be employed -  the data will be");
      outlog << util::parse_line(cl," multiplied with a Gaussian of sigma " + std::to_string(settings.initial_apodization_sigma) + " 1/Å.", 1);  // Split up manually to account for the silly angstrom char   
    }
    
    if (settings.n_apodization_steps > 0)
    {
      outlog << util::parse_line(cl, " The apodization function will change in " + std::to_string(settings.n_apodization_steps) + " steps of even length");
      outlog << util::parse_line(cl," The apodization function will start changing at iterate " + std::to_string(settings.n_apodization_first_iterate) + ", and finish changing at iterate " + std::to_string(settings.n_apodization_last_iterate) + ".");
      outlog << util::parse_line(cl," At the initial apodization step, the data will be multiplied with a");
      outlog << util::parse_line(cl," Gaussian of sigma " + std::to_string(settings.initial_apodization_sigma) + " 1/Å.", 1); // Split up manually to account for the silly angstrom char   

      
      if (settings.upper_apodization_limit == false)
      {
        outlog << util::parse_line(cl," At the final apodization step, the unmodified data will be used");
      }
      else
      {
        outlog << util::parse_line(cl," At the final apodization step, the data will will be multiplied");
        outlog << util::parse_line(cl," with a Gaussian of sigma " + std::to_string(settings.final_apodization_sigma) + " Å.", 1); // Split up manually to account for the silly angstrom char
      }
    }
  }

  if (settings.free_lunch_scale_factor != 1.0 && settings.free_lunch)
  {
    outlog  << util::parse_line(cl,"The dataset will be expanded by " + std::to_string((settings.free_lunch_scale_factor - 1) * 100) + "\% of the spherical volume in Fourier space, to a resolution of " + std::to_string(settings.hkl_wilson.resolution().limit())) << " Å" << std::endl;
  }
  
  if (settings.found_beta_params)
  {
    outlog << util::parse_line(cl, " Custom Beta Parameter will be used.");
  }
  if (settings.found_filter_params)
  {
    outlog << util::parse_line(cl, " Custom Filter Radius will be used.");
  }
  if (settings.found_solvent_fraction_multiplier_params)
  {
    outlog << util::parse_line(cl, " Custom Solvent Fraction Multiplier will be used.");
  }
  
  if (settings.compute_mask_from_local_mean == true)
  {
    outlog << util::parse_line(cl, " Mask will be computed using local mean");
  }
  if (settings.compute_mask_from_local_variance == true)
  {
      outlog << util::parse_line(cl, " Mask will be computed using local variance");
  }
  if (settings.compute_mask_from_local_mean_and_variance == true)
  {
        outlog << util::parse_line(cl, " Mask will be computed using local mean and variance");
  }
  
  if (settings.input_working_mask == true)
  {
          outlog << util::parse_line(cl, " An initial working mask has been provided");
  }
  
  return true;
}

bool exp_manager::read_mask_consensus_setting_into_struct(std::vector<std::string>& info, std::ostream& outlog)
{
  if ( info[1] == "epsilon-absolute" ) {
    envelope_consensus_mask_settings.epsilon_is_absolute = true;
    if ( 2 < info.size() ) 
    {
      //envelope_consensus_mask_settings.epsilon = clipper::String(info[2]).f64(); // the threshold distance for clustering in absolute terms
      envelope_consensus_epsilon.push_back(clipper::String(info[2]).f64()); // Append this epsilon value to our eps container.
    }
    return true;
  }
  
  if ( info[1] == "epsilon-as-percentile" ) {
    envelope_consensus_mask_settings.epsilon_is_absolute = false;
    if ( 2 < info.size() ) 
    {
      //envelope_consensus_mask_settings.epsilon = clipper::String(info[2]).f64(); // the threshold distance for clustering in terms of a percentile of the k-distance distribution
      envelope_consensus_epsilon.push_back(clipper::String(info[2]).f64()); // Append this epsilon value to our eps container.
    }
    return true;
  }

  
  if ( info[1] == "min-points" ) {
    if ( 2 < info.size() ) envelope_consensus_mask_settings.min_pts = clipper::String(info[2]).i(); // defines the minimum cluster size
    return true;
  }
  
  if ( info[1] == "erase-islands-threshold" ) {
    if ( 2 < info.size() ) envelope_consensus_mask_settings.erase_islands_threshold = clipper::String(info[2]).f64(); // defines the threshold (as a fraction of largest set connected set size) for removal of small connected sets during mask connectivity enforcement
    return true;
  }
  
  if ( info[1] == "fill-voids-threshold" ) {
    if ( 2 < info.size() ) envelope_consensus_mask_settings.fill_voids_threshold = clipper::String(info[2]).f64(); // defines the threshold (as a fraction of largest set connected set size) for removal of small connected sets during mask connectivity enforcement
    return true;
  }
  
  if ( info[1] == "mskin-known" ) {
    envelope_consensus_mask_settings.input_known_mask = true;
    if ( 2 < info.size() ) envelope_consensus_mask_settings.input_filename_known_mask = info[2]; // Name of the input known mask file for the target structure. Useful for checking algorithm performance
    return true;
  }
  
  if ( info[1] == "fracsolvent" ) {
    if ( 2 < info.size() ) envelope_consensus_mask_settings.fraction_solvent_target = clipper::String(info[2]).f64(); // Fractional solvent content. Required input.
    return true;
  }
  
  if ( info[1] == "apodization-limit" ) {
    envelope_consensus_mask_settings.apodization = true;
    if ( 2 < info.size() ) envelope_consensus_mask_settings.apodization_limit = clipper::String(info[2]).f64(); // If specified then an apodization function will be used to calculate weights, which will be applied when correlating the Fourier amplitudes from the mask with the experimental amplitudes. This is the effective resolution limit = the resolution at which the apodization function is 1/20th of the maximumum.
    return true;
  }
  
  if ( info[1] == "mtzin-target" ) {
    envelope_consensus_mask_settings.input_Fourier_amplitudes = true;
    if ( 2 < info.size() ) envelope_consensus_mask_settings.input_filename_target_data = info[2]; // Name of the mtz file carrying the Fourier amplitudes for the target structure
    return true;
  }
  
  if ( info[1] == "distance-measure" ) {
    if ( 2 < info.size() ) envelope_consensus_mask_settings.distance_measure = clipper::String(info[2]).f64(); // integer indicating the distance measure to be used.
    return true;
  }
  
  if ( info[1] == "target-fo" ) {
    envelope_consensus_mask_settings.input_Fourier_amplitudes = true;
    if ( 2 < info.size() ) envelope_consensus_mask_settings.target_col_fo = info[2]; // Labels for the Fourier Amplitudes and their estimated standard deviations
    return true;
  }
  
  if ( info[1] == "shift-to-known-origin" ) {
    envelope_consensus_mask_settings.shift_to_known_origin = true; // If true, the origin shifts required to maximize agreement with the known mask will be applied to the concensus masks 
    return true;
  }

  if ( info[1] == "allow-inversion" ) {
    envelope_consensus_mask_settings.allow_inversion = true; //  If true, inversion symmetry will be used  when calculating distances for clustering (if the space group is achiral ... ). During averaging to produce a cluster consensus the input masks will be inverted as required
    return true;
  }

  if (info[1] == "verbose") // allow verbose setting for envelope consensus only
  {
    bool verbose_setting = true; // defaults to true,  , if there is no explicit true/false assignment 
    if ( 2 < info.size() ) {
      verbose_setting = convert_string_into_bool(info[2]);
    }
    envelope_consensus_mask_settings.verbose = verbose_setting;
    return true;
  }
  
  if ( info[1] == "list-file" ) {
    if ( 2 < info.size() ) envelope_consensus_mask_settings.mask_list_file = info[2]; // Custom list-file independant of job name.
    return true;
  }
  
  {
    outlog << "Unrecognized argument:\t" << info[1] << "\n";
    return false;
  }
}

bool exp_manager::read_phase_consensus_setting_into_struct(std::vector<std::string>& info, std::ostream& outlog)
{
  if ( info[1] == "col-fp" ) {
    if ( 2 < info.size() ) phase_consensus_phase_settings.col_fp = info[2]; // labels for the amplitudes and phases in the input mtz files
    return true;
  }
  if ( info[1] == "disable-origin-shifts" ) {
    phase_consensus_phase_settings.allow_origin_shifts =false;
    return true;
  }
  if ( info[1] == "epsilon" ) {
    phase_consensus_phase_settings.threshold_specified = true;
    if ( 2 < info.size() ) phase_consensus_phase_settings.epsilon_input = clipper::String(info[2]).f64(); // defines the threshold distance (mean absolute phase difference) used for clustering
    return true;
  }
  if ( info[1] == "min-points" ) {
    if ( 2 < info.size() ) phase_consensus_phase_settings.min_pts = clipper::String(info[2]).i(); // defines the minimum cluster size
    return true;
  }
  if ( info[1] == "effective-res-limit" ) {
    phase_consensus_phase_settings.apodization = true;
    if ( 2 < info.size() ) phase_consensus_phase_settings.effective_resolution_limit = clipper::String(info[2]).f64(); // If specified then an apodization function will be used to calculate weights for computation of the phase agreement statistics. This is the effective resolution limit = the resolution at which the apodization function is 1/20th of the maximumum.
    return true;
  }
  if ( info[1] == "known-solution-phase-filename" ) {
    phase_consensus_phase_settings.input_known_phases = true;
    if ( 2 < info.size() ) phase_consensus_phase_settings.input_filename_known_phase_set = info[2]; // Name of the known phase set for the target structure. Useful for checking algorithm performance
    return true;
  }

  if ( info[1] == "shift-to-known-origin" ) {
    phase_consensus_phase_settings.shift_to_known_origin = true; // If true, the origin shifts required to maximize agreement with the known phase set will be applied to the concensus phase set
    return true;
  }
  if ( info[1] == "known-solution-fp-column-labels" ) {
    if ( 2 < info.size() ) phase_consensus_phase_settings.known_phase_set_fp_column_labels = info[2]; // labels for the amplitudes and phases ... known phase set ...
    phase_consensus_phase_settings.input_known_phases = true;
    return true;
  }
  if ( info[1] == "list-file" ) {
    if ( 2 < info.size() ) phase_consensus_phase_settings.phase_list_file = info[2]; // Custom list-file independant of job name.
    return true;
  }
  
  if (info[1] == "verbose") // allow verbose setting for phase consensus only
  {
    bool verbose_setting = true; // defaults to true,  , if there is no explicit true/false assignment 
    if ( 2 < info.size() ) {
      verbose_setting = convert_string_into_bool(info[2]);
    }
    phase_consensus_phase_settings.verbose = verbose_setting;
    return true;
  }
  
  
  {
    std::cout << "Unrecognized argument:\t" << info[2] << "\n";
    return false;
  }
}

bool exp_manager::convert_string_into_bool(std::string input)
{
  /*
   This helper function checks common ways people may input a boolean as a string, to prevent user input hiccups.
   Also provides a standardised place for use to parse strings into bools.
   */
  bool output;
  lower_string(input); // Remove capitalization, to simplify checking
  
  // Systematically check for common ways to implement a bool:
  if (input == "true" || input == "t" || input == "1" || input == "tru" || input == "true.") output = true;
  else
  if (input == "false" || input == "f" || input == "0" || input == "fals" || input == "false.") output = false;
  else n_bad_inputs_count++; // increment bad inputs if there's something not caught.
  
  return output;
}

void exp_manager::initialise_hkl_info_and_data_objects_for_settings(ipa_settings& settings)
{
   // With this function we ensure that the correct hkl info and data objects are created within the ipa_settings structs for the envelope and/or phase determination steps.
   // Importantly, because it is encapsulated within the struct object, the whole object can be copied, without breaking the hkl_info and hkl_data connection during direct copying.
   // This step requires that the scale and B-factor have been estimated, and the mean intensities for the data set estimated, which happens when we execute Rogers Analysis.
  
  // Member variables we should have initialised:
  //clipper::HKL_info hkl_target; // Hkl info objects to handle the data below.
  //clipper::HKL_info hkl_wilson;
  //clipper::HKL_data <clipper::data64::F_sigF> measured_f_all; // All of the provided data from the user.
  //clipper::HKL_data <clipper::data64::F_sigF>  mean_i_wilson; // This needs to be copied to the ipa_manager and worker.
  
  // We start off by copying the respective hkl info objects to the hkl infos within the struct:
  // Initial testing confirmed we can simply copy these objects as such...
  settings.hkl_target = hkl_target;
  settings.hkl_wilson = hkl_wilson;
  
  // Now we need to initialise the settings.dataobjects with the appropriate hkl info.
  settings.measured_f_all.init(settings.hkl_target, settings.hkl_target.cell());
  settings.mean_i_wilson.init(settings.hkl_wilson, settings.hkl_wilson.cell());
  
  // Now copy in all the data for both. As they may be at different resolution (meaning the reflection lists will be different) they must be copied in separate loops.
  clipper::HKL_data_base::HKL_reference_index ih;
  for (ih = measured_f_all.first(); !ih.last(); ih.next() )
  {
    settings.measured_f_all[ih] = measured_f_all[ih];
  }
  
  for (ih = mean_i_wilson.first(); !ih.last(); ih.next() )
  {
    settings.mean_i_wilson[ih] = mean_i_wilson[ih];
  }
  
  // The settings struct should be updated, and ready to go.
}

inline b_factor_data exp_manager::estimate_B_from_diffraction_data(clipper::HKL_data<clipper::data64::F_sigF>& measured_f, std::ofstream& bfactorlog)
{
  // The function estimate_B_from_diffraction_data is inline as it uses Fortran subroutines (module atom_count, module rogers_analysis) to perform all the core functions, and hence can't be defined externally
  
  double estimate_B_from_diffraction_data;
  
  // Set up logging of the output.
  // Because we have to interleave output from the C++ calling function and the Fortran subroutines, it is a bit more complicated than usual
    
  char output_filename[128];
  snprintf(output_filename, sizeof(output_filename), "%s", summary_log_filename.c_str()); // The file where all output from this function will be logged
  
  //std::ofstream bfactorlog = summary_log; // We do this, so we don't need to rename things...

  //bfactorlog.open(output_filename, std::ofstream::app); // Associate our filename with the stream object.

  bfactorlog << "Estimation of scale, the overall anisotropic B-factor, and the equivalent isotropic B-factor by the method of Rogers" << std::endl;
  bfactorlog << "            As implemented by Robert Blessing and David Langs Acta Crystallogr A. 1988 Sep 1;44:729–35.             " << std::endl;
  bfactorlog << "--------------------------------------------------------------------------------------------------------------------\n" << std::endl;

  bfactorlog << "**Estimating atomic composition of the crystallographic asymmetric unit**\n" << std::endl;
          
  bfactorlog.close(); // temporarily dissociate the filename from the stream object so the Fortran routines can open and append to the file. Seems safest to close() the stream here before allowing this
   
     
  // Function for estimating B using Patterson Origin analysis (Method of Rogers/Blessing and Lang)
  
  // Step 1 - Atom counting
  // TODO  ... split off atom counting into a seperate routine for clarity
  //--------------------------------------------------------
  
  // Estimate the elemental composition of the asymmetric unit (excluding solvent) based on protein or nucleic acid sequence information
  // In the future could implement a method for imput of "additional" atoms - associated with bound ligands, post-translational modifications etc etc
    
  // Variables needed interior to the function for atom counting.
  // We will work with all elements to make things extensible - obviously only 6 elements are actually relevant if we are sticking with protein and nucleic acid
  
  //  Create an array of string objects carrying the element symbols, up to atomic number 98 (Californium !!) That should do it :-)
  
  std::string e_symbols[] =
  {"H ","HE","LI","BE","B ","C ","N ","O ","F ","NE",
    "NA","MG","AL","SI","P ","S ","CL","AR","K ","CA",
    "SC","TI","V ","CR","MN","FE","CO","NI","CU","ZN",
    "GA","GE","AS","SE","BR","KR","RB","SR","Y ","ZR",
    "NB","MO","TC","RU","RH","PD","AG","CD","IN","SN",
    "SB","TE","I ","XE","CS","BA","LA","CE","PR","ND",
    "PM","SM","EU","GD","TB","DY","HO","ER","TM","YB",
    "LU","HF","TA","W ","RE","OS","IR","PT","AU","HG",
    "TL","PB","BI","PO","AT","RN","FR","RA","AC","TH",
    "PA","U ","NP","PU","AM","CM","BK","CF"};
  
  const std::vector<std::string> element_symbols(e_symbols, e_symbols + sizeof(e_symbols) / sizeof(std::string)); // Initialize a vector carrying the element symbols
  
  std::vector<int> protein_atom_counts(n_elements,0); // here is the correspondent array for the protein atom counts, all values initially set to zero
  std::vector<int> dna_atom_counts(n_elements,0); // here is the correspondent array for the dna atom counts, all values initially set to zero
  std::vector<int> rna_atom_counts(n_elements,0); // here is the correspondent array for the rna atom counts, all values initially set to zero
  std::vector<int> total_atom_counts(n_elements,0); // here is the correspondent array for the total atom counts, all values initially set to zero
  
  
  if (protein_sequence_is_provided || dna_sequence_is_provided || rna_sequence_is_provided) // user has provided protein or nucleic acid sequence information for the target
  {
    
    if (protein_sequence_is_provided)
    {
      count_protein_atoms(protein_sequence_filename, &protein_sequence_nchains, &protein_atom_counts[0], output_filename); // pass everything to the fortran subroutine by reference (character variables are passed this way by default)
    }
    
    if (dna_sequence_is_provided)
    {
      bool is_dna = true;
      count_na_atoms(dna_sequence_filename, &is_dna, &dna_sequence_nchains, &dna_atom_counts[0], output_filename); // pass everything to the fortran subroutine by reference (character variables are passed this way by default)
    }
    
    if (rna_sequence_is_provided)
    {
      bool is_dna = false;
      count_na_atoms(rna_sequence_filename, &is_dna, &rna_sequence_nchains, &rna_atom_counts[0], output_filename); // pass everything to the fortran subroutine by reference (character variables are passed this way by default)
    }
    
    
    bfactorlog.open(output_filename, std::ofstream::app); // We are done with Fortran for a bit so reassociate our filename with the stream object. But now we *append* to the existing file contents.
    
    bfactorlog << "\nSumming Contributions from all Macromolecular Components" << std::endl;
    bfactorlog << "\nAtomic composition of the crystallographic asymmetric unit: \n" << std::endl;
    
    for ( int i= 0; i < n_elements; i++ )
    {
      total_atom_counts[i] = protein_atom_counts[i]*protein_sequence_ncopies_in_asu +
      dna_atom_counts[i]*dna_sequence_ncopies_in_asu     +
      rna_atom_counts[i]*rna_sequence_ncopies_in_asu;
      if (total_atom_counts[i] > 0)
      {
        bfactorlog << i+1 << " " << element_symbols[i] << " " << total_atom_counts[i] << std::endl;
      }
    }
    
  }
  else // no sequence information for the target is provided, estimate the elemental composition 
    
  {
    // TODO: ...  estimate the actual number of atoms (rather than the base atomic ratio) by making use of the solvent fraction and protein density.
    // Then we will always get a decent estimate of the scale factor k, which we can use to check the scaling of the input data 

    bfactorlog.open(output_filename, std::ofstream::app); // reassociate our filename with the stream object. But now we *append* to the existing file contents.
    
    bfactorlog << "No input data allowing for exact atom counts - approximating under the assumption that the crystal contains only protein" << std::endl;
    bfactorlog << "\nBase Atomic composition within the crystallographic asymmetric unit: \n" << std::endl;

    // 1, 1, 3 O, N, C is the approximate elemental ratio for a protein
    for ( int i = 0; i < n_elements; i++ )
    {
      if ( i == 5) total_atom_counts[i] = 3; // C
      if ( i == 6) total_atom_counts[i] = 1; // N
      if ( i == 7) total_atom_counts[i] = 1; // O
      if (total_atom_counts[i] > 0)
      {
        bfactorlog << i+1 << " " << element_symbols[i] << " " << total_atom_counts[i] << std::endl;
      }
    }
  }
  
  //------------------------------------------------------------------------
  // End of Step 1 atom counting - output is total_atom_counts - the atom counts in the asymmetric unit for all elements
  // total_atom_counts is indexed 0,n_elements-1, you get the atomic number by adding 1 to the array index ...
  
  // Step 2 : Estimate scale and overall B using the method of Rogers
  
  // This is a C++ replacement for the original main program written by Blessing and Lang, which calls the needed Fortran subroutines
  
  // Note carefully that in the original code Blessing and Lang define as s = sin(theta)/lambda, whereas conventionally s = the magnitude of the scattering vector = 2*sin(theta)/lambda = 1/resolution
  // Here we have reverted to the more common definition
  
  // Variables associated with control of the algorithm. Locally defined, but we could consider putting these under user control
  
  double rogers_cutoff = -3.0; // The I/sigma(I) cutoff used in the procedure. Default is 1.0
  double rogers_fraction = 0.5; // a fractional multiplier by which to multiply shift values between least squares cycles
  double rogers_smin =  0.0; //  if s = 2sin(theta)/lambda < smin, observation is omitted from the least-squares analysis
  double rogers_smax = 18.0; //   if s = 2sin(theta)/lambda > smax, observation is omitted from the least-squares analysis

  bfactorlog << "\n***Beginning Rogers analysis***\n" << std::endl;
  bfactorlog << "Parameters controlling the algorithm Rogers:" << std::endl;
  bfactorlog << "Threshold F/sigma(F) for including observations in the Patterson summation: " << rogers_cutoff << std::endl;
  bfactorlog << "Fractional shift multiplier for least-squares refinement: " << rogers_fraction << std::endl;
  bfactorlog << "Minimum 2sin(theta)/lambda for including observation in the Patterson summation (1/Å): " << rogers_smin << std::endl;
  bfactorlog << "Maximum 2sin(theta)/lambda for including observation in the Patterson summation (1/Å): " << rogers_smax << "\n" << std::endl;
  
  
  bfactorlog << "TARGET STRUCTURE -- CRYSTALLOGRAPHIC INFO:\n";
  
  bfactorlog << " Space group " << measured_f.hkl_info().spacegroup().symbol_xhm() << std::endl;
  char target_cf = util::crystal_family( measured_f.hkl_info().spacegroup().spacegroup_number() );
  bfactorlog << " Crystal Family " << target_cf << std::endl;
  bfactorlog << " Lattice Multiplicity " << util::lattice_multiplicity( measured_f.hkl_info().spacegroup().symbol_xhm()) << std::endl;
  bfactorlog << measured_f.hkl_info().cell().format() << std::endl;
  bfactorlog << " Number of Symmetry Operators " << measured_f.hkl_info().spacegroup().num_symops() << std::endl;
  
  float target_min_resn = 1.0/sqrt(measured_f.invresolsq_range().min());
  float target_max_resn = 1.0/sqrt(measured_f.invresolsq_range().max());
  
  bfactorlog << " Resolution range " <<  target_min_resn << " - " << target_max_resn << std::endl;
  bfactorlog << " Number of data - possible " << measured_f.hkl_info().num_reflections() << " - observed " << measured_f.num_obs() << std::endl;
  
  
  // Convert the atomic composition of the asymmetric unit into the atomic composition of the unit cell
  
  float num_symops = measured_f.hkl_info().spacegroup().num_symops();
  for ( int i= 0; i < n_elements; i++ ) { total_atom_counts[i] *= num_symops; }
  
  bfactorlog << "\nAtomic composition of the unit cell: \n" << std::endl;
  
  for ( int i= 0; i < n_elements; i++ )
  {
    if (total_atom_counts[i] > 0) { bfactorlog << i+1 << " " << element_symbols[i] << " " << total_atom_counts[i] << std::endl;}
  }
  
  // Get the Atomic scattering factors (Cromer-Mann coefficients) and relative atomic masses for each element
  
  struct Cromer_Mann_coefficients
  {
    double a1;
    double b1;
    double a2;
    double b2;
    double a3;
    double b3;
    double a4;
    double b4;
    double c0;
  };
  
  struct anomalous_scattering_components
  {
    double fp;
    double fpp;
  };
  
  std::vector<Cromer_Mann_coefficients> cmc(n_elements);
  std::vector<anomalous_scattering_components> anom(n_elements);
  std::vector<double> rm(n_elements);
  
  
  double F000=0.0; // The amplitide of the zeroeth order Fourier term in electrons
  double RMCELL=0.0; // The summed relative atomic masses of all atoms in the unit cell
  double SUMZSQ = 0.0; // Required for the fitting of the Patterson function in subroutine PFIT.
  
  for (int i = 0; i < n_elements; i++)
  {
    int iz = i+1; // the atomic number of the element being referenced at index i in the C vector
    FTABLE(&iz, &rm[i], &cmc[i].a1, &cmc[i].b1, &cmc[i].a2, &cmc[i].b2, &cmc[i].a3, &cmc[i].b3, &cmc[i].a4, &cmc[i].b4, &cmc[i].c0); // Retrieves the relative atomic mass and Cromer-Mann Coefficients describing the normal X-ray scatteting of each element
    
    // TODO ... insert a subroutine here to return the anomalous scattering corrections f' and f''. For now we simply set these equal to zero
    // The original version of FTABLE would supply the corrections relevant to either Cu or Mo Kalpha radiation.
    // Interestingly these corrections were not used in the calculation of F000 and SUMZSQ, though technically, I can't see why not
    
    anom[i].fp = 0.0; // anomalous scattering component f'
    anom[i].fpp = 0.0; // anomalous scattering component f''
    
    F000   += (double)total_atom_counts[i]*(double)iz; // scattering contribution at zero angle is equal to the atomic number for a neutral atom
    RMCELL += (double)total_atom_counts[i]*rm[i];
    SUMZSQ += std::pow( (double)total_atom_counts[i]*(double)iz, 2);
  }
  
  
  // We use C++ vectors to reserve contiguous blocks of memory for multi-dimensional arrays used by the Fortran subroutines and functions.
  // So long as they are the same size in both C++ and Fortran, and we never operate on these vectors on the C++ side, nothing can go wrong :-)
  
  std::vector<double> CELL(6); // The unit cell dimensions
  std::vector<double> G(3*3); // The metric matrix (metric tensor) G - transforms (a*, b*, c*) to (a, b, c)
  std::vector<double> GINV(3*3); // The inverse of the metric matrix (metric tensor) - transforms (a, b, c) to (a*, b*, c*)
  double VCELL; // The Volume of the Unit Cell
  
  // Note that cell angles are stored as radians within hkl_info, while the fortran subroutines expect values in degrees
  
  CELL[0] =  measured_f.hkl_info().cell().a();
  CELL[1] =  measured_f.hkl_info().cell().b();
  CELL[2] =  measured_f.hkl_info().cell().c();
  CELL[3] =  (180.0/clipper::Util::pi())*measured_f.hkl_info().cell().alpha();
  CELL[4] =  (180.0/clipper::Util::pi())*measured_f.hkl_info().cell().beta();
  CELL[5] =  (180.0/clipper::Util::pi())*measured_f.hkl_info().cell().gamma();
  
  // outlog << CELL[0] <<  " " << CELL[1] <<  " " << CELL[2] <<  " " << CELL[3] <<  " " << CELL[4] <<  " " << CELL[5]  << std::endl;
  
  METRIC(&CELL[0], &VCELL, &G[0], &GINV[0]); // Calculate VCELL, G and GINV from the unit cell parameters.
  
  bfactorlog << "\nCell Volume (Å^3): " << VCELL << std::endl;
  
  // Compute density of the crystal by dividing relative atomic mass in the unit cell by the volume of the unit cell.
  // 0.60221 is a conversion factor to get the correct units (g/cm^3 = mg/mm^3)
  // See e.g. X-Ray Crystallography and Crystal Packing Analysis in "Solid State Properties of Pharmaceutical materials"
  
  double DX = (RMCELL/VCELL)/0.60221;
  
  bfactorlog << "F(000) (electrons): " << F000 << std::endl;
  bfactorlog << "Summed Relative Atomic Masses of the Atoms in the Unit Cell: " << RMCELL << std::endl;
  bfactorlog << "Crystal Density, considering the counted atoms (g/cm^3): " << DX <<  std::endl;
  
  // Move functionality of former subroutines DATA and READ1 into the main program. Basic manipulations of the input data are more robustly handled via the clipper library.
  
  clipper::HKL_info hkl =  measured_f.base_hkl_info();
  clipper::HKL_data<clipper::data32::Flag>  rejection_flag( hkl ); // create a data object of type Flag to handle data selection
  clipper::HKL_data<clipper::data64::I_sigI>    measured_I( hkl ); // create a data object of type I_SIGI to store the intensities I
  clipper::HKL_data<clipper::data64::I_sigI>    transformed_I( hkl ); // create a data object of type I_SIGI to store the transformed intensities I/ε
  
  clipper::HKL_data_base::HKL_reference_index ih;
  
  // The aim here is to treat the data exactly as done by Blessing and Lang in subroutines DATA and READ1
  // What we do here could subsequently adjusted, but for testing purposes it needs to work equivalently
  
  bfactorlog << "\nCreating transformed intensity data (I/ε) from the input Fourier amplitudes |F|, and expanding to a full hemisphere (space group P1)." << std::endl;
  bfactorlog << "Number of observations in the input mtz: " << measured_f.num_obs() << std::endl;

  // Note that we compute and store the untransformed intensities, I, to enable subsequent comparison with the Wilson model
  
  int n_data_valid = 0;
  
  for ( ih = hkl.first(); !ih.last(); ih.next() )
  {
    if ( !measured_f[ih].missing() )
    {
      
      // Invalid estimate for the standard deviation results in rejection, following Blessing and Lang
      if (measured_f[ih].sigf() < 0)
      {
        rejection_flag[ih].flag() = 1;
        continue; // loop immediately to the next observation
      }
      
      n_data_valid += 1;
      
      // Generate estimates for I and Sig(I), following Blessing and Lang
      
      double    f   = measured_f[ih].f();
      double sigf   = measured_f[ih].sigf();
      double weight = 1.0/ih.hkl_class().epsilon();
      
      // outlog << ih.hkl().h() << " " << ih.hkl().k() << " " << ih.hkl().l() << " " << weight << std::endl;
      
      if (f < 0)
      {
           measured_I[ih].I() = 0;
        transformed_I[ih].I() = 0;
      }
      else
      {
           measured_I[ih].I() = std::pow( f, 2 );
        transformed_I[ih].I() = weight*std::pow( f, 2 ); // transform intensity using statistical correction factors epsilson, prior to expansion to P1
        
      }
      
         measured_I[ih].sigI() = std::max( 2*f*sigf, 4*std::pow( sigf, 2 ) );
      transformed_I[ih].sigI() = weight*std::max( 2*f*sigf, 4*std::pow( sigf, 2 ) ); // transform standard deviation using statistical correction factors epsilson, prior to expansion to P1
      
      // Now work through some rejection criteria, following Blessing and Lang
      
      if (transformed_I[ih].I() < rogers_cutoff*transformed_I[ih].sigI()) // eliminated by the I/Sig(I) cutoff
      {
        rejection_flag[ih].flag() = 1;
        continue; // loop immediately to the next observation
      }
      if ( (sqrt(ih.invresolsq()) < rogers_smin) || ( sqrt(ih.invresolsq()) > rogers_smax) ) // eliminated by the resolution cutoffs
      {
        rejection_flag[ih].flag() = 1;
        continue; // loop immediately to the next observation
      }
      
    }
  }
  
  bfactorlog << "Number of valid observations in the input mtz (): " << n_data_valid << std::endl;
  
  transformed_I.mask( rejection_flag != 1 );  // now apply the rejection criteria to the transformed dataset
  
  bfactorlog << "Number observations following application  of the rejection criteria: " << transformed_I.num_obs()  << std::endl;
  
  // Finally expand the transformed Intensity data (In = I/ε) to the full hemisphere (space group p1 ) and format for access by the remaining Fortran procedures
  // The expansion is done slightly differently in the original code (subroutine data), where I'm not sure the treatment of statistical weights ε is strictly correct
  // Hence the numerical results of this re-implementation differ very very slightly from the original code.
  
  clipper::HKL_info hkl_p1( clipper::Spacegroup( clipper::Spgr_descr( 1 ) ), measured_f.base_hkl_info().cell(), measured_f.base_hkl_info().resolution() ); // New reflection list. Space group p1 with original cell and resolutioon
  hkl_p1.generate_hkl_list(); // initially there are no observations hkl in the list - better make some
  clipper::HKL_data<clipper::data64::I_sigI> transformed_I_p1(hkl_p1);
  
  
  //clipper::HKL_data_base::HKL_reference_index ih;
  for ( ih = hkl_p1.first(); !ih.last(); ih.next() ) {
    transformed_I_p1[ih] = transformed_I[ih.hkl()];
  }
  
  bfactorlog << "Number of observations following expansion to space group P1: " << transformed_I_p1.num_obs()  << std::endl;
  
  int n_intensity_data = transformed_I_p1.num_obs();
  
  std::vector<int> index_h(n_intensity_data);
  std::vector<int> index_k(n_intensity_data);
  std::vector<int> index_l(n_intensity_data);
  std::vector<double> normalized_intensity(n_intensity_data);
  
  int counter = 0;
  for ( ih = hkl_p1.first(); !ih.last(); ih.next() )
  {
    if ( !transformed_I_p1[ih].missing() )
    {
      // Compute the sum of the squared scattering factors for all atoms in the unit cell,  and use it to normalize the experimental intensities
      // The data now covers a full hemisphere, and the intensities have previously been divided through by the statistical correction factors, ε.
      // This subsumes the functionality of the original subroutine FCALC - it's more convenient to handle that here.
      
      double s_squared =  ih.invresolsq(); // (amplitude of the scattering vector)^2
                                           //double s = sqrt(s_squared); // amplitude of the scattering vector
      double sinthl_squared = s_squared/4; // // (sin(theta)/lambda)^2
                                           //double sinthl = s/2; // sin(theta)/lambda
      
      double SUMFSQ=0;
      
      for (int i = 0; i < n_elements; i++)
      {
        SUMFSQ +=  (double)total_atom_counts[i]
        * (  std::pow( cmc[i].a1*exp(-cmc[i].b1*sinthl_squared ) +
                      cmc[i].a2*exp(-cmc[i].b2*sinthl_squared ) +
                      cmc[i].a3*exp(-cmc[i].b3*sinthl_squared ) +
                      cmc[i].a4*exp(-cmc[i].b4*sinthl_squared ) +
                      cmc[i].c0 + anom[i].fp , 2)
           + std::pow( anom[i].fpp, 2) ) ;
      }
      
      //FCALC(&n_elements, &total_atom_counts[0], &cmc[0].a1, &cmc[0].b1, &cmc[0].a2, &cmc[0].b2, &cmc[0].a3, &cmc[0].b3, &cmc[0].a4, &cmc[0].b4, &cmc[0].c0, &anom[0].fp, &anom[0].fpp, &S, &SUMFSQ);
      
      index_h[counter]   = ih.hkl().h();
      index_k[counter]   = ih.hkl().k();
      index_l[counter]   = ih.hkl().l();
      normalized_intensity[counter]  = transformed_I_p1[ih].I()/SUMFSQ;
      
      counter += 1;
      
    }
  }
  
  //-------------------------------------------------------------------------------------------------------------------------------
  // Constants controlling dimensioning of arrays associated with the Patterson function
  // Any changes made to n_ax, NR and NA *must* be mirrored in the shared declaration section of the Fortran module rogers_analysis
  
  const int n_ax=200;  // number of grid points for the calculation of the Patterson along the axial directions a, b, c
  const int NR=25;     // radial grid spacing
  const int NA=12;     // angular grid spacing, this must be a multiple of 4.
  const int n_pc = 1 + NR + (NR*NA*NA)/4; // total number of points in the polar coordinate grid
  //-------------------------------------------------------------------------------------------------------------------------------
  
  std::vector<double> Patterson_axial(3*n_ax); // Stores the Patterson function evaluated along the axial directions
  double PMIN; // The minimum Patterson density along the axial directions
  double P0; // An estimate for the Amplitude of the trivariate Gaussian function - see eqn 6 of Blessing and Lang
  std::vector<double> PP(3*3); // An estimate for the Matrix specifying the trivariate Gaussian function - see eqn 6 of Blessing and Lang
  double RADIUS; // The half-height radius (Angstroms) of the Patterson origin peak.
  
  bfactorlog.close(); // temporarily dissociate the filename from the stream object so the Fortran routines can open and write to the file. Seems safest to close() the stream here before allowing this
  
  
  PRELIM(&n_intensity_data, &index_h[0], &index_k[0], &index_l[0], &normalized_intensity[0], &CELL[0], &VCELL, &GINV[0], &Patterson_axial[0], &PMIN, &P0, &PP[0], &RADIUS);   // Evaluate the Patterson function along the axial directions
  
  
  PLOT1(&Patterson_axial[0], output_filename); // Plot the Patterson along the axial directions
  
  std::vector<double> U(n_pc); // Coordinates of the grid for evaluation of the Patterson function - U axis
  std::vector<double> V(n_pc); // Coordinates of the grid for evaluation of the Patterson function - V axis
  std::vector<double> W(n_pc); // Coordinates of the grid for evaluation of the Patterson function - W axis
  std::vector<double> Patterson(n_pc); // The Patterson function evaluated at the grid points specified in U,V,W
  
  GRID(&RADIUS, &GINV[0], &U[0] , &V[0], &W[0]); // Calculate the grid for evaluation of the Patterson function
  PSUM(&n_intensity_data, &index_h[0], &index_k[0], &index_l[0], &normalized_intensity[0], &VCELL, &U[0], &V[0] ,&W[0], &Patterson[0]); // Evaluate the Patterson function at the grid points
  
  std::vector<double> Patterson_sculpted(n_pc); //
  
  // Make a copy of the Patterson that can be modified
  
  for (int i = 0; i < n_pc; i++)
  {
    Patterson_sculpted[i] = Patterson[i];
  }
  
  PADJ(&Patterson_sculpted[0], &PMIN, &P0); // Adjust the Origin peak to compensate for the missing F(000) term in the Patterson synthesis
  
  int NCYCLE; // The number of least squares cycles employed
  double R; // Measure of fit between the caclulated and observed Patterson origin peak
  
  // Now Fit the (modified) Patterson Origin peak with a Gaussian
  PFIT(&U[0] , &V[0], &W[0], &Patterson_sculpted[0], &PMIN, &P0, &PP[0] , &VCELL, &F000, &SUMZSQ, &rogers_fraction, &NCYCLE, &R);
  
  // Plot the Results
  PLOT2(&U[0] , &V[0], &W[0], &Patterson[0], &PMIN, &P0, &PP[0], &RADIUS, &G[0], output_filename);
  
  double SCALEK; // Absolute Scale Factor
  std::vector<double> U_CIF(3*3);  // Symmetric tensor U_CIF describing atomic displacement - for definition see see Grosse-Kunstleve & Adams J Appl Cryst. 2002;35(4):477–80.
  std::vector<double> U_STAR(3*3); // Symmetric tensor U_STAR describing atomic displacement - for definition see see Grosse-Kunstleve & Adams J Appl Cryst. 2002;35(4):477–80.
  std::vector<double> B(3*3); //  Symmetric tensor B describing atomic displacement B = 2π^2 * U_STAR
  double BISO; // Equivalent isotropic displacement parameter B
  double UISO; // Equivalent isotropic displacement parameter U
  
  // Recover the ADPs etc from the fitted Gaussian
  SOLVE(&P0, &PP[0], &VCELL, &G[0], &GINV[0], &target_cf, &SCALEK, &B[0], &U_CIF[0], &U_STAR[0], &BISO, &UISO, output_filename);
  
  // Print the Results
  PRINT_SOLUTION(&RADIUS, &NCYCLE, &R, &SCALEK,  &B[0], &U_CIF[0], &BISO, &UISO, output_filename);
  
  bfactorlog.open(output_filename, std::ofstream::app); // All done with the Fortran so reassociate our filename with the stream object. We *append* to the existing file contents.
    
  // Resolution limit for model based estimation of <I> may need to be greater than the experimental limit, to enable the free lunch algorithm
  // How to set the resolution extension sensibly and automatically?
  // 
  // If s_max_expt is the maximum |s| observed experimentally, we could augment by a constant volume in Fourier space.  
  // So the extended resolution limit s_max_extd would be  defined by  4/3 π s_max_extd ^3 - 4/3 π s_max_expt^3 = a constant (c) in Å^-3
  // hence s_max_extd = ∛ ((3/4π)c + s_max_expt^3)
  // 
  // If c were 0.05 Å^-3 this would give:
  // 5 Å extended to 3.69 Å
  // 3 Å extended to 2.73 Å
  // 2 Å extended to 1.94 Å
  //
  // This seems intuitively reasonable
  // 
  // As an alternative we could determine the extension limit by looking at where the isotropic average of the fitted <I/ε> fell below some threshold value.
  // For datasets which had been integrated out to the noise limit this would probably result in very little resolution extension
  // For datasets which had been "truncated" because of experimental issues (finite detector size) this might result in substantial resolution extension
  // This might be the desired behavior, as it is more closely tied to the actual properties of the data than a blanket resolution extension. 
  //
  

  double resolution;
  if (free_lunch)
  {
    resolution = util::calculate_new_resolution_limit_from_volume_scale_factor(measured_f.hkl_info().resolution().limit(), phasing_settings.free_lunch_scale_factor);
  }
  else
  {
    //resolution = 1.0/sqrt(measured_f.invresolsq_range().max());
    resolution = measured_f.hkl_info().resolution().limit(); // Pretty sure this is a neater way of retrieving this information...
  }
    
  hkl_wilson.init(measured_f.base_hkl_info().spacegroup(), measured_f.base_hkl_info().cell(), clipper::Resolution(resolution));
  hkl_wilson.generate_hkl_list(); // Make the hkl reflections.
  mean_i_wilson.init(hkl_wilson, hkl_wilson.cell()); // Initialise an I-SigI object, to store the expected values for the intensities, based on the Wilson model  
  
  double ml = util::lattice_multiplicity( hkl_wilson.spacegroup().symbol_xhm()); // the lattice multiplicity - will be non-unitary if cell is centered. 
  
  // bfactorlog << "lattice multiplicity: " << ml << std::endl;
  
  for ( ih = hkl_wilson.first(); !ih.last(); ih.next() )
  {
    // Compute the sum of the squared atomic scattering factors for all atoms in the unit cell using the stored Cromer-Mann coefficients and atom counts
    
    double s_squared =  ih.invresolsq(); // (amplitude of the scattering vector)^2
    double sinthl_squared = s_squared/4; // (sin(theta)/lambda)^2
    
    double sumfsq = 0;
    
    for (int i = 0; i < n_elements; i++)
    {
      sumfsq +=  (double)total_atom_counts[i]
      * (  std::pow( cmc[i].a1*exp(-cmc[i].b1*sinthl_squared ) +
                    cmc[i].a2*exp(-cmc[i].b2*sinthl_squared ) +
                    cmc[i].a3*exp(-cmc[i].b3*sinthl_squared ) +
                    cmc[i].a4*exp(-cmc[i].b4*sinthl_squared ) +
                    cmc[i].c0 + anom[i].fp , 2)
         + std::pow( anom[i].fpp, 2) ) ;
    }
    
    // Evaluate the Debye-Waller Factor from tensor U* (Could also use B since B = 2π^2 U*)
    double wa = util::calculate_debye_waller_factor(ih.hkl().h(), ih.hkl().k(), ih.hkl().l(), U_STAR); // Pretty sure this'll work.
        
    // retrieve the standard statistical weights epsilon for the current observation hkl
    double epsilon = ih.hkl_class().epsilon(); 

    double wa_squared = wa*wa;
    double one_over_ksquared = 1.0/(SCALEK*SCALEK);

    // mean_i_wilson[ih].I() = one_over_ksquared * wa_squared * ml * epsilon * sumfsq; // This is the expectation value for I
    // mean_i_wilson[ih].I() = one_over_ksquared * wa_squared * ml * sumfsq; // This is the expectation value for I/ε  
    
    // Note that if the lattice multiplicity ml were accounted for in Rogers estimation of k and B, its effects would be subsumed into the linear scale factor k
    // As it currently is not, we must include it as a multiplier here ... 
     
    mean_i_wilson[ih].I() = one_over_ksquared * wa_squared * ml * epsilon * sumfsq;
  }
  
  // Write out the Wilson model intensities, the experimental intensities, and the model residuals to a file, for diagnostic purposes 
  
  // Get the data ready for output
  
  clipper::HKL_data<clipper::data64::I_sigI>    model_residuals( hkl ); // create a data object of type I_SIGI to store the signed residuals. 
  clipper::HKL_data<clipper::data64::I_sigI>    model( hkl ); // create a data object of type I_SIGI to store the Wilson model at the experimental resolution
  
  for ( ih = hkl.first(); !ih.last(); ih.next() )
  {
   // Transfer Wilson model intensities into an object having the same resolution as the experimental data
    model[ih].I()    = mean_i_wilson[ih.hkl()].I(); 
    model[ih].sigI() = 0; // Standard deviations are irrelevant in this setting 

   // Calculate the signed model residuals
    if (!measured_I[ih].missing())
    {
      //model_residuals[ih].I()    = abs(measured_I[ih].I() - model[ih].I()); // Unsigned residuals ... less useful
      model_residuals[ih].I()    = measured_I[ih].I() - model[ih].I(); 
      model_residuals[ih].sigI() = 0; // Standard deviations are irrelevant in this setting 
    }
    
  }
  
  // Now output the data
  
  clipper::CCP4MTZfile mtzout;
  char output_mtz_filename[128];
  snprintf(output_mtz_filename, sizeof(output_mtz_filename), "%s", "Wilson-model.mtz"); 
  mtzout.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
  mtzout.open_write( output_mtz_filename);
  mtzout.export_hkl_info( hkl );
  
  clipper::String output_col_labels;
  
  output_col_labels = "/*/*/[I_measured,SIGI_measured]";
  mtzout.export_hkl_data(measured_I, output_col_labels ); 
  
  output_col_labels = "/*/*/[I_model,SIGI_model]";
  mtzout.export_hkl_data(model, output_col_labels ); 
  
  output_col_labels = "/*/*/[I_resid,SIGI_resid]";
  mtzout.export_hkl_data(model_residuals, output_col_labels ); 
  
  mtzout.close_write();
  
  //bfactorlog.close(); // close the stream object for good.

  mean_i_wilson_is_valid = true; // Flag to indicate that we have calculated < I >
  estimate_B_from_diffraction_data = BISO;

  return b_factor_data(UISO,BISO,SCALEK,U_STAR,U_CIF,B);  
}

