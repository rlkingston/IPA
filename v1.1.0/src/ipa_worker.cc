// IPA: Iterative Projection Algorithms for protein crystallography
// IPA worker class

#include "ipa_worker.h"
#include "ipa-lib.h"
#include "util-lib.h"

ipa_worker::ipa_worker(ipa_settings new_settings)
{
  /*
   Constructor for ipa_worker.
   */
  settings = new_settings;
}

bool ipa_worker::run_experiment(std::atomic<float> &ThreadProgress,
                                std::atomic<int> &IterateCounter,
                                std::atomic<float> &CUSUM,
                                const int run_id)
{
  // ------------------------------------------------------- //
  // ------------------ Initializing Run ------------------- //
  // ------------------------------------------------------- //
  my_run_id = run_id;

  ThreadProgress = 2.0; // Update Thread progress to Initializing.

  if (!initialise_run())
  {
    return false;
  };

  // ------------------------------------------------------- //
  // ------------------- Algorithm Loop -------------------- //
  // ------------------------------------------------------- //

  int total_rules = settings.UpdateRuleRegime.size();
  for (int rule = 0; rule < total_rules; rule++)
  {
    // For every UpdateRuleRegime, we perform requested number of iterates from each update rule container.
    ipa_setting_algorithm_container UpdateRegime = settings.UpdateRuleRegime[rule];
    int iterates_at_this_step = UpdateRegime.n_iterates;

    // Calculate LUT to quickly know which update rules and constraints to use.
    std::vector<IPA_updaterule_command> UpdateRule_LUT;
    std::vector<IPA_constraint_command> Constraint_LUT;

    // TODO Take the loops below and refactor, passing in the LUT probably.
    // These while loops assume that the UpdateRule is valid. If the n_commands or the n_constraints are 0, you'll get an infinite loop. Yay.
    while (UpdateRule_LUT.size() < iterates_at_this_step)
    {
      for (int i = 0; i < UpdateRegime.commands.size(); i++)
      {
        for (int j = 0; j < UpdateRegime.n_commands[i]; j++)
        {
          UpdateRule_LUT.push_back(UpdateRegime.commands[i]);
        }
      }
    }
    while (Constraint_LUT.size() < iterates_at_this_step)
    {
      for (int i = 0; i < UpdateRegime.constraints.size(); i++)
      {
        for (int j = 0; j < UpdateRegime.n_constraints[i]; j++)
        {
          Constraint_LUT.push_back(UpdateRegime.constraints[i]);
        }
      }
    }

    // Both LUTs are now fully populated, let's carry out the program ...
    for (int i = 0; i < iterates_at_this_step; i++)
    {
      run_loop(ThreadProgress, UpdateRule_LUT[i], Constraint_LUT[i]);
      total_iterate_counter++;  // Increment the total iterate counter.
      IterateCounter++;         // Tell the Thread Manager we performed an iterate too!
      CUSUM = stat_phase_cusum; // Update our CUSUM predictor of success (Experimental ...)
    }
  }

  // ------------------------------------------------------- //
  // --------------------- Closing Run --------------------- //
  // ------------------------------------------------------- //
  total_iterate_counter -= 1; // Decrement by one, as we're going to add a final iterate line regardless to the final run.
                              // We could set a bool if we already wrote to the summary file on the iterate, and thus avoid this duplicating behaviour... but a summary line might be nicer?
  if (!finish_run())
  {
    return false;
  }

  return true; // Successful run
}

bool ipa_worker::initialise_run()
{
  /*
   Initialises all settings and memory for a looping run. This includes importing and parsing all mtz data i.e. raw data, reference data, working and test masks etc.
   Calculates important looping Look Up Table values for settings such as filter_radius and beta.
   Further will check settings to ensure they are safe to run with.
   */

  // Create custom log file for this threaded run.
  char log_filename[64];
  snprintf(log_filename, sizeof(log_filename), "%slog_run%d.txt", settings.working_dir.c_str(), my_run_id);
  outlog.open(log_filename);
  if (!outlog.good())
  {
    std::cerr << "Bad log file generation!";
    return false; // Unsuccesful run
  }

  util::version_fingerprint(outlog);

  outlog << std::setprecision(8) << std::fixed;

  outlog << "\nBeginning Run Initialization" << std::endl;
  outlog << "----------------------------\n"
         << std::endl;

  total_iterate_counter = 0; // Set this nice and early.

  // PRINT ALL SETTINGS OUT FOR DEBUGGING.
  if (settings.verbose)
  {
    outlog << "\nA dump of all settings\n"
           << std::endl;
    util::print_all_settings_to_command_line(outlog, settings);
  }

  outlog << "\nTotal Number of iterates = " << settings.n_total_iterates << std::endl;

  // Initialise Random Number Generator: (Pseudo-random is done by default)
  double range_from = 0;
  double range_to = 1;
  uniform_distribution = std::uniform_real_distribution<double>(range_from, range_to);
  std::random_device rd;
  std::default_random_engine generator(rd());         // Random numbers for deployment.
  static std::default_random_engine generator_static; // Pseudo Random-numbers for DEBUG

  if (!settings.really_random) // Store whichever generator we want to use.
  {
    MyGenerator = generator_static;

    outlog << "Using a static RN generator." << std::endl;
    for (int m = 0; m < my_run_id; m++)
    {
      // Generate my_run_id's number of random numbers, ensuring static random number generation is not completely the same across threads.
      float randomint = uniform_distribution(MyGenerator);
    }
  }
  else
  {
    MyGenerator = generator;
    outlog << "Using a non-static RN generator." << std::endl;
  }

  //Report Magic number 
  outlog << "Magic Number: " << uniform_distribution(MyGenerator) << std::endl; // Runs with the same magic number when using the same random number generation.

  // The parameter LUT are passed directly through the settings object now.
  // Initialising the LUT for parameter beta.
  if (settings.found_beta_params)
  {
    calculate_ipa_beta_for_iterate(settings.n_total_iterates, true);
  }
  else
  {
    outlog << "Beta Iterate Inputs are bad." << std::endl;
    return false;
  }

  // Initialising the LUT for parameter filter_radius.
  if (settings.found_filter_params)
  {
    calculate_filter_radius_for_iterate(settings.n_total_iterates, true);
  }

  // Initialising the LUT for parameter solvent_fraction multiplier.
  calculate_solvent_fraction_multiplier_for_iterate(settings.n_total_iterates,
                                                    settings.found_solvent_fraction_multiplier_params,
                                                    true);

  // Populate the reference_structure[] with appropriate data.
  populate_reference_data();

  // read input mtz for target, populates HKL data pointers, populates resolution numbers for target. Also read in known phases if we have them (test cases etc)
  if (!initialise_data_objects())
  {
    outlog << "Bad Input mtz" << std::endl; // ... TODO: "Bad Input mtz" could result from a failure to read the target data or the known phase set. Split these operation apart for clarity ?
    return false;
  } // Bad input, abort.

  outlog << "\nIntensity statistics for the Target data\n"
         << std::endl;
  
  // Calculate and store <I/ε>, and number of observations, for centric and acentric data.
  // Deprecated - functionality subsumed by print_intensity_statistics_data_vs_model
  /*
  mean_I_acentric_measured_data.resize(settings.res_bins + 1, 0.0);
  mean_I_centric_measured_data.resize(settings.res_bins + 1, 0.0);
  observation_counts_acentric_measured_data.resize(settings.res_bins + 1, 0.0);
  observation_counts_centric_measured_data.resize(settings.res_bins + 1, 0.0);

  ipa_functions::compute_intensity_statistics(measured_f,
                                              settings.res_bins,
                                              mean_I_acentric_measured_data,
                                              mean_I_centric_measured_data,
                                              observation_counts_acentric_measured_data,
                                              observation_counts_centric_measured_data,
                                              outlog, true); // force verbose behavior here
  */
  
  if (settings.verbose)
  {
    outlog << "HIGH BINNED DATA AND MODEL OUTPUT: " << std::endl;
    util::print_intensity_statistics_data_vs_model(measured_f, mean_i_wilson, 500, outlog); // Extreme verbose output of data and model.
  }

  // Starting our working f as a copy of measured_f
  copy_F_sigF_data(measured_f, working_f);

  // Generate integer flags that define the work and test sets, store in test_set
  if (settings.use_test_set)
  {
    current_test_set = -1; // Will increment to 0 on first iterate.
    clipper::HKL_data_base::HKL_reference_index ih;
    for (ih = test_set_flags.first(); !ih.last(); ih.next())
    {
      test_set_flags[ih].flag() = lround((double(settings.n_test_sets - 1) * uniform_distribution(MyGenerator)));
    }
  }

  // Ensure all the test objects are filled when not using test sets.
  clipper::HKL_data_base::HKL_reference_index ih;
  for (ih = hkl_target_eff_res.first(); !ih.last(); ih.next())
  {
    work_set_f[ih] = measured_f[ih];
    test_set_f[ih] = measured_f[ih];
  }

  // Import known mask if provided (test cases only)
  // We don't interpolate onto the correct grid here - this is done within the run loop
  if (settings.input_known_mask)
  {
    outlog << "\nAn known Mask Is Input For The Target Structure: \n\n";
    if (!safely_import_mask(input_target_mask_known,
                            settings.input_filename_known_mask,
                            hkl_target_exp_res))
    {
      return false;
    } // Bad known mask input, stop run.
  }

  // Import a working mask if provided.
  // We don't interpolate onto the correct grid here - this is done within the run loop
  if (settings.input_working_mask)
  {
    outlog << "\nA Working Mask Is Input For The Target Structure: \n\n";
    if (!safely_import_mask(input_target_mask_working,
                            settings.input_filename_target_mask_working,
                            hkl_target_exp_res)) // hkl target is fine here as only space group is used.
    {
      return false;
    } // Bad working mask input, stop run.
  }

  // This was never done, so do it now...
  define_reference_density_histograms();
  // Set up the density value histogram for the target structure
  density_histogram.number_of_intervals = settings.n_density_bins;
  density_histogram.interval = new float[density_histogram.number_of_intervals];
  // Set up the gradient magnitude histogram for the target structure
  gradient_histogram.number_of_intervals = settings.n_gradient_bins;
  gradient_histogram.interval = new float[gradient_histogram.number_of_intervals];

  if (!select_reference_structure()) // Select and import a reference structure
  {
    return false;
  } // bad reference structure

  // Set the maximum resolution for this run:
  maximum_resolution = hkl_target_max_res.resolution().limit(); // The maximum resolution is carried in hkl object hkl_target_max_res.
  if (settings.free_lunch)
  {
    outlog << "Maximum resolution has been extended to " << maximum_resolution << " Å." << std::endl;
  }
  else
  {
    outlog << "Maximum resolution is " << maximum_resolution << " Å." << std::endl;
  }

  // Operationally, settings.apodize should always be true. So we might remove this when things calm down.
  if (settings.apodize) // we are apodizing the data - apodization may be fixed, or it may vary with iteration
  {
    outlog << "\nCalculating the apodization scheme" << std::endl;
    calculate_apodization_look_up_table(); // new way - could call this from below...
    initialise_apodization_scheme();
  }

  prev_rlimit = -1; // This will be set below.
  // ap_step_iterate_counter = 0; // Might be able to be moved now, as it occurs within prepare_data_objects.check_apodization
  outlog << "Preparing data objects (including Apodization) for the first iterate." << std::endl;
  check_if_histogram_or_gradient_matching();
  ipa_worker::prepare_data_objects(0); // Prepare all objects for the first iterate, the 0th apodization step.

  if (!settings.input_starting_phases)
  { // Generate Random Phases if no input phases were found.
    if (!generate_a_starting_phase_set())
    {
      return false;
    }
  }

  // Experimental platform for manipulating the centric phases, prior to a run ...
  // ############################################
  // #        TESTING CENTRIC PHASE             #
  // ############################################
  // Update with the centric starting phases if this is a thing.
  if (settings.keep_known_centric_phases)
  {
    outlog << "Copying the centric phases from the known phase set to the starting phase set. " << std::endl;

    int centric_count = 0;
    int centric_changed = 0;

    clipper::HKL_data_base::HKL_reference_index ih;
    for (ih = hkl_target_eff_res.first(); !ih.last(); ih.next())
    {
      // starting_phase_set[ih].fom() = 1.0;
      if (ih.hkl_class().centric() == true)
      {
        centric_count++;
        float rand_float = uniform_distribution(MyGenerator); // MyGenerator will be static or not depending on settings.
        if (rand_float < settings.centric_random_percent && !settings.control_for_phase_copying)
        {
          trusted_phase_set[ih.hkl()].phi() = known_phase_set[ih.hkl()].phi(); // We also copy to the trusted_phase_set, to allow updates over iterates.
          starting_phase_set[ih].phi() = known_phase_set[ih.hkl()].phi();      // We copy to the starting phase set.
                                                                               // custom_flag_01[ih].flag() = 1; // Remember which ones were changed using custom flags.
          centric_changed++;
        }
      }
    }
    outlog << "Of " << centric_count << " centric reflections, " << centric_changed << " were updated to their known phase value." << std::endl;

    settings.fix_trusted_phases = settings.n_fix_centric_phases;
    outlog << "Trusted phase set object initialized with centric phases. These phases will be fixed for the first " << settings.fix_trusted_phases << "iterations. " << std::endl;
  }
  // ############################################

  initialise_starting_map(); // Now that everything is initialized we can calculate the starting xn real space map.

  // Output some Run information:
  outlog << "\nAll done with run initialization, run ID: " << my_run_id << "." << std::endl;

  // Generate the file name for the Summary File and open an outstream to write to it.
  { // We do this in a custom scope to be rid of the temp vars line and filename.
    char line[512];
    char summary_filename[64];
    snprintf(summary_filename, sizeof(summary_filename), "%ssummary_run%d.txt", settings.working_dir.c_str(), my_run_id);
    summary_file.open(summary_filename); // Remember to close this at end.
  }
  // init complete.

  return true;
}

bool ipa_worker::run_loop(std::atomic<float> &ThreadProgress, const IPA_updaterule_command my_update_rule, const IPA_constraint_command my_constraints)
{
  /*
   Here we perform the basics of a single iteration, this will soon cleanly entail the following steps:

   1a. On the very first iteration, or any iteration where we change the apodization function we need to:
   Calculate the number of iterations for the apodization step
   Calculate histograms (density value and gradient magnitude) for the reference structure, with the relevant apodization function applied to the reference data.

   1b On any iteration where we change the apodization function we need to:
   Determine if we need to regrid the maps and restructure/repopulate the hkl list
   -If needed, regrid the map and interpolate both the current iterate xn and solution estimate xn_b onto the new grid. xn is needed for the update rule, xn_b for envelope generation. Also regenerate the weights array
   -If needed, regenerate the hkl list, and repopulate with the experimental data, applying any needed apodization

   2. On every iteration we have to:
   - Check if we need to update the envelope, and do so when required (envelope is generated from solution estimate xn_b, except on the very first iterate)
   - Perform the UpdateRule and associated constraint order, taking xn and generating xn+1
   - Check if we should be calculating and reporting agreement statistics, and do so when required

   TODO: ... apodization step counter doesn't work properly when we have a final apodization limit (doesn't get reset). investigate and fix
   */

  outlog << "\n*----------------------------------------------------------*" << std::endl;
  outlog << "Begining run-loop initialization: " << std::endl;

  clock_t iterate_start_time = clock();
  check_for_apodization(total_iterate_counter);
  check_if_histogram_or_gradient_matching();

  if (apodize_on_this_step)
  {
    ThreadProgress = 3 + (float(total_iterate_counter) / float(settings.n_total_iterates));
    prepare_data_objects(apodization_step);
  }

  // Generate a envelope if requested on this iterate. Decisions about the need for envelope generation are made internal to the function.
  // Envelope generation is based on the current solution estimate xn_b
  stat_xn_radius = filter_radius_look_up_table[total_iterate_counter]; // Retrieve the filter radius that is to be used for mask generation, and store it somewhere convenient for reporting purposes
  generate_a_mask(filter_radius_look_up_table[total_iterate_counter]);

  // if (apodize_on_this_step) {ap_step_iterate_counter = 0;}

  ap_step_iterate_counter++; // We do this immediately after to try not break functionality... could just set it to 1, and do this in the previous stage along with the total iterate counter.

  outlog << "\n*----------------------------------------------------------*\n"
         << std::endl;
  outlog << "Total time since start of iterate: " << (double)(clock() - iterate_start_time) / CLOCKS_PER_SEC << " seconds\n"
         << std::endl;
  // Print some information regarding this iterate:
  beta = beta_look_up_table[total_iterate_counter]; // retrieve the value of beta  from the look up table

  if ((my_update_rule == DifMap) || (my_update_rule == ReReRe) || (my_update_rule == RevRRR)) // All update rules that require parameter beta
  {
    beta = beta_look_up_table[total_iterate_counter]; // retrieve the value of beta from the look up table
    outlog << "Beta = " << beta << " " << std::endl;
  }
  if (my_update_rule == ErrRed) // All update rules where beta is irrelevant
  {
    beta = 0; // Set to 0 for simplicity
    outlog << "Beta unspecified" << std::endl;
  }

  // outlog << "Filter Radius = " << stat_xn_radius << " Å" << std::endl; // Reported when generating a mask.
  // outlog << "Working Solvent Fraction = " << stat_xn_fraction_solvent << " " << std::endl; // Reported when generating a mask.

  outlog << "Run: " << my_run_id << std::endl;
  if (settings.apodize)
  {
    outlog << "Apodization Step: " << apodization_step << std::endl;
  }
  outlog << "IPA Iteration: " << total_iterate_counter << std::endl;
  outlog << std::endl;
  if (settings.apodize)
  {
    outlog << "This is " << ap_step_iterate_counter << " of " << n_iterations_ap_step << " iterations at this apodization step" << std::endl;
  }
  else
  {
    outlog << "This is " << ap_step_iterate_counter << " of " << n_iterations_ap_step << " iterations" << std::endl;
  }
  outlog << total_iterate_counter << " iterations have been completed so far this run" << std::endl;
  if (settings.use_test_set)
  {
    outlog << "Work/Free set statistics will be calculated for this iterate" << std::endl;
  }
  else
  {
    outlog << "Work set statistics only will be calculated for this iterate" << std::endl;
  }
  define_working_and_test_sets();

  outlog << "\nResolution information:" << std::endl;
  outlog << std::setprecision(2) << " Experimental resolution (Å): " << target_max_resn << std::endl;
  outlog << std::setprecision(2) << "    Effective resolution (Å): " << rlimit << std::endl;
  outlog << std::setprecision(2) << "      Maximum resolution (Å): " << maximum_resolution << std::endl;

  outlog << "\nUpdate rule for this iterate: ";
  outlog << util::return_name_from_enum(my_update_rule) << std::endl;
  outlog << "Real Space constraints for this iterate: ";
  outlog << util::return_name_from_enum(my_constraints) << std::endl;
  outlog << "\n*----------------------------------------------------------*\n"
         << std::endl;

  // All solutions become invalid when we execute the update rule. Set the appropriate flags
  xn_a_is_valid = false;
  xn_b_is_valid = false;
  xn_a_fp_is_valid = false;
  xn_b_fp_is_valid = false;
  xn_fp_is_valid = false;

  // Update the fractional part of Threadprogress...
  ThreadProgress = (float(total_iterate_counter) / float(settings.n_total_iterates));

  // Perform the correct Update Rule, pass the correct constraint to that rule.
  if (my_update_rule == DifMap)
  {
    ThreadProgress = ThreadProgress + 4.0; // Int code for DifMap is 4.
    difference_map_algorithm(my_constraints);
  }

  if (my_update_rule == ErrRed)
  {
    ThreadProgress = ThreadProgress + 6.0; // Int code for ErrRed is 6.
    error_reduction_algorithm(my_constraints);
  }

  if (my_update_rule == ReReRe)
  {
    ThreadProgress = ThreadProgress + 5.0; // Int code for ReReRe is 5.
    relax_relax_reflect_algorithm(my_constraints);
  }

  if (my_update_rule == RevRRR)
  {
    ThreadProgress = ThreadProgress + 7.0; // Int code for RevRRR is 7.
    reversed_relax_relax_reflect_algorithm(my_constraints);
  }

  // Check the concistency of the phase sets if required for predicting IPA success:
  if (total_iterate_counter % settings.phase_consistency_interval == 0)
  {
    // Quickly validate xn_a_fp
    if (!xn_a_fp_is_valid)
    {
      xn_a.fft_to(xn_a_fp);
      xn_a_fp_is_valid = true;
    }
    // Update our phase consistency statistic
    update_phase_consistency(xn_a_fp, apodization_weights);
    outlog << "The phase consistency over the previous " << settings.phase_consistency_interval << " iterates has been calculated as: " << stat_phase_consistency << std::endl;
  }

  // If we are evaluating and writing statistics this iterate, do it
  write_stats_this_iterate = ((total_iterate_counter % settings.report_statistics_interval) == 0);
  if (write_stats_this_iterate)
  {
    outlog << "\nStatistics will be calculated for this iterate..." << std::endl;
    ThreadProgress = 8.0 + (float(total_iterate_counter) / float(settings.n_total_iterates));
    // Calculate and write iterate statistics to the log file.
    summarise_iterate_statistics(summary_file);
  }

  // DEBUG OUTPUT:
  // outlog << "FINAL ITERATE FINGERPRINT: " << std::endl;
  // write_fingerprint_to_log(xn_fp, xn_a_fp, xn_b_fp, apodized_f_work_set, xn, xn_a, xn_b, outlog);

  return true;
}

bool ipa_worker::finish_run()
{
  /*
   Here we finish up any loose ends regarding the run, closing summary log files etc, and dealing with the output of statistics and output of the final map and phases for the next step in the program, whether that be phase consensus, or envelope consensus...

   The final masks/phases are stored in the x_determination/outputs/ folder, and are labelled according to the job_name and run number. We also add the names of this file to a file_list.txt that resides in the same location. These are important for consensus jobs to automtically run at the next step of the program.

   */

  char io_map_filename[256];

  outlog << "\nThis is the end of the update loop" << std::endl;
  outlog << "Total iterations performed: " << total_iterate_counter << std::endl;
  outlog << "\nCalculating final statistics" << std::endl;

  summarise_iterate_statistics(summary_file);

  if ((settings.job_is_envelope_determination) || (settings.output_envelopes_and_phases)) // always output binary envelopes to the outputs folder during envelope determination, optionally output phases
  {
    clipper::CCP4MAPfile map_out;
    snprintf(io_map_filename, sizeof(io_map_filename), "%s%s_mask_run%d.ccp4", settings.working_dir.c_str(), settings.job_name.c_str(), my_run_id);
    map_out.open_write(io_map_filename);
    map_out.export_xmap(current_mask); // write mask - pretty sure current mask is fine from the final iterate.
    map_out.close_write();

    // Add the name of the output to the filelist.
    std::ofstream filelist;
    filelist.open(settings.working_dir + "mask_list.txt", std::ios_base::app);
    if (filelist.is_open())
    {
      filelist << io_map_filename << std::endl;
      outlog << "File written to " << settings.working_dir << "mask_list.txt" << std::endl;
    }
    filelist.close(); // Close read quickly, technically this is NOT THREAD SAFE. But MacOS is the worst for doing things re directories with standard C++ libraries, so here we are, just hope that none of the threads finish within the same couple of milliseconds...
  }

  if ((settings.job_is_phase_determination) || (settings.output_envelopes_and_phases)) // always output phases sets to the outputs folder during phase determination, optionally output binary envelopes
  {
    char io_mtz_filename[256];
    clipper::CCP4MTZfile mtzout;

    //           clipper::String output_col_fp         = "/*/*/PB_PA_xn";

    clipper::String output_col_fp = "/*/*/[F,PHIC]";
    snprintf(io_mtz_filename, sizeof(io_mtz_filename), "%s%s_phase_run%d.mtz", settings.working_dir.c_str(), settings.job_name.c_str(), my_run_id);
    mtzout.set_column_label_mode(clipper::CCP4MTZfile::Legacy);
    mtzout.open_write(io_mtz_filename);
    mtzout.export_hkl_info(hkl_target_eff_res);

    if (!xn_b_fp_is_valid)
    {
      xn_b.fft_to(xn_b_fp); // Quickly calculate this if we haven't already...
      xn_b_fp_is_valid = true;
    }

    mtzout.export_hkl_data(xn_b_fp, output_col_fp); // Pretty sure xn_b_fp is also appropriate here, as even for Error Reduction, xn() is just copied from xn_b, and xn_b will always be the most recent fourier projection.
    mtzout.close_write();

    // Add the name of the output to the filelist.
    std::ofstream filelist;
    filelist.open(settings.working_dir + "phase_list.txt", std::ios_base::app);
    if (filelist.is_open())
    {
      filelist << io_mtz_filename << std::endl;
      std::cout << "File written to " << settings.working_dir << "phase_list.txt" << std::endl;
    }
    filelist.close(); // Close read quickly, technically this is NOT THREAD SAFE. But MacOS is the worst for doing things re directories with standard C++ libraries, so here we are, just hope that none of the threads finish within the same couple of milliseconds...
  }

  outlog << "\nThis is the end of the Run." << std::endl;

  // Close any open output files...
  summary_file.close();
  outlog.close();

  delete[] density_histogram.interval;
  delete[] density_histogram_reference_raw.interval;
  delete[] density_histogram_reference.interval;
  delete[] density_histogram_reference_rescaled.interval;

  delete[] gradient_histogram.interval;
  delete[] gradient_histogram_reference_raw.interval;
  delete[] gradient_histogram_reference.interval;

  // Successful termination.
  outlog << "\nSuccessful Termination." << std::endl;
  return true;
}

// ------------------------------------------------------- //
// ------------------ INITIALISATION --------------------- //
// ------------------------------------------------------- //

bool ipa_worker::initialise_data_objects()
{
  /*
   Former approach:
   This function imports the measured diffraction data (|F|, sig|F|, stored in a mtz file)
   This is the primary input from the user.
   The data is read into the appropriate containers - referenced as local variables within the worker instance - for easy access within all the worker classes routines.
   Each worker has its own copy of the data in memory, and independently works from this baseline within each thread ...
   The import of the mtz may not be thread safe if multiple objects are all opening and closing the file at once ... this needs investigation

   New approach:
   measured_f_all, and mean_I_wilson along with associated hkl_info objects are passed in via the settings.
   Within this function we now copy these rather than re-read the data from the mtz.
   We now have two copies of the primary data (measured_f_all), one in the settings. and one as a member variable,
   In the future it might just be easier to set measured_f_all as a reference, and reference the one in settings...

   */

  outlog << "\nTransferring |F|, sig|F|  into the appropriate objects ... " << std::endl;

  // Copying hkl info objects is easy.
  hkl_target_exp_res = settings.hkl_target; // Constant, Experimental resolution limit
  hkl_target_eff_res = settings.hkl_target; // Variable, Effective resolution limit
  hkl_target_max_res = settings.hkl_wilson; // Constant, Maxmimum resolution limit, may differ from experimental limit

  // Initialise Data Objects to the hkl_target_eff_res info object.
  measured_f.init(hkl_target_exp_res, hkl_target_exp_res.cell()); // NOTE: INPUT DATA HKL

  test_set_flags.init(hkl_target_exp_res, hkl_target_exp_res.cell()); // NOTE: INPUT DATA HKL

  resolution_cutoff_flags.init(hkl_target_exp_res, hkl_target_exp_res.cell()); // Low resolution cut off.

  known_phase_set.init(hkl_target_exp_res, hkl_target_eff_res.cell()); // NOTE: INPUT DATA HKL

  trusted_phase_set.init(hkl_target_exp_res, hkl_target_exp_res.cell());

  // apodized_f_all.init(hkl_target_exp_res,hkl_target_exp_res.cell()); // NOTE: INPUT DATA HKL
  // measured_f.init(hkl_target_eff_res,hkl_target_exp_res.cell());
  // apodized_f.init(hkl_target_eff_res,hkl_target_exp_res.cell());

  working_f.init(hkl_target_eff_res, hkl_target_exp_res.cell());
  test_set_f.init(hkl_target_eff_res, hkl_target_exp_res.cell());
  work_set_f.init(hkl_target_eff_res, hkl_target_exp_res.cell());

  apodization_weights.init(hkl_target_eff_res, hkl_target_eff_res.cell());
  starting_phase_set.init(hkl_target_eff_res, hkl_target_eff_res.cell());

  xn_fp.init(hkl_target_eff_res, hkl_target_eff_res.cell());
  xn_a_fp.init(hkl_target_eff_res, hkl_target_eff_res.cell());
  xn_b_fp.init(hkl_target_eff_res, hkl_target_eff_res.cell());

  // measured_f_work_set.init(hkl_target_exp_res,hkl_target_eff_res.cell()); // NOTE: INPUT DATA HKL
  // apodized_f_work_set.init(hkl_target_exp_res,hkl_target_eff_res.cell()); // NOTE: INPUT DATA HKL
  // measured_f_test_set.init(hkl_target_exp_res,hkl_target_eff_res.cell()); // NOTE: INPUT DATA HKL
  // apodized_f_test_set.init(hkl_target_exp_res,hkl_target_eff_res.cell()); // NOTE: INPUT DATA HKL

  mean_i_wilson.init(hkl_target_max_res, hkl_target_max_res.cell());
  working_mean_i_wilson.init(hkl_target_eff_res, hkl_target_eff_res.cell());

  // Copy the data from the settings into the local data object.
  copy_F_sigF_data(settings.measured_f_all, measured_f);

  // Copy mean_I_wilson from the settings into the local data object.
  for (clipper::HKL_data_base::HKL_reference_index ih = mean_i_wilson.first(); !ih.last(); ih.next())
  {
    mean_i_wilson[ih] = settings.mean_i_wilson[ih];
  }

  if (settings.input_starting_phases)
  {
    // Read starting phases, if given
    if (settings.verbose)
    {
      outlog << "Importing starting phase set" << std::endl;
    }

    clipper::CCP4MTZfile mtzin;                                        // create a CCP4MTZfile object
    mtzin.open_read(settings.input_filename_target_data);              // open the mtz file for reading
    mtzin.import_hkl_data(starting_phase_set, settings.target_col_pw); // Here are our starting phases if they are given. Any input weights are ignored - within the program these weights are manipulated to effect apodization of the data - they are not used to represent the reliability of the phase estimates
    mtzin.close_read();                                                // close the mtz file after reading.
  }

  if (settings.input_known_phases)
  { // Read known phases, if they are given ... Test cases only
    if (settings.verbose)
    {
      outlog << "Importing known phase set" << std::endl;
    }

    clipper::HKL_data<clipper::data64::F_phi> pw_known(hkl_target_exp_res); // temporary object to hold Phi fom, we'll copy the phi into calculated_phase_set

    clipper::CCP4MTZfile mtz_knownphases;                               // create a CCP4MTZfile object
    mtz_knownphases.open_read(settings.input_filename_known_phase_set); // open the mtz file for reading
    mtz_knownphases.import_hkl_data(pw_known, settings.known_phase_set_fp_column_labels);
    mtz_knownphases.close_read(); // close the mtz file after reading.

    // Copy the phases, Ignore fom from input.
    clipper::HKL_data_base::HKL_reference_index ih;
    for (ih = hkl_target_exp_res.first(); !ih.last(); ih.next())
    {
      known_phase_set[ih].phi() = pw_known[ih].phi();
      known_phase_set[ih].fom() = 1.0;
    }
  }

  // Store information about space group chirality in member variable.
  space_group_is_chiral = util::space_group_is_chiral(hkl_target_exp_res.spacegroup().spacegroup_number());

  // the minimum and maximum resolution of the input data set
  target_min_resn = 1.0 / sqrt(measured_f.invresolsq_range().min());
  target_max_resn = 1.0 / sqrt(measured_f.invresolsq_range().max());

  // the maximum resolution of the model based estimates for <I>  ... these may extend to higher resolution than the data set itself, enabling the free lunch algorithm
  wilson_max_resn = 1.0 / sqrt(mean_i_wilson.invresolsq_range().max());

  // Report some basic information for the target data set

  outlog << "\nTARGET STRUCTURE: " << std::endl;
  report_hkl_info_and_data(hkl_target_exp_res, measured_f, outlog);

  // TODO: ... Should we check for non-standard space group settings which are unsupported by certain library functions ?
  // e.g. space group P 21 where the unique axis is not b  etc etc
  // Or can we assume that no modern data processing program would do that to us?
  // We should certainly check for the use of rhombohedral axes for space groups 146 (H3/R3) and 155 (H32/R32) which is similarly unsupported.

  return true;
}

bool ipa_worker::read_reference_mtz()
{
  /*
   This function reads the appropriate reference .mtz and creates the appropriate maps.
   Should be called on the reference mtz ONLY.
   */

  // Build filenames for the structure factors and the mask
  src_code_location = std::getenv("IPA_DIR"); // NOTE: THIS REQUIRES THE BASH SCRIPT TO MAKE RUN PROPERLY.
  ref_mtz_filename = src_code_location + "/reference_structures/" + reference_structures[ref_id].pdb_id + "-sigmaa.mtz";
  ref_msk_filename = src_code_location + "/reference_structures/" + reference_structures[ref_id].pdb_id + "-mask.ccp4";
  // The labels for Fourier amplitude and phase are always the same in these files - These coefficients produce a Sigmaa-Weighted 2mFo-DFc synthesis ...

  clipper::String reference_col_fp = "/*/*/[FWT,PHIC]";

  // Get the data from the reference mtzfile
  clipper::CCP4MTZfile mtzin2; // create a CCP4MTZfile object

  mtzin2.open_read(ref_mtz_filename);    // open the mtz file for reading
  mtzin2.import_hkl_info(hkl_reference); // read sg, cell, reso

  measured_fp_reference.init(hkl_reference, hkl_reference.cell()); // populate a data object of type  F + phi to hold the Fourier coefficients for the reference structure map
  apodized_fp_reference.init(hkl_reference, hkl_reference.cell()); // populate a data object of type  F + phi to hold the weighted Fourier coefficients for the reference structure map

  mtzin2.import_hkl_data(measured_fp_reference, reference_col_fp); // Here are the Fourier coefficients for the reference structure - this set will be left alone
  mtzin2.import_hkl_data(apodized_fp_reference, reference_col_fp); // Here are the Fourier coefficients for the reference structure - this set will be modified as required by the apodization scheme

  mtzin2.close_read(); // close the mtz file after reading.

  // Output some info to the log
  outlog << "\nREFERENCE STRUCTURE: " << std::endl;
  report_hkl_info_and_data(hkl_reference, measured_fp_reference, outlog);

  return true;
}

bool ipa_worker::select_reference_structure()
{
  // Select the reference structure for histogram matching
  // We look at data sets of approximately the same resolution, and select the instance with closest overall isotropic B-factor
  // To make the needed comparisons we convert resolution 1/|s| into the magnitude of the scattering vector |s|.
  // The numerical identifier for the reference structure is stored in the ref_id member variable.

  float max_magnitude_s_diff = 0.06;
  int n_accepted = 0;
  float b_difference = 10000.0;

  for (int i = 0; i < n_reference_structures; i++) // Could change this to reference_structures.size() to minimise the possibility of future error. Be mindful of off by 1 errors.
  {
    if (std::abs((1.0 / reference_structures[i].resolution_limit) - (1.0 / target_max_resn)) < max_magnitude_s_diff)
    {
      n_accepted += 1;
      if (std::abs(reference_structures[i].mean_bfactor - settings.target_b_factor_data.B_ISO) < b_difference)
      {
        ref_id = i;
        b_difference = std::abs(reference_structures[i].mean_bfactor - settings.target_b_factor_data.B_ISO);
      }
    }
  }

  // Uh oh - No reference structure was found!
  if (n_accepted == 0)
  {
    std::cerr << "No reference structures that are acceptably close in resolution - stopping now";
    return false;
  }

  outlog << "\nA reference structure (" << ref_id << ") is selected for histogram matching" << std::endl;
  outlog << "PDB ID: " << reference_structures[ref_id].pdb_id << std::endl;
  outlog << "Overall Isotropic B-factor (Å^2): " << reference_structures[ref_id].mean_bfactor << std::endl;
  outlog << "Resolution Limit (Angstroms): " << reference_structures[ref_id].resolution_limit << std::endl;
  outlog << "\nOn the basis of the following values for the target structure" << std::endl;
  outlog << "Overall Isotropic B-factor (Å^2): " << settings.target_b_factor_data.B_ISO << std::endl;
  outlog << "Resolution Limit (Å): " << target_max_resn << std::endl;

  if (!read_reference_mtz()) // Import the chosen reference structure
  {
    return false;
  }

  // Now that we have selected and imported the reference mtz, we can initialise the grid and map objects...,
  grid_reference.init(hkl_reference.spacegroup(), hkl_reference.cell(), hkl_reference.resolution(), 2.0); // Initialise grid object.

  outlog << "\n MAP INFO:\n\n";
  outlog << " Grid Dimensions: nu " << grid_reference.nu() << " nv " << grid_reference.nv() << " nw " << grid_reference.nw() << std::endl;
  outlog << " Total Size of the Grid array: " << grid_reference.size() << std::endl;

  reference_map.init(hkl_reference.spacegroup(), hkl_reference.cell(), grid_reference);
  reference_mask.init(hkl_reference.spacegroup(), hkl_reference.cell(), grid_reference);
  reference_weights.init(hkl_reference.spacegroup(), hkl_reference.cell(), grid_reference);
  reference_map_rescaled.init(hkl_reference.spacegroup(), hkl_reference.cell(), grid_reference);

  // Store the statistical weighting factors for the reference structure in a map object, as interactive lookup is very slow, and we don't want to do that repeatedly

  clipper::Xmap_base::Map_reference_index ix;
  for (ix = reference_map.first(); !ix.last(); ix.next())
  {
    reference_weights[ix] = 1.0 / reference_map.multiplicity(ix.coord());
  }

  // Import the reference mask and linearly interpolate onto the current grid

  if (safely_import_mask(input_reference_mask, ref_msk_filename, hkl_reference))
  {
    outlog << "\nPerforming linear interpolation of the input reference mask onto the current grid" << std::endl;
    for (ix = reference_mask.first(); !ix.last(); ix.next())
    {
      clipper::Coord_frac cf = ix.coord().coord_frac(grid_reference);
      reference_mask[ix] = lroundf(input_reference_mask.interp<clipper::Interp_linear>(cf)); // Round result to the nearest integer
    }
  }
  else
  {
    outlog << "Failed to linearly interpolate input reference mask onto the current grid." << std::endl;
    return false;
  }

  return true; // Successfuly selected and imported data associated with the reference structure.
}

void ipa_worker::populate_reference_data()
{
  // Detail the reference structures, this is all hard coded for now, but could be changed to read out of a folder or similar in the future if we need to dynamically add more reference data.

  // initialize an appropriately sized vector to hold the information about the reference structures

  n_reference_structures = 11;
  reference_structures.resize(n_reference_structures); // reference_structures is a member variable.

  // and here is the information ...

  reference_structures[0].pdb_id = "4kyc";
  reference_structures[0].mean_bfactor = 0.2773E+02;
  reference_structures[0].resolution_limit = 1.95;

  reference_structures[1].pdb_id = "1xf9";
  reference_structures[1].mean_bfactor = 0.4927E+02;
  reference_structures[1].resolution_limit = 2.70;

  reference_structures[2].pdb_id = "5wor";
  reference_structures[2].mean_bfactor = 0.6356E+02;
  reference_structures[2].resolution_limit = 2.77;

  reference_structures[3].pdb_id = "5g1u";
  reference_structures[3].mean_bfactor = 0.4496E+02;
  reference_structures[3].resolution_limit = 2.57;

  reference_structures[4].pdb_id = "4d82";
  reference_structures[4].mean_bfactor = 0.5527E+02;
  reference_structures[4].resolution_limit = 3.20;

  reference_structures[5].pdb_id = "5foy";
  reference_structures[5].mean_bfactor = 0.3994E+02;
  reference_structures[5].resolution_limit = 2.25;

  reference_structures[6].pdb_id = "4d4k";
  reference_structures[6].mean_bfactor = 0.6895E+02;
  reference_structures[6].resolution_limit = 3.24;

  reference_structures[7].pdb_id = "4gex";
  reference_structures[7].mean_bfactor = 0.6634E+02;
  reference_structures[7].resolution_limit = 2.80;

  reference_structures[8].pdb_id = "4c0p";
  reference_structures[8].mean_bfactor = 0.8461E+02;
  reference_structures[8].resolution_limit = 2.95;

  reference_structures[9].pdb_id = "3zxu";
  reference_structures[9].mean_bfactor = 0.1900E+03;
  reference_structures[9].resolution_limit = 3.70;

  reference_structures[10].pdb_id = "4bpx";
  reference_structures[10].mean_bfactor = 0.2495E+03;
  reference_structures[10].resolution_limit = 3.40;
}

void ipa_worker::define_working_resolution(float low_resolution_cutoff, float high_resolution_cutoff)
{
  /*
   This function does two things:
  1 ) Alters hkl_target_eff_res so that it contains only the observations inside the high resolution cutoff.
  2 ) Populates a data object that will enable subsequent flagging of data inside the low resolution cutoff as missing.
  */

  // Output some information regarding the dataset and the resolution cut-off to be applied
  outlog << "\nThe working resolution is changing:" << std::endl;
  outlog << " - A low resolution data cutoff of " << std::setprecision(2) << low_resolution_cutoff << " Å will be applied to the data." << std::endl;
  outlog << " - A high resolution data cutoff of " << std::setprecision(2) << high_resolution_cutoff << " Å will be applied to the data." << std::endl;
  outlog << "Data inside the low resolution cutoff are treated as missing." << std::endl;

  // Alter the resolution of the hkl_list, propagate to data children.

  // outlog << " Number of data at previous resolution: " << hkl_target_eff_res.num_reflections() << "." << std::endl;

  // hkl_target_eff_res.init(hkl_target_exp_res.spacegroup(), hkl_target_exp_res.cell(), clipper::Resolution((double)high_resolution_cutoff));
  hkl_target_eff_res = clipper::HKL_info(hkl_target_exp_res.spacegroup(), hkl_target_exp_res.cell(), clipper::Resolution((double)high_resolution_cutoff));
  hkl_target_eff_res.generate_hkl_list(); // Regenerates the hkl list.

  // Generate a mask to be used subsequently to apply the low resolution cut-off.
  clipper::HKL_data_base::HKL_reference_index ih;
  resolution_cutoff_flags.init(hkl_target_eff_res, hkl_target_eff_res.cell());

  for (ih = hkl_target_eff_res.first(); !ih.last(); ih.next())
  {
    if (1.0 / sqrt(ih.invresolsq()) < low_resolution_cutoff)
    {
      resolution_cutoff_flags[ih].flag() = 1; // the datum is not subject to a low resolution cutoff.
    }
  }

  // outlog << " Number of data at new resolution: " << hkl_target_eff_res.num_reflections() << "." << std::endl;
}

void ipa_worker::rebuild_hkl_target_objects()
{
  /*
   This is called in tandem with regrid working maps, and regenerates all hkl_target_eff_res data objects.
   */

  // These are simply reinitialised. The information in them is never carried over between iterates.
  xn_fp.init(hkl_target_eff_res, hkl_target_eff_res.cell());
  xn_a_fp.init(hkl_target_eff_res, hkl_target_eff_res.cell());
  xn_b_fp.init(hkl_target_eff_res, hkl_target_eff_res.cell());

  xn_fp_is_valid = false;
  xn_a_fp_is_valid = false;
  xn_b_fp_is_valid = false;

  working_f.init(hkl_target_eff_res, hkl_target_eff_res.cell());
  test_set_f.init(hkl_target_eff_res, hkl_target_eff_res.cell());
  work_set_f.init(hkl_target_eff_res, hkl_target_eff_res.cell());

  working_mean_i_wilson.init(hkl_target_eff_res, hkl_target_eff_res.cell());

  outlog << std::fixed << std::setprecision(5);

  outlog << " Number of possible data inside the high-resolution cutoff: " << hkl_target_eff_res.num_reflections() << "." << std::endl;

  // Regenerate all the data objects  by accessing the original data and associated Wilson model estimates for <I> and applying the apodization function
  clipper::HKL_data_base::HKL_reference_index ih;
  int counter = 0;
  for (ih = hkl_target_eff_res.first(); !ih.last(); ih.next())
  {
    // Note that considering data objects on the RHS below
    //
    //           apodization weights has hkl_info hkl_target_eff_res
    //           measured_f          has hkl_info hkl_target_exp_res = settings.hkl_target
    //           mean_i_wilson       has hkl_info hkl_target_max_res = settings.hkl_wilson
    //
    // Hence we need to access data in measured_f and mean_i_wilson via hkl rather than the position in the reflection list

    working_f[ih].f() = measured_f[ih.hkl()].f() * apodization_weights[ih].fom();
    working_f[ih].sigf() = measured_f[ih.hkl()].sigf() * apodization_weights[ih].fom();

    test_set_f[ih].f() = measured_f[ih.hkl()].f() * apodization_weights[ih].fom();
    test_set_f[ih].sigf() = measured_f[ih.hkl()].sigf() * apodization_weights[ih].fom();

    work_set_f[ih].f() = measured_f[ih.hkl()].f() * apodization_weights[ih].fom();
    work_set_f[ih].sigf() = measured_f[ih.hkl()].sigf() * apodization_weights[ih].fom();

    // Intensities must be multiplied by the square of the apodization function
    double w_squared = apodization_weights[ih].fom() * apodization_weights[ih].fom();

    working_mean_i_wilson[ih].I() = mean_i_wilson[ih.hkl()].I() * w_squared;
    working_mean_i_wilson[ih].sigI() = mean_i_wilson[ih.hkl()].sigI() * w_squared;
  }

  outlog << " Number of observed data before application of low-resolution cutoff : " << working_f.num_obs() << "." << std::endl;

  // Now apply the low resolution data cutoff - data falling below the cutoff get set to missing.
  // Note that all data objects used here have hkl_info hkl_target_eff_res

  working_f.mask(resolution_cutoff_flags == 1);
  test_set_f.mask(resolution_cutoff_flags == 1);
  work_set_f.mask(resolution_cutoff_flags == 1);

  outlog << " Number of observed data after application of low-resolution cutoff : " << working_f.num_obs() << "." << std::endl;

  outlog << "Intensity Statistics for the data following Apodization: " << std::endl;

  util::print_intensity_statistics_data_vs_model(working_f, working_mean_i_wilson, settings.res_bins, outlog);
}

bool ipa_worker::define_reference_density_histograms()
{
  // holds the raw reference density histogram
  density_histogram_reference_raw.number_of_intervals = settings.n_density_bins;
  density_histogram_reference_raw.interval = new float[density_histogram_reference_raw.number_of_intervals];

  // holds the smoothed and renormalized reference density histogram
  // histogram density_histogram_reference; // moved to member
  density_histogram_reference.number_of_intervals = settings.n_density_bins;
  density_histogram_reference.interval = new float[density_histogram_reference.number_of_intervals];

  // holds a rescaled reference density histogram
  // histogram density_histogram_reference_rescaled;
  density_histogram_reference_rescaled.number_of_intervals = settings.n_density_bins;
  density_histogram_reference_rescaled.interval = new float[density_histogram_reference_rescaled.number_of_intervals];

  // Set up the gradient magnitude histograms

  // holds the raw reference gradient magnitude histogram (not currently in use)
  // histogram gradient_histogram_reference_raw;
  gradient_histogram_reference_raw.number_of_intervals = settings.n_gradient_bins;
  gradient_histogram_reference_raw.interval = new float[gradient_histogram_reference_raw.number_of_intervals];

  // holds the reference gradient magnitude histogram
  // histogram gradient_histogram_reference; //moved to member
  gradient_histogram_reference.number_of_intervals = settings.n_gradient_bins;
  gradient_histogram_reference.interval = new float[gradient_histogram_reference.number_of_intervals];

  return true;
}

bool ipa_worker::safely_import_mask(clipper::Xmap<float> &mask, const std::string &mask_file, const clipper::HKL_info &hkl_check)
{
  /*
   Imports a real-valued mask from a file, storing it in the map object passed by reference.
   Also checks if the imported mask has appropriate symmetry.
   Used to import both the known mask if it is known for testing, but also a starting mask, which may be obtained from previous runs or experimentally. i.e. a working mask.

   Note the hkl object passed is only to check space group number, thus resolution independant, making hkl_target_eff_res and hkl_target_max_res and hkl_target_exp_res the same usability here.
   */
  clipper::CCP4MAPfile mapin1;
  mapin1.open_read(mask_file);
  mapin1.import_xmap(mask);
  mapin1.close_read();

  // Check that the input mask has the correct symmetry
  // TODO ... Check the cell dimensions and extent of the mask for safety ...

  outlog << " Grid Dimensions of the mask: nu " << mask.grid_sampling().nu() << " nv " << mask.grid_sampling().nv() << " nw " << mask.grid_sampling().nw() << std::endl;
  outlog << " Symmetry of the mask: " << mask.spacegroup().symbol_hm() << std::endl;

  if (mask.spacegroup().spacegroup_number() != hkl_check.spacegroup().spacegroup_number())
  {
    outlog << "\n Space group symmetry of the mask '" << mask_file << "' and the target data are different" << std::endl;
    return false;
  }
  else
  {
    outlog << "\n Looks okay. Proceeding ..." << std::endl;
    return true;
  }
}

void ipa_worker::calculate_apodization_look_up_table()
{
  /*
   Calculates the iterates at which the apodization function will change.

   Depending on the requested start, end and number of steps, it may not be possible to exactly meet the user's request

   As configured, the total number of apodization steps is weighted more important than the position of the 'last' step. This guarantees the total number of iterates.

   In the future we may consider dividing by num + 1? Thus guaranteeing that both the last iterate and total number of iterates is satisfied.

   Either way iterate is a quantized int, and thus division may lead to edge cases.

   There have been some updates in this space, the first apodisation is now a misnomer - its more of an offset - if ANY apodisation occurs, then it WILL occur on the 0th iterate, but then the next series of apodisation is even distributed between start and last iterate from the settings. This should allow for the offset behaviour we want to mimic the old parameters.
   */

  if (settings.n_apodization_steps > 0) // safety for /0
  {
    int start = settings.n_apodization_first_iterate;
    int last = settings.n_apodization_last_iterate;
    int num = settings.n_apodization_steps;

    if (start > 0) // If the first Apodisation step is delayed, we still want to apodise on the 0th iterate.
    {
      apodization_iterate_look_up_table.push_back(0);
      settings.n_apodization_steps += 1; // We increment this here, AFTER the num has been defined, thus affecting the rest of the program, but not the LUT generation.
    }

    int increments = (last - start) / (num); // Define the size of an increment as "number of iterates", We no longer add one to num here.

    if (increments == 0)
    {
      outlog << "WARNING: The requested number of Apodization steps does not fit within the defined starting and final iterates within which Apodization is requested. Thus apodization will change on EVERY iteration ... This is likely an experiment.param error input by the user..." << std::endl;
    }

    for (int i = 0; i < num; i++)
    {
      int apo_iter = start + (i * increments);
      apodization_iterate_look_up_table.push_back(apo_iter);
    }
  }

  if (settings.n_apodization_steps == 0 && settings.apodize)
  {
    // There are no apodization steps, however there has been an initial apodization call made, which means a single step is requested that never changes
    apodization_iterate_look_up_table.push_back(0); // Add a single LUT for iterate 0;
  }

  // If we are using verbose settings, we report the iterates at which apodization will change.
  if (settings.verbose)
  {
    outlog << "Total size of Apodization Iterate LUT: " << apodization_iterate_look_up_table.size() << std::endl;
    outlog << "Apodization function will change on iterates: [";
    for (int i = 0; i < apodization_iterate_look_up_table.size(); i++)
    {
      outlog << apodization_iterate_look_up_table[i];

      // Place a comma between every number, but not at the end.
      if ((i + 1) < apodization_iterate_look_up_table.size())
      {
        outlog << ", ";
      }
    }
    outlog << "]." << std::endl;
  }
}

bool ipa_worker::initialise_apodization_scheme()
{

  /*
   This function initializes the apodization scheme, based on the total number of apodization steps,
   initialising the apodization scheme entails filling up the sigma_apodize[] with sigma values,
   to be used as a LUT for each apodization step.
   */

  // The scalar defaults to 2.4477, 1/20th of the guassian maximum at half width. But can be changed in input settings.
  // double half_width_scalar = ipa_functions::calculate_scale_factor_for_half_width_at_nth_max(settings.guassian_half_width_apodization);

  sigma_apodize.resize(settings.n_apodization_steps + 1, 0.0); // Sigma_apodize is now populated as a std::vector.

  double smax = sqrt(measured_f.invresolsq_range().max());

  double sigma_initial = settings.initial_apodization_sigma; // 1.0/(half_width_scalar*settings.initial_res_at_twentieth_of_max); // Parameterise.

  sigma_apodize[0] = sigma_initial;
  double initial_integral = gaussian_integral(&sigma_initial, &smax);

  outlog << "\nApodization Scheme :\n\n";

  outlog << " Step (iterate): 0 (" << apodization_iterate_look_up_table[0] << ") "
         << " integral: " << initial_integral << " Sigma (Å^-1): " << sigma_apodize[0] << std::endl;

  if (settings.n_apodization_steps > 0) // apodization varies, work out the details
  {

    double final_integral;

    if (settings.upper_apodization_limit) // an apodization function is applied at the final step
    {
      double sigma_final = settings.final_apodization_sigma;
      sigma_apodize[settings.n_apodization_steps - 1] = sigma_final;
      final_integral = gaussian_integral(&sigma_final, &smax);
    }
    else // No apodization is applied at the final step - use full extent of the data
    {
      // Force the final step to have no apodization
      sigma_apodize[settings.n_apodization_steps - 1] = float(-1.0); // -1 is no apodization applied.
      final_integral = smax;
    }
    // outlog << "Final Integral = " << final_integral << std::endl;

    double target_integral = initial_integral;
    for (int apo_step = 1; apo_step < settings.n_apodization_steps - 1; apo_step++)
    {
      double sigma;
      target_integral += ((final_integral - initial_integral) / (double(settings.n_apodization_steps - 1)));
      get_gaussian_width(&smax, &target_integral, &sigma);
      sigma_apodize[apo_step] = (float)sigma; // use to cast float.
      outlog << " Step (iterate): " << apo_step << " (" << apodization_iterate_look_up_table[apo_step] << ") "
             << " integral: " << target_integral << " Sigma (Å^-1): " << sigma_apodize[apo_step] << std::endl;
    }

    if (settings.upper_apodization_limit) // an upper apodization limit is applied
    {
      outlog << " Step (iterate): " << (settings.n_apodization_steps - 1) << " (" << apodization_iterate_look_up_table[settings.n_apodization_steps - 1] << ") "
             << " integral: " << final_integral << " Sigma (Å^-1): " << sigma_apodize[settings.n_apodization_steps - 1] << std::endl;
    }
    else
    {
      outlog << " Step (iterate): " << (settings.n_apodization_steps - 1) << " (" << apodization_iterate_look_up_table[settings.n_apodization_steps - 1] << ") "
             << "  No apodization applied " << std::endl;
    }

    /*
     outlog << "######### DEBUG SIGMA APODIZE #########" << std::endl;
     for (int i = 0 ; i < apodization_iterate_look_up_table.size() ; i ++ )
     {
     outlog << "Sigma_Apodize[" << i << "]: " << sigma_apodize[i] << std::endl;
     }

     for (int i = 0 ; i < apodization_iterate_look_up_table.size() ; i ++ )
     {
     outlog << "ApoLUT[" << i << "]: " << apodization_iterate_look_up_table[i] << std::endl;
     }
     */
  }

  return true;
}

void ipa_worker::calculate_ipa_beta_for_iterate(int n_iterates, bool lerp)
{
  /*
   This function calculates the parameter beta for the DM, RRR, and RevRRR IPAs, based on linear interpolation of the function parameters found in the beta_params_look_up_table[]

   Beta is calculated based on one of two functions

   1. A chirped rectangular pulse wave, with the pulse duration under user control
   2. A chirped and "flattened" sinusoid, with the degree of flattening under user control

   The "chirping" of each function can be based on a linear change in frequency or a linear change in period, as specified by the bool linear_freq_change
   */

  outlog << "\nCalculating Beta for all iterates... " << std::endl;
  if (settings.verbose)
  {
    outlog << "iterate  beta" << std::endl;
  }

  beta_look_up_table.resize(n_iterates); // resize the vector containing the beta values according to the specified number of iterations.

  int number_of_parameter_sets = settings.beta_params_look_up_table.size();
  if (number_of_parameter_sets == 0)
  {
    std::cerr << "The parameter array has not been populated, call parse_ipa_beta_settings_file(), and double check the file beta.params is populated";
  }

  double integral_summation = 0;

  for (int iterate = 0; iterate < n_iterates; iterate++) // Loop over all iterates
  {

    double output = 0;

    beta_param current;

    for (int iset = 1; iset < number_of_parameter_sets; iset++)
    {
      if (iterate <= settings.beta_params_look_up_table[iset].iterate)
      {

        // Identify the relevant settings
        beta_param start = settings.beta_params_look_up_table[iset - 1];
        beta_param end = settings.beta_params_look_up_table[iset];

        current.phase = start.phase; // the phase offset for each linear segment, taken as the specified starting value

        // Calculate the LERP constant t, force floats to prevent rounding
        float t = float(iterate - start.iterate) / float(end.iterate - start.iterate);

        if (!lerp)
        {
          t = 0;
        } // Don't linearly interpolate, just use the starting settings

        // calculate the settings for everyting except the period and phase by linear interpolation, and assign to current settings param.
        // std::lerp(a,b,t) doesn't work cause !C++20, so doing lerp manually.  c = a + t(b - a)

        double newmidline = start.midline + t * (end.midline - start.midline);
        double newamplitude = start.amplitude + t * (end.amplitude - start.amplitude);
        float newfunctionparam_a = start.functionparam_a + t * (end.functionparam_a - start.functionparam_a);
        float newfunctionparam_b = start.functionparam_b + t * (end.functionparam_b - start.functionparam_b);

        // Beta parameter settings for the current iterate
        current.midline = newmidline;
        current.amplitude = newamplitude;
        current.functionparam_a = newfunctionparam_a;
        current.functionparam_b = newfunctionparam_b;

        // outlog << iterate  << std::endl;
        // outlog << current.midline << std::endl;
        // outlog << current.amplitude << std::endl;
        // outlog << current.phase << std::endl;
        // outlog << current.functionparam_a << std::endl;
        // outlog << current.functionparam_b << std::endl;

        /*
         The change in Frequency/Period needs to be handled differently to the other variables

         We wish to increase or decrease the wave frequency with x.

         In signal processing this is known as a chirp.

         Let's first suppose we want the *frequency* to change linearly with x (a linear chirp)

         i.e. f(x) = f0 + c*x where f0 is the initial frequency and c is the constant rate of frequency change (the chirp rate) (both written in terns of x)

         chirp rate = (frequncy_end - frequency_start) / length of the transition

         If we have a sinusoidal function we need to supply as an argument the *integral* of the function describing the frequency change
         ( see e.g. https://en.wikipedia.org/wiki/Chirp#Linear_chirp). This can be conceptualized as a phase correction

         φ(t) = φ(0) + 2*π*∫F(x)dx (where the integral is evaluated between 0 and x)

         This will evaluate to:

         φ(t) = φ(0) + 2*π*(f0*x + 1/2*c*x^2)

         Let's now suppose we want the *period* to change linearly with x

         i.e. T(x) = T0 + d*x  where T0 is the initial period and d is the constant rate of period change  (the chirp rate)

         chirp rate = (period_end - period_start) / length of the transition

         Of course we can also express this as a *non-linear* change in the frequency

         f(x) = 1/(T0 + d*x)

         As above we calculate the required phase correction as

         φ(t) = φ(0) + 2*π*∫F(x)dx (where the integral is evaluated between 0 and x)
         = φ(0) + 2*π*(ln(T0+d*x)/d -ln(T0)/d)
         = φ(0) + 2*π*(1/d)*ln(1+(d*x/T0)
         = φ(0) + 2*π*(1/d)*ln(1+d*x*f0)

         */

        double chirp_rate;

        // Work out the linear change in either frequency of period ... aka the chirp rate

        if (settings.linear_freq_change)
        {
          chirp_rate = double(1.0 / end.period - 1.0 / start.period) / double(end.iterate - start.iterate); // dimensionless chirp rate specifies a linear change in wave frequency
        }
        else
        {
          chirp_rate = double(end.period - start.period) / double(end.iterate - start.iterate); // dimensionless chirp rate specifies a linear change in wave period
        }

        double initial_period = start.period;
        double initial_frequency = 1 / start.period;

        double x = double(iterate - start.iterate); // the relative position for this linear segment
        double frequency;

        // Work out the current frequency, based on our relative position and the chirp rate

        if (settings.linear_freq_change)
        {
          frequency = (initial_frequency + chirp_rate * x);
        }
        else
        {
          frequency = 1.0 / (initial_period + chirp_rate * x);
        }

        // work out the integral for this linear segment, required for evaluation of the functions

        double integral;

        if (settings.linear_freq_change)
        {
          integral = x * initial_frequency + 0.5 * chirp_rate * std::pow(x, 2);
        }
        else
        {
          if (!std::isinf(1 / chirp_rate))
          {
            integral = (1 / chirp_rate) * std::log(1.0 + chirp_rate * x * initial_frequency);
          }
          else // chirp_rate must be zero and the expression above is invalid. Use the expression for the zero change case
          {
            integral = x * initial_frequency;
          }
        }

        // Add to the evaluated integrals for all previous linear segments.

        double total_integral = integral + integral_summation;

        if (iterate == settings.beta_params_look_up_table[iset].iterate) // We are at the final iterate in this particular linear segment, so add the integral to the grand total
        {
          integral_summation += integral;
        }

        // outlog << current.iterate  << " " << integral << " " << integral_summation << " " << total_integral << " " << std::endl;

        double nonscaled;

        // double phase = clipper::Util::pi(); // start at the center of a rectangular pulse / the peak of the sinusoidal function
        // double phase = clipper::Util::pi()/2; // start at the beginning of a rectangular pulse / the zero point of the sinusoidal function

        double phase = current.phase * clipper::Util::twopi();

        switch (settings.beta_function)
        {
        case 1: // square pulse train
        {
          double pulse_duration = current.functionparam_b;
          nonscaled = ipa_functions::pulse_wave(pulse_duration, phase, total_integral);
        }
        break;
        case 2: // flattened sinusoidal wave
        {
          double squareness = current.functionparam_a;
          nonscaled = ipa_functions::flattened_sinusoid(squareness, phase, total_integral);
        }
        break;
        }

        output = nonscaled * current.amplitude + current.midline; // Linear scaling of the output function

        break; // exit the loop over available parameter sets - no reason to continue, as we have found the right settings for this iterate.
      }
    }

    beta_look_up_table[iterate] = output; // add the calculated beta value to the look up table
    if (settings.verbose)
    {
      outlog << iterate << " " << beta_look_up_table[iterate] << std::endl;
    }
  }

  outlog << n_iterates << " Beta values have been calculated." << std::endl;
}

void ipa_worker::calculate_filter_radius_for_iterate(const int n_iterates, const bool lerp)
{
  outlog << "\nCalculating filter radius for all iterates... " << std::endl;
  if (settings.verbose)
  {
    outlog << "iterate  radius" << std::endl;
  }

  filter_radius_look_up_table.resize(n_iterates); // resize the vector containing the filter radii according to the specified number of iterations.

  int number_of_parameter_sets = settings.filter_radius_params_look_up_table.size();
  if (number_of_parameter_sets == 0)
  {
    std::cerr << "The parameter array has not been populated, call parse_filter_radius_settings_file(), and double check the file filter_radius.params is populated";
  }

  double filter_radius = 0.0; // Note the scope, not to be confused with settings.filter_radius;

  for (int iterate = 0; iterate < n_iterates; iterate++) // Loop over all iterates
  {
    for (int iset = 1; iset < number_of_parameter_sets; iset++) // Loop over available parameter sets
    {
      if (iterate <= settings.filter_radius_params_look_up_table[iset].iterate)
      {
        // Identify the start and end settings
        filter_param start = settings.filter_radius_params_look_up_table[iset - 1];
        filter_param end = settings.filter_radius_params_look_up_table[iset];

        // Calculate the LERP consant t, force floats to prevent rounding
        double t = double(iterate - start.iterate) / double(end.iterate - start.iterate);

        if (!lerp)
        {
          t = 0;
        } // Don't linearly interpolate, just use the starting settings

        // Linear interpolation between the start and end settings
        filter_radius = start.radius + t * (end.radius - start.radius);
        break; // exit the loop over available parameter sets - no reason to continue, as we have found the right settings for this iterate.
      }
    }

    filter_radius_look_up_table[iterate] = filter_radius; // add the calculated filter radius to the look up table
    if (settings.verbose)
    {
      outlog << iterate << " " << filter_radius_look_up_table[iterate] << std::endl;
    }
  }

  outlog << n_iterates << " Filter radii have been calculated." << std::endl;
}

void ipa_worker::calculate_solvent_fraction_multiplier_for_iterate(const int n_iterates, const bool found_solvent_fraction_multiplier_params, const bool lerp)
{

  solvent_fraction_multiplier_look_up_table.resize(n_iterates); // resize the vector containing the solvent fraction multiplier according to the specified number of iterations.

  double solvent_fraction_multiplier;

  if (found_solvent_fraction_multiplier_params)
  {
    outlog << "\nCalculating solvent fraction multiplier for all iterates... " << std::endl;
    if (settings.verbose)
    {
      outlog << "iterate  multiplier" << std::endl;
    }

    int number_of_parameter_sets = settings.solvent_fraction_multiplier_params_look_up_table.size();
    if (number_of_parameter_sets == 0)
    {
      std::cerr << "The parameter array has not been populated, call parse_solvent_fraction_multiplier_settings_file(), and double check the file solvent_fraction_multiplier.params is populated";
    }

    for (int iterate = 0; iterate < n_iterates; iterate++) // Loop over all iterates
    {
      for (int iset = 1; iset < number_of_parameter_sets; iset++) // Loop over available parameter sets
      {
        if (iterate <= settings.solvent_fraction_multiplier_params_look_up_table[iset].iterate)
        {
          // Identify the start and end settings
          solvent_fraction_param start = settings.solvent_fraction_multiplier_params_look_up_table[iset - 1];
          solvent_fraction_param end = settings.solvent_fraction_multiplier_params_look_up_table[iset];

          // Calculate the LERP consant t, force floats to prevent rounding
          double t = double(iterate - start.iterate) / double(end.iterate - start.iterate);

          if (!lerp)
          {
            t = 0;
          } // // Don't linearly interpolate, just use the starting settings

          // Linear interpolation between the start and end settings
          solvent_fraction_multiplier = start.multiplier + t * (end.multiplier - start.multiplier);
          break; // exit the loop over available parameter sets - no reason to continue, as we have found the right settings for this iterate.
        }
      }

      solvent_fraction_multiplier_look_up_table[iterate] = solvent_fraction_multiplier; // add the calculated solvent fraction multiplier to the look up table
      if (settings.verbose)
      {
        outlog << iterate << " " << solvent_fraction_multiplier_look_up_table[iterate] << std::endl;
      }
    }

    outlog << n_iterates << " Solvent fraction multipliers have been calculated." << std::endl;
  }
  else // solvent radius multiplier is not specified, so silently set all the values in the lookup table to unity
  {
    for (int iterate = 0; iterate < n_iterates; iterate++) // Loop over all iterates
    {
      solvent_fraction_multiplier_look_up_table[iterate] = 1.0;
    }
  }
}

void ipa_worker::remove_low_resolution_data(clipper::HKL_data<clipper::data64::F_sigF> &data, clipper::HKL_data<clipper::data64::Flag> &datacutoff)
{
  // TO BE MADE OBSELETE>
  /*
   Simple function which flags the very low resolution data as missing, based on a user supplied cut-off which is defined as a member variables.
   data is the unmodified input data object, datacutoff is the modified output data object
   measured_f is typically passed here, which is strictly unnecessary as it's also a member variable. However passing the dataset as an argument makes the function generally applicable.
   */

  // Output some information regarding the dataset and the resolution cut-off to be applied
  outlog << "\nA low resolution data cutoff of " << settings.low_res_cutoff << " Angstroms will be applied" << std::endl;
  int n_start = data.num_obs();
  outlog << "Number of observations before applying the cutoff: " << n_start << std::endl;

  // Generate a mask to apply the resolution cut-off to the data.
  clipper::HKL_data_base::HKL_reference_index ih;
  for (ih = hkl_target_eff_res.first(); !ih.last(); ih.next())
  {
    if (1.0 / sqrt(ih.invresolsq()) > settings.low_res_cutoff)
    {
      datacutoff[ih].flag() = 1;
    }
  }

  data.mask(datacutoff != 1); // apply the mask

  // Report some output information
  int n_end = data.num_obs();
  outlog << "Number of observations after applying the cutoff: " << n_end << std::endl;
  outlog << "Number of observations eliminated: " << n_start - n_end << std::endl;
}

bool ipa_worker::generate_a_starting_phase_set()
{
  starting_phase_set.init(hkl_target_eff_res, hkl_target_eff_res.cell());

  if (settings.fourier_space_phase_generation)
  {
    outlog << "\nGenerating a random starting phase set\n"
           << std::endl;

    // Fourier space method for phase generation - simply assign random phases for each observation
    clipper::HKL_data_base::HKL_reference_index ih;
    for (ih = hkl_target_eff_res.first(); !ih.last(); ih.next())
    {
      starting_phase_set[ih].fom() = 1.0;
      starting_phase_set[ih].phi() = 0.0;

      if (ih.hkl_class().centric())
      // centric observations ( Must be 0 or 180 degrees, i.e. 0 or pi.)
      {
        if (uniform_distribution(MyGenerator) < 0.5)
          starting_phase_set[ih].phi() = clipper::Util::pi();
      }
      else
      // acentric observations
      {
        starting_phase_set[ih].phi() = uniform_distribution(MyGenerator) * clipper::Util::twopi();
      }
    }
  }
  else
  {
    outlog << "\nGenerating a starting phase set by random positioning of spheres in the unit cell \n"
           << std::endl;

    // real space method for phase generation - based on the random positioning of appropriately sized spheres, and backtransformation to generate the starting phase set.

    // clipper::HKL_data<clipper::data64::F_phi> work_fp ( hkl_target_eff_res ); // Temporary.

    clipper::Grid_sampling grid(hkl_target_eff_res.spacegroup(), hkl_target_eff_res.cell(), hkl_target_eff_res.resolution(), 1.5); // Create grid object. Linear Shannon rate of 1.5 -> Grid Spacing of Resolution/3
    clipper::Xmap<float> xstart(hkl_target_eff_res.spacegroup(), hkl_target_eff_res.cell(), grid);                                 // Dimension a map object to store the starting electron density

    std::vector<clipper::Coord_frac> sphere_centres;
    clipper::Coord_frac cfa;
    clipper::Coord_frac cfb;
    clipper::Coord_orth coa;
    clipper::Coord_orth cob;

    float map_frac_solvent = 0.0;

    int n_attempts = 0;

    //         float lower_limit = fraction_solvent - 0.05*fraction_solvent;
    //         float upper_limit = fraction_solvent + 0.05*fraction_solvent;

    float lower_limit = solvent_fraction_multiplier_look_up_table[0] * settings.fraction_solvent * 0.95;
    float upper_limit = solvent_fraction_multiplier_look_up_table[0] * settings.fraction_solvent * 1.05;

    while ((map_frac_solvent < lower_limit) || (map_frac_solvent > upper_limit))
    {
      n_attempts += 1;

      if (n_attempts > settings.max_attempts)
      {
        std::cerr << "Too many unsuccessful attempts - Alter the sphere radius or the number of spheres" << std::endl;
        return false;
      }

      outlog << "Attempt: " << n_attempts << std::endl;

      // Randomly generate the centroids of the spheres

      sphere_centres.clear();
      for (int i = 0; i < settings.n_sphere_phase_generation; i++)
      {
        cfa = clipper::Coord_frac(uniform_distribution(MyGenerator),
                                  uniform_distribution(MyGenerator),
                                  uniform_distribution(MyGenerator));
        sphere_centres.push_back(cfa);

        outlog << "Fractional coordinates of sphere: " << i << std::endl;
        outlog << sphere_centres[i].u() << " " << sphere_centres[i].v() << " " << sphere_centres[i].w() << std::endl;

        /*             for ( int j = 1; j < hkl_target_eff_res.spacegroup().num_symops(); j++ )
         {
         cfb = cfa.transform(hkl_target_eff_res.spacegroup().symop(j));
         sphere_centres.push_back(cfb);
         }
         */
      }

      /*
       outlog << "Fractional Coordinates of symmetry equivalent sphere centroids" << std::endl;
       for ( int i = 0; i < sphere_centres.size(); i++ ) // loop over the symmetry equivalent sphere centroids << std::endl;
       {
       outlog << sphere_centres[i].u() << " " << sphere_centres[i].v() << " " << sphere_centres[i].w() << std::endl;
       }
       */

      // now figure out if we are protein or solvent
      int n_protein = 0;
      int n_total = 0;
      clipper::Xmap_base::Map_reference_index ix;
      for (ix = xstart.first(); !ix.last(); ix.next()) // loop over the asymmetric unit of the starting map
      {
        n_total += 1;
        xstart[ix] = settings.expected_mean_density_solvent;         // By default the point is solvent
        cfa = ix.coord().coord_frac(grid);                           // get the fractional coordinates of the current grid point
        coa = cfa.coord_orth(hkl_target_eff_res.cell());             // get the orthogonal coordinates of the current grid point
                                                                     //               for ( int i = 0; i < sphere_centres.size(); i++ ) // loop over the sphere centroids
        for (int i = 0; i < settings.n_sphere_phase_generation; i++) // loop over the sphere centroids
        {
          // find the position of the sphere centroid nearest to the current grid point, allowing for translational symmetry
          //                 cfb = sphere_centres[i].lattice_copy_near(cfa);
          // find the position of the sphere centroid nearest to the current grid point, allowing for symmetry
          cfb = sphere_centres[i].symmetry_copy_near(hkl_target_eff_res.spacegroup(), hkl_target_eff_res.cell(), cfa);

          cob = cfb.coord_orth(hkl_target_eff_res.cell());        // get the corresponding orthogonal coordinates of the current grid point
          float distance = clipper::Coord_orth::length(coa, cob); // compute the distance between sphere centroid and current point
          if (distance <= settings.rad_sphere_phase_generation)
          {
            xstart[ix] = settings.expected_mean_density_protein;
            n_protein += 1;
            //                      outlog << distance << std::endl;
            //                      outlog << cfa.u() << " " << cfa.v() << " " << cfa.w() << " " << coa.x() << " " << coa.y() << " " << coa.z() << " " << std::endl;
            //                      outlog << cfb.u() << " " << cfb.v() << " " << cfb.w() << " " << cob.x() << " " << cob.y() << " " << cob.z() << " " << std::endl;
            //                      outlog << "yuppers" << std::endl;
            break;
          }
        }
      }

      map_frac_solvent = 1 - (float(n_protein) / float(n_total));
      outlog << "Calculated Fraction Solvent: " << map_frac_solvent << std::endl;
    }

    // Now Fourier transform the map to generate the starting phase set
    // clipper::HKL_data<clipper::data64::F_phi> work_fp ( hkl_target_eff_res );
    xstart.fft_to(xn_b_fp); // FFT of the map to get the Crystallographic structure factors

    clipper::HKL_data_base::HKL_reference_index ih;
    for (ih = hkl_target_eff_res.first(); !ih.last(); ih.next())
    {
      starting_phase_set[ih].fom() = 1.0;
      starting_phase_set[ih].phi() = xn_b_fp[ih].phi();
    }
    /*
     clipper::CCP4MAPfile map_out;
     map_out.open_write( "starting_phases.ccp4" );
     map_out.export_xmap(xstart); // write starting map
     map_out.close_write();
     */
  }

  return true;
}

bool ipa_worker::set_starting_phase_set_from_phi(const clipper::HKL_data<clipper::data64::F_phi> input_phi)
{
  /*
   This function takes the phases from the input, and uses them to update the starting phase set.

   Assumes HKL target has been imported. ( This is not checked - but could be? )
   */

  clipper::HKL_data_base::HKL_reference_index ih;
  for (ih = hkl_target_eff_res.first(); !ih.last(); ih.next())
  {
    starting_phase_set[ih].phi() = input_phi[ih].phi();
  }

  return true;
}

// ------------------------------------------------------- //
// ---------------------- RUN LOOP ----------------------- //
// ------------------------------------------------------- //

void ipa_worker::perform_fourier_space_projection(const clipper::Xmap<float> &map, clipper::HKL_data<clipper::data64::F_phi> &fp, clipper::Xmap<float> &map_out)
{
  /*
   Performs a projection onto the constraint set in Fourier space. All Fourier space projections are performed via this function - regardless of the update rule.
   To perform the projection, we Fourier transform the input map, set the amplitudes to their measured values (accounting for data apodization), and backtransform, generating map_out

   When there is a set of trusted phases, these are treated as additional constraints in Fourier space, and projection onto these constraints proceeds analogously

   The function is passed an HKL_data object where the Fourier coefficients of the output map get stored, in case the user wishes to use these amplitudes and phases for some purpose after calling the function

   */

  map.fft_to(fp); // FFT of the input map to get the Fourier coefficients

  substitute_Fourier_amplitudes(fp,
                                work_set_f,
                                settings.p_threshold,
                                settings.replace_largeI_with_square_of_expectedF,
                                settings.intensity_multiplier,
                                settings.remove_low_res_data,
                                settings.low_res_cutoff,
                                settings.res_bins,
                                outlog,
                                settings.verbose);

  // If required for this iterate, fix those phases which are from a trusted source.
  if (settings.fix_trusted_phases > total_iterate_counter)
  {
    substitute_trusted_phases(fp);
    outlog << " Substituted trusted phases during Fourier space projection." << std::endl;
  }

  map_out.fft_from(fp); // Backtransform into the output map, completing the projection operation

  /*
   if (settings.verbose)
   {
   outlog << "\n Fourier space representation after Fourier projection - First 20 Observations:" << std::endl;
   outlog << "    n     h     k     l             |F|          Phi             |Fobs|(after apodization)"  << std::endl;

   clipper::HKL_data_base::HKL_reference_index ih;

   for ( ih = map_fp.first(); ih.index() <= 20; ih.next())
   {
   outlog << std::boolalpha << std::fixed <<
   std::setw(5) << ih.index()   << " " <<
   std::setw(5) << ih.hkl().h() << " " <<
   std::setw(5) << ih.hkl().k() << " " <<
   std::setw(5) << ih.hkl().l() << " " <<
   std::setw(20) << map_fp[ih].f() <<  " "  <<
   std::setw(10) << map_fp[ih].phi() << " " <<
   std::setw(20) << measured_f_work_set[ih].f()*apodization_weights[ih].fom() << " " << std::endl;
   }
   outlog << "\n" << std::endl;
   }*/
}

void ipa_worker::perform_real_space_projection(const clipper::Xmap<float> &map_in, clipper::Xmap<float> &map_out, const IPA_constraint_command my_constraints)
{
  /*
   Current real space constraints are solvent flatness (in the solvent region)  and histogram specification (in the protein region)

   Additional constraints can be added to the enumeration IPA_constraint_command, and passed to here and acted on accordingly.

   this will subsume all the functionality of ipa_functions::density_modification, which can become the default histogram matching with solvent flattening.
   */

  if (my_constraints == pHistogramMatch_sSmooth)
  {
    clipper::Xmap<float> map_in_copy(map_in);
    clipper::Xmap_base::Map_reference_index ix;
    for (ix = map_in.first(); !ix.last(); ix.next())
    {
      map_in_copy[ix] = map_in[ix];
    }

    double new_solvent_mean = 0;
    ipa_functions::smooth_filter_map_region(map_in_copy, current_mask, 3.0, new_solvent_mean, solvent_region_flag); // Smooth with a filter radius of 3.0 angstroms.
    ipa_functions::constant_bias_scale_map_region(map_in_copy, current_mask, new_solvent_mean, settings.expected_mean_density_solvent - new_solvent_mean, 1.0, solvent_region_flag);

    // We finish off with a hard coded old school density modification, with solvent flattening turned off.
    ipa_functions::density_modification(false, // Solvent flattenning is off.
                                        true,  // Hard coded Histogram Matching
                                        false, // Hard coded off gradient matching.
                                        settings.expected_mean_density_protein,
                                        settings.expected_mean_density_solvent,
                                        density_histogram,
                                        density_histogram_reference,
                                        protein_density_stats_reference,
                                        gradient_histogram,
                                        gradient_histogram_reference,
                                        protein_gradient_stats_reference,
                                        map_in_copy,
                                        current_mask,
                                        target_weights,
                                        hkl_target_eff_res,
                                        map_out,
                                        outlog,
                                        settings.verbose);
  }

  if (my_constraints == pHistogramMatch_sFlatten)
  {
    ipa_functions::density_modification(settings.solvent_flattening,
                                        histogram_matching_on_current_iterate,
                                        gradient_matching_on_current_iterate,
                                        settings.expected_mean_density_protein,
                                        settings.expected_mean_density_solvent,
                                        density_histogram,
                                        density_histogram_reference,
                                        protein_density_stats_reference,
                                        gradient_histogram,
                                        gradient_histogram_reference,
                                        protein_gradient_stats_reference,
                                        map_in,
                                        current_mask,
                                        target_weights,
                                        hkl_target_eff_res,
                                        map_out,
                                        outlog,
                                        settings.verbose);
  }

  if (my_constraints == pGlobicize_pHistogramMatch_sFlatten)
  {
    clipper::Xmap<float> map_in_copy(map_in);
    clipper::Xmap_base::Map_reference_index ix;
    for (ix = map_in.first(); !ix.last(); ix.next())
    {
      map_in_copy[ix] = map_in[ix];
    }

    clipper::Atom_list my_atoms;
    clipper::Atom globule_atom;
    globule_atom.set_element("Xbb");
    globule_atom.set_occupancy(1.0);
    globule_atom.set_u_aniso_orth(settings.target_b_factor_data.get_clipper_U_aniso_orth());

    for (int i = 0; i < 9999; i++) // More atoms just seems to work well..
    {
      my_atoms.push_back(clipper::Atom(globule_atom));
    }

    ipa_functions::atomize_map_region(map_in_copy, current_mask, my_atoms, hkl_target_eff_res, 0.1, 1, 2.0);

    ipa_functions::density_modification(settings.solvent_flattening,
                                        histogram_matching_on_current_iterate,
                                        gradient_matching_on_current_iterate,
                                        settings.expected_mean_density_protein,
                                        settings.expected_mean_density_solvent,
                                        density_histogram,
                                        density_histogram_reference,
                                        protein_density_stats_reference,
                                        gradient_histogram,
                                        gradient_histogram_reference,
                                        protein_gradient_stats_reference,
                                        map_in_copy,
                                        current_mask,
                                        target_weights,
                                        hkl_target_eff_res,
                                        map_out,
                                        outlog,
                                        settings.verbose);
  }
}

// This function is inline as it uses Fortran subroutines (module gamma_distribution) to calculate the Incomplete Gamma Integral, and hence can't be defined externally
inline void ipa_worker::substitute_Fourier_amplitudes(clipper::HKL_data<clipper::data64::F_phi> &fp,
                                                      const clipper::HKL_data<clipper::data64::F_sigF> &measured_f,
                                                      const double &p_threshold,
                                                      const bool &replace_largeI_with_square_of_expectedF,
                                                      const double &intensity_multiplier,
                                                      const bool &remove_low_res_data,
                                                      const double low_res_cutoff,
                                                      const int &res_bins,
                                                      std::ostream &outlog,
                                                      const bool verbose)
{
  /*
   This function performs the key operation involved in projecting onto the Fourier space constraints

   Critical Input
   A set of Fourier coefficients - both amplitudes and phases (fp),
   The measured Fourier amplitudes (measured_f) with any apodization function being used already applied.

   The function substitutes the amplitudes carried in fp with the measured Fourier amplitudes
   The modified Fourier coefficients are directly returned

   Notes on the treatment of missing amplitude data
   -------------------------------------------------

   First we convert Fourier amplitudes to intensities.
   For any missing data we check if the reconstructed intensities have a reasonable probablity of being observed, based on Wilson model statistics
   If they appear to have reasonable magnitude we leave them as is.
   If they appear to be excessively large we either
      • Replace them with the square of the expected value for |F|, multiplied by a scale factor (replace_largeI_with_square_of_expectedF = true).
      • Multiply them directly by a scale factor (replace_largeI_with_square_of_expectedF = false)
    We then convert the intensities back to amplitudes

   If we allow the reconstructed amplitudes/intensities for the missing data to freely evolve, they can take on highly unrealistic values, which can severely distort the reconstructed image

   p_threshold is the probability threshold that will trigger these actions

   Probabilities are approximated using "Wilson Statistics". That is, the intensities I should be distributed as gamma random variables with

   acentric observations
   shape parameter α = 1 and scale parameter Θ =  α x μ = μ
   shape parameter α = 1 and rate parameter β = 1/Θ = 1/μ
   shape parameter α = 1 and mean μ = α x Θ = μ

   centric observations
   shape parameter α = 1/2 and scale parameter Θ = 2μ
   shape parameter α = 1/2 and rate parameter β = 1/Θ = 1/2μ
   shape parameter α = 1/2 and mean μ = α x Θ = μ

   where μ is the expected value for I

   In shape/rate parameterization the pdf of gamma distributed random variable X is
   g(x; α,β) = β^α * x^(α-1) * exp(-βx) / Γ(α)
   Let G(y; α) = 1/Γ(α) * ∫(0 <= t <= y ) t^(α-1) * exp(-t) dt   - this is the incomplete Gamma integral
   Then for any t > 0  P(X < t) = G(βt; α)

   Subroutine gammad (y, α, ifault) will evaluate the incomplete gamma integral


   ---

   Expectation values ...

   I α |F|^2. Replacing I with the expected value of I is not the same thing as replacing |F| with the expected value of |F|, which is what we really want to do

   (For any random variable X, E[X] is not necessarily equal to E[X^2]).

   How does this work out in our setting?

   Let's suppose we have random variables X and Y with X = √Y

   If Y is Gamma distributed with shape parameter α and scale parameter Θ then X is Nakagami distributed with shape parameter m, and spread parameter Ω,

   where m =  α and Ω = αΘ

   See https://en.wikipedia.org/wiki/Nakagami_distribution

   This means the acentric |F|'s are Nakagami distributed with  m = 1   and Ω = 1   * Θ = E[I]
        while the centric  |F|'s are Nakagami distributed with  m = 1/2 and Ω = 1/2 * Θ = E[I]


   The mean of a Nakagami distributed random variable is ( Γ(m+1/2)/ Γ(m) ) * √(Ω/m)

   So for the acentric observations mean |F| = Γ(3/2)/Γ(1) * √(E[I]/1) = (1/2*√π) / 1  *  √(  E[I]) = (√π)/2  * √(E[I])
      for the  centric observations mean |F| = Γ(1)/Γ(1/2) * √(2*E[I]) =       1 / √π  *  √(2*E[I]) = (√2)/√π * √(E[I])

  Hence what we need to do is multiply √E[I] by the appropriate corrective factor, to get  E[F]

   */

  clipper::HKL_data_base::HKL_reference_index ih;

  // Perform Fourier Space Projection onto the Amplitude constraints. This always performed when the function is called

  int n_flagged = 0; // a counter for the number of unmeasured (missing) amplitudes for which the reconstructed amplitudes appear improbably large.

  for (ih = fp.first(); !ih.last(); ih.next())
  {

    if (!measured_f[ih].missing()) // |F(hkl)| is measured  ... life is simple
    {
      fp[ih].f() = measured_f[ih].f(); // reset |F(hkl)| to the measured value
    }
    else if (!((ih.hkl().h() == 0) && (ih.hkl().k() == 0) && (ih.hkl().l() == 0))) // |F(hkl)| is unmeasured and is not |F(000)|, so decide what to do
    {

      double mean_I;
      double alpha;
      double beta;

      // get the expected value for I(hkl) = |F(hkl)|^2 , and construct the parameters of the relevant gamma distribution
      // Note that we can loop over the reference index here, as fp and working_mean_I_wilson will share the hkl_target_eff_res list.
      // Note also that E[F] ≠ √ E[I], so we have to compute the appropriate corrective factors for the acentric and centric case

      mean_I = working_mean_i_wilson[ih].I(); // Take < I > from the Wilson model

      double EI_to_EF_corrective_factor;

      if (ih.hkl_class().centric()) // observation is centric
      {
        alpha = 0.5;
        beta = 1.0 / (2.0 * mean_I);
        // EI_to_EF_corrective_factor = sqrt(2.0) / sqrt(clipper::Util::pi()); // The corrective factor which converts √E[I] into E[|F|]
        EI_to_EF_corrective_factor = 2.0 / clipper::Util::pi(); // The corrective factor which converts E[I] into E[|F|]^2
      }
      else // observation is acentric
      {
        alpha = 1.0;
        beta = 1.0 / (mean_I);
        // EI_to_EF_corrective_factor = sqrt(clipper::Util::pi()) / 2.0; // The corrective factor which converts √E[I] into E[|F|]
        EI_to_EF_corrective_factor = clipper::Util::pi() / 4.0; // The corrective factor which converts E[I] into E[|F|]^2
      }

      int ifault;

      double I = pow(fp[ih].f(), 2); // Calculate the intensity I(hkl) for this observation
      double y = beta * I;           // The argument for evaluation of the incomplete gamma integral

      double pgamma = gammad(&y, &alpha, &ifault);

      if (pgamma > p_threshold) // the reconstructed value appears highly improbable
      {

        n_flagged += 1;

        // operate on the intensity estimate as requested

        double I_corrected;

        if (replace_largeI_with_square_of_expectedF) // we replace the intensity estimate with the square of the expected value for |F|, multiplied by a user-specified constant
        {
          I_corrected = EI_to_EF_corrective_factor * mean_I * intensity_multiplier;
        }
        else // we simply multiply the intensity estimate by a user-specified constant
        {
          I_corrected = I * intensity_multiplier;
        }

        fp[ih].f() = sqrt(I_corrected); // Convert the intensity into the corrected value for |F|

        if (n_flagged == 1) // First time around, output a header
        {
          outlog << "\nMissing observations with improbably large reconstructed amplitudes for which P(X > I_reconstructed) < " << 1 - p_threshold << std::endl;

          if (replace_largeI_with_square_of_expectedF)
          {
            outlog << "For all listed observations I_reconstructed is replaced with E(|F|)^2 * " << intensity_multiplier << std::endl;
          }
          else
          {
            outlog << "For all listed observations I_reconstructed is replaced with I_reconstructed * " << intensity_multiplier << std::endl;
          }

          outlog << "* = observations which fall inside the low-resolution cutoff " << std::endl;
          outlog << "\n        h      k      l     I_recons      E(I)      E(|F|)^2     P(X > I_recons)   I_corrected" << std::endl;
        }

        if (remove_low_res_data && (1.0 / sqrt(ih.invresolsq()) > low_res_cutoff))
        {
          // if we are omitting the ultralow resolution data, and this observation hkl satisfies the resolution cutoff, note this in the output
          outlog << std::fixed << " * " << std::setw(6) << ih.hkl().h() << " " << std::setw(6) << ih.hkl().k() << " " << std::setw(6) << ih.hkl().l() << " " <<
              // std::setw(6) << bin << " " <<
              std::setw(12) << I << " " << std::setw(12) << mean_I << " " << std::setw(12) << EI_to_EF_corrective_factor * mean_I << " " << std::setw(16) << 1 - pgamma << " " << std::setw(14) << I_corrected << std::endl;
        }
        else
        {
          outlog << std::fixed << "   " << std::setw(6) << ih.hkl().h() << " " << std::setw(6) << ih.hkl().k() << " " << std::setw(6) << ih.hkl().l() << " " <<
              // std::setw(6) << bin << " " <<
              std::setw(12) << I << " " << std::setw(12) << mean_I << " " << std::setw(12) << EI_to_EF_corrective_factor * mean_I << " " << std::setw(16) << 1 - pgamma << " " << std::setw(14) << I_corrected << std::endl;
        }
      }
    }
  }
}

void ipa_worker::difference_map_algorithm(const IPA_constraint_command &my_constraints)
{
  /*
   This function executes the Difference Map Algorithm on the current map (xn) producing an updated map (xn+1). It assumes that apodization and initialisation of all relevant variables has been performed.

   xn_A and xn_B, are updated to contain the two solution estimates arising from the algorithm. However they are also used as convenient storage for some calculations en route.
   */

  // Not needed anymore, hkl_target_eff_res is stored as a member variable from the get go?
  // const clipper::HKL_info hkl_target_eff_res = apodized_f_work_set.base_hkl_info();

  // STEP 1: Calculate PA(xn) store in xn_A
  perform_real_space_projection(xn, xn_a, my_constraints);
  /*
  outlog << "\nComputing PA(xn) for the Difference Map Algorithm" << std::endl;
  ipa_functions::density_modification(
                                      settings.solvent_flattening,
                                      histogram_matching_on_current_iterate,
                                      gradient_matching_on_current_iterate,
                                      settings.expected_mean_density_protein,
                                      settings.expected_mean_density_solvent,
                                      density_histogram,
                                      density_histogram_reference,
                                      protein_density_stats_reference,
                                      gradient_histogram,
                                      gradient_histogram_reference,
                                      protein_gradient_stats_reference,
                                      xn,
                                      current_mask,
                                      target_weights,
                                      hkl_target_eff_res,
                                      xn_a,
                                      outlog,
                                      settings.verbose);
  */

  // STEP 2: Calculate PB(xn) store in xn_B
  outlog << "\nComputing PB(xn) for the Difference Map Algorithm" << std::endl;
  perform_fourier_space_projection(xn, xn_fp, xn_b);

  // STEP 3: Co-calculate FA(xn) and FB(xn)
  outlog << "\nComputing Relaxed Projections FA(xn) and FB(xn) for the Difference Map Algorithm" << std::endl;

  // Create map objects to store the relaxed projections FA(xn) and FB(x)
  clipper::Xmap<float> xn_fa(xn.spacegroup(), xn.cell(), xn.grid_sampling());
  clipper::Xmap<float> xn_fb(xn.spacegroup(), xn.cell(), xn.grid_sampling());

  // Optimisation: Pre-calculate loop variable constants
  clipper::Xmap_base::Map_reference_index ix; // Iterator obj
  double beta_inverse = 1 / beta;             // Precalc (1/b)
  double beta_neg_inverse = 1 - beta_inverse; // Precalc 1-(1/b)
  double beta_pos_inverse = 1 + beta_inverse; // Precalc 1+(1/b)

  for (ix = xn_a.first(); !ix.last(); ix.next())
  {
    // Here we generate both FA(xn) and FB(xn) within the same loop.
    xn_fa[ix] = (beta_neg_inverse * xn_a[ix]) + (beta_inverse * xn[ix]); // FA(xn)
    xn_fb[ix] = (beta_pos_inverse * xn_b[ix]) - (beta_inverse * xn[ix]); // FB(xn)
  }

  // STEP 4: Calculate PA(FB(xn)) and store in xn_A, overwriting PA(xn)
  outlog << "\nComputing PA(FB(xn)) for the Difference Map Algorithm" << std::endl;
  perform_real_space_projection(xn_fb, xn_a, my_constraints);

  /*
  ipa_functions::density_modification(settings.solvent_flattening,
                                      histogram_matching_on_current_iterate,
                                      gradient_matching_on_current_iterate,
                                      settings.expected_mean_density_protein,
                                      settings.expected_mean_density_solvent,
                                      density_histogram,
                                      density_histogram_reference,
                                      protein_density_stats_reference,
                                      gradient_histogram,
                                      gradient_histogram_reference,
                                      protein_gradient_stats_reference,
                                      xn_fb,
                                      current_mask,
                                      target_weights,
                                      hkl_target_eff_res,
                                      xn_a,
                                      outlog,
                                      settings.verbose);
  */

  // xn_A solution is valid, its Fourier transform is not
  xn_a_fp_is_valid = false;
  xn_a_is_valid = true;

  // STEP 5:  Calculate PB(FA(xn)) and store in xn_B, overwriting PB(xn)
  outlog << "\nComputing PB(FA(xn)) for Difference Map Algorithm" << std::endl;
  perform_fourier_space_projection(xn_fa, xn_b_fp, xn_b);

  // xn_B solution is valid, and its Fourier transform is valid
  xn_b_fp_is_valid = true;
  xn_b_is_valid = true;

  // STEP 6: Finalise the next iterate: xn+1

  for (ix = xn_a.first(); !ix.last(); ix.next())
  {
    xn[ix] = xn[ix] + beta * (xn_a[ix] - xn_b[ix]);
  }

  xn_fp_is_valid = false; // invalid after updated iterate.

  // xn now contains the starting point for the next iterate.
  // xn_a contains the first solution estimate. (Real space constraints applied last)
  // xn_b contains the second solution estimate. (Fourier space constraints applied last)
  // xn_b_fp contains the fourier transform of xn_b.

  outlog << "Difference Map Algorithm Completed." << std::endl;
}

void ipa_worker::relax_relax_reflect_algorithm(const IPA_constraint_command &my_constraints)
{
  // const clipper::HKL_info hkl_target_eff_res = apodized_f_work_set.base_hkl_info();

  /*
   Takes the current map (xn), and applies the RRR algorithm. Two solution estimates are stored in xn_A and xn_B. The Fourier coefficients of the maps are all stored in xn_fp*** respectively. The next iterate (xn+1) is stored in xn_next.

   UPDATE RULE:
   xn+1= xn+β(PB(2PAxn−xn)−PAxn)

   SOLUTION ESTIMATES:
   xn_A = PA_xn
   xn_B = PB(2PA_xn − xn)

   */

  // STEP 1: Calculate PA(xn)
  outlog << "\nComputing PA(xn) for RRR" << std::endl;
  perform_real_space_projection(xn, xn_a, my_constraints);

  // STEP 2: Calculate 2PAxn−xn, store in xn_b
  outlog << "\nComputing (2PA(xn) - xn) for RRR" << std::endl;

  clipper::Xmap_base::Map_reference_index ix;
  for (ix = xn_a.first(); !ix.last(); ix.next())
  {
    xn_b[ix] = 2.0 * xn_a[ix] - xn[ix];
  }

  // STEP 3: Calculate PB(2PAxn−xn), by projecting Fourier constraints on itself.
  outlog << "\nComputing PB(2*PA(xn) − xn) for RRR" << std::endl;
  perform_fourier_space_projection(xn_b, xn_b_fp, xn_b);

  // STEP 4: Calculate next iterate: xn+β(PB(2PAxn−xn)−PAxn)
  outlog << "Computing the next iterate: xn+1 = xn+ Beta(PB - PA)" << std::endl;
  for (ix = xn_b.first(); !ix.last(); ix.next())
  {
    xn[ix] = xn[ix] + beta * (xn_b[ix] - xn_a[ix]); // Need to store the current beta better - or its potentially more appropriate to pass into the function...
  }

  // Set validity tags of stored data.
  xn_b_is_valid = true;
  xn_b_fp_is_valid = true;

  xn_fp_is_valid = false;

  xn_a_is_valid = true;
  xn_a_fp_is_valid = false;

  outlog << "RRR iterate completed." << std::endl;
}

void ipa_worker::reversed_relax_relax_reflect_algorithm(const IPA_constraint_command &my_constraints)
{

  // const clipper::HKL_info hkl_target_eff_res = apodized_f_work_set.base_hkl_info();

  /*
   Takes the current map (xn), and applies the RRR algorithm. Two solution estimates are stored in xn_A and xn_B. The fourier coeficients of the maps are all stored in xn_fp*** respectively. The next iterate (xn+1) is stored in xn_next.

   UPDATE RULE:
   xn+1= xn+β(PB(2PAxn−xn)−PAxn)

   SOLUTION ESTIMATES:
   xn_A = PA_xn
   xn_B = PB(2PA_xn − xn)

   */

  // STEP 1: Calculate PB(xn)
  outlog << "\nComputing PB(xn) for Reversed RRR" << std::endl;

  perform_fourier_space_projection(xn, xn_b_fp, xn_b);

  // STEP 2: Calculate 2PBxn−xn, store in xn_work
  outlog << "\nComputing (2PB(xn) - xn) for Reversed RRR" << std::endl;

  clipper::Xmap<float> xn_work(hkl_target_eff_res.spacegroup(), hkl_target_eff_res.cell(), xn.grid_sampling());

  clipper::Xmap_base::Map_reference_index ix;

  for (ix = xn_b.first(); !ix.last(); ix.next())
  {
    xn_work[ix] = 2.0 * xn_b[ix] - xn[ix];
  }

  // STEP 3: Calculate PA(2PBxn−xn) - i.e. apply real space restraints to xn_a (which contains 2PBxn-xn)
  outlog << "\nComputing PA(2*PB(xn) − xn) for Reversed RRR" << std::endl;

  perform_real_space_projection(xn_work, xn_a, my_constraints);

  // STEP 4: Calculate next iterate: xn+β(PB(2PAxn−xn)−PAxn)
  outlog << "Computing the next iterate: xn+1 = xn+ Beta(PA - PB)" << std::endl;

  for (ix = xn_b.first(); !ix.last(); ix.next())
  {
    xn[ix] = xn[ix] + beta * (xn_a[ix] - xn_b[ix]);
  }

  outlog << "Reversed RRR iterate completed." << std::endl;

  // Set validity tags of stored data.
  xn_b_is_valid = true;
  xn_b_fp_is_valid = true;
  xn_fp_is_valid = false;
  xn_a_is_valid = true;
  xn_a_fp_is_valid = false;
}

void ipa_worker::error_reduction_algorithm(const IPA_constraint_command &my_constraints)
{
  // const clipper::HKL_info hkl_target_eff_res = apodized_f_work_set.base_hkl_info();

  // clipper::HKL_data<clipper::data64::F_phi> xn_fp;

  // xn.fft_to(xn_fp);

  // outlog << "ERROR REDUCTION PRE CALCULATION FINGERPRINT:" << std::endl;
  // write_fingerprint_to_log(xn_fp, xn_a_fp, xn_b_fp, apodized_f_work_set, xn, xn_a, xn_b, starting_phase_set, outlog);

  // STEP 1: Calculate xn_PA, store in xn_a
  outlog << "\nComputing PA(xn) for Error reduction algorithm" << std::endl;
  ipa_functions::density_modification(settings.solvent_flattening,
                                      histogram_matching_on_current_iterate,
                                      gradient_matching_on_current_iterate,
                                      settings.expected_mean_density_protein,
                                      settings.expected_mean_density_solvent,
                                      density_histogram,
                                      density_histogram_reference,
                                      protein_density_stats_reference,
                                      gradient_histogram,
                                      gradient_histogram_reference,
                                      protein_gradient_stats_reference,
                                      xn,
                                      current_mask,
                                      target_weights,
                                      hkl_target_eff_res,
                                      xn_a,
                                      outlog,
                                      settings.verbose);

  // STEP 2: Calculate PB(xn_PA)
  outlog << "\nComputing PB(PA(xn)) for Error reduction algorithm" << std::endl;

  // xn_a.fft_to(xn_a_fp);
  // outlog << "ERROR REDUCTION INTERMEDIATE Fingerprint: " << std::endl;
  // write_fingerprint_to_log(xn_fp, xn_a_fp, xn_b_fp, apodized_f_work_set, xn, xn_a, xn_b, starting_phase_set, outlog);

  perform_fourier_space_projection(xn_a, xn_b_fp, xn_b);

  // The solution estimate for ER is also the next step in the algorithm, so we copy xn_b to xn ready for the next iterate.

  outlog << "\nCopying xn_b into xn for next iterate." << std::endl;
  clipper::Xmap_base::Map_reference_index ix;
  for (ix = xn_b.first(); !ix.last(); ix.next())
  {
    xn[ix] = xn_b[ix];
  }

  xn_a_is_valid = true;
  xn_a_fp_is_valid = false;
  xn_b_is_valid = true;
  xn_b_fp_is_valid = true;
  xn_fp_is_valid = false;

  // outlog << "ERROR REDUCTION POST CALCULATION FINGERPRINT:" << std::endl;
  // write_fingerprint_to_log(xn_fp, xn_a_fp, xn_b_fp, apodized_f_work_set, xn, xn_a, xn_b,starting_phase_set, outlog);
}

bool ipa_worker::define_working_and_test_sets()
{
  /*
   This function generates working and test sets from the original input data. However, it utilises a mask which is generated to determine which datasets should contain which reflections. This may not be the fasted way to perform this operation. Will look into ways to making it speedier - as it occurs every iterate.
   */

  // Place data into working and test sets:
  // current_test_set = 0; // turns out this IS a member variable...

  // Define working and test sets from the low-resolution limited data.
  if (settings.use_test_set)
  {
    outlog << "Defining the working and test sets. " << std::endl
           << std::endl;

    // update the test set ID
    current_test_set++; // Increment forward.
    if (current_test_set > (settings.n_test_sets - 1))
      current_test_set = 0; // if we have cycled through the available test sets, start again

    // Ensure all the test objects are complete
    clipper::HKL_data_base::HKL_reference_index ih;
    for (ih = hkl_target_eff_res.first(); !ih.last(); ih.next())
    {
      work_set_f[ih].f() = working_f[ih].f() * (test_set_flags[ih.hkl()].flag() == current_test_set ? 1 : 0);
      test_set_f[ih].f() = working_f[ih].f() * (test_set_flags[ih.hkl()].flag() == current_test_set ? 0 : 1);
    }

    // Apply masks to all test and work sets.
    // work_set_f.mask( test_set_flags != current_test_set );
    // test_set_f.mask( test_set_flags == current_test_set );

    outlog << "Current test set: " << current_test_set << std::endl;
    outlog << "Size of the work set: " << work_set_f.num_obs() << std::endl;
    outlog << "Size of the test set: " << test_set_f.num_obs() << std::endl;
  }
  else
  {
    // Note both work_set_f and test_set_f are reset during apodization or initialization to contain measured_f at required resolution.
    if (settings.verbose)
    {
      outlog << "No test set will be used based on settings." << std::endl;
    }
  }

  return true;
}

void ipa_worker::check_for_apodization(int iterate)
{
  /*
   This function checks if the apodization scheme changes on the supplied iterate, and toggles the boolean apodize_on_this_step accordingly.
   */

  apodize_on_this_step = false;
  if (settings.apodize == true) // Operationally, this is ALWAYS true.
  {
    for (int i = 1; i < apodization_iterate_look_up_table.size(); i++) // We start at 1, as the 0th Apodization LUT is ALWAYS present, and is dealt with during Initialisation.
    {
      if (apodization_iterate_look_up_table[i] == iterate)
      {
        // iterate found in LUT, thus apodize on this iterate.
        apodize_on_this_step = true;
        apodization_step = i;
        break;
      }
      else if (apodization_iterate_look_up_table[i] > iterate)
      {
        break; // LUT is always in order, thus stop looking if we've passed the iterate.
      }
    }
  }

  if (settings.verbose)
  {
    outlog << "Iteration : " << iterate << std::endl;
    if (apodize_on_this_step)
    {
      outlog << "Apodization will change on this iterate. " << std::endl;
      outlog << "Beginning apodization initialisation for step: " << apodization_step << std::endl;
    }
    else if (iterate == 0)
    {
      outlog << "Apodization already calculated during setup." << std::endl;
    }
    else
    {
      outlog << "Apodization is unchanged on this iterate." << std::endl;
    }
  }
}

void ipa_worker::calculate_how_many_iterations_for_this_apodization_step()
{
  /*
   This function calculates how many iterations will occur at this apodization step ... mainly for logging purposes.
   It is also important for correctly evaluating the progression through an apodization step.
   The function utilizes the apodisation LUT, calculated in the init() step.
   */

  if (settings.n_apodization_steps > 0) // Safety check, this can get called when you are not apodizing... would be nice to detangle the whole 0th iterate apodization thing...
  {

    int start_iter = apodization_iterate_look_up_table[apodization_step];
    int end_iter;

    if ((apodization_step + 1) < apodization_iterate_look_up_table.size())
    {
      // We are not on the last apodization step...
      end_iter = apodization_iterate_look_up_table[apodization_step + 1];
    }
    else
    {
      // We are on the final Apodization step.
      end_iter = settings.n_apodization_last_iterate; // We go to the end. TODO: This is not calculating correctly, bug fix why n_apo_last_iterate isn't being set properly.
    }

    // Store the total number of apodization steps in appropriate member var.
    n_iterations_ap_step = end_iter - start_iter;
  }
  else
  {
    n_iterations_ap_step = settings.n_total_iterates;
  }

  ap_step_iterate_counter = 0; // Reset the counter to 0, as why else would we be recalculating the total number if not on a new apodization step.
}

void ipa_worker::prepare_data_objects(const int &apodization_step)
{
  /*
   This functions acts to re-initialise all the appropriate data objects used by the IPA, under a new apodization / resolution scheme.
   The function is called on the 0th iterate, and any iterate where the apodization function changes

   When the effective resolution changes we need to:
   - restructure the list of observations in Fourier space, appropriate to the resolution.
   - regrid the working maps, appropriate to the resolution.
   - update the reference histogram, appropriate to the resolution.

   */

  calculate_how_many_iterations_for_this_apodization_step(); // Just a handy little subroutine to help keep track of how many iterates have occured in each apo step.

  outlog << "\nApodization step: " << apodization_step << std::endl;
  outlog << "Apodization function - sigma (1/Å): " << sigma_apodize[apodization_step] << std::endl;

  rlimit = calculate_resolution_for_apodization_step(apodization_step, settings.resolution_gaussian_width_cutoff); // Calculate the effective resolution at this step

  outlog << "Effective resolution (Å): " << std::setprecision(2) << rlimit << std::endl;
  // outlog << "Effective Resolution at 1/" << std::setprecision(0) << settings.resolution_gaussian_width_cutoff << "th of maximum height (Å): " << std::setprecision(2) << rlimit << std::endl;

  if (rlimit != prev_rlimit)
  {
    updated_resolution = true; // The effective resolution has changed.
    prev_rlimit = rlimit;
  }

  if (updated_resolution)
  {
    outlog << "\nThe effective resolution has changed since the last apodization step, redefining the list of observations in Fourier space" << std::endl;
    define_working_resolution(settings.low_res_cutoff, rlimit); // Update the working resolution for objects.
  }

  // The actual apodization and resolution changes are applied here:
  outlog << "\nCalculating Apodization Weights for sigma: " << sigma_apodize[apodization_step] << std::endl;
  define_apodization_weights(sigma_apodize[apodization_step]); // Calculate the apodization weights which will be unitary if no apodization is being applied
  outlog << "Rebuilding Fourier space data objects" << std::endl;
  rebuild_hkl_target_objects(); // Rebuild the working objects with appropriate apodization and resolution.

  if (updated_resolution)
  {
    outlog << "The effective resolution has changed since the last apodization step, regridding the maps as required!" << std::endl;
    regrid_working_maps();
  }

  // Calculate histogram || gradient matching
  outlog << "\nCalculating reference histograms: " << std::endl;
  if (updated_resolution)
  {
    calculate_reference_histograms(); // This currently performs apodization on the reference within it, we could try moving this into the standard rebuild_hkl_objects system.
  }

  updated_resolution = false; // Resolution has been updated.
}

void ipa_worker::define_apodization_weights(const double sigma)
{
  /*
   Initializes the apodization weights based on a sigma value, typically given from the apodization LUT in sigma_apodize.
   If sigma is -1, no apodization is to be applied, and weights are set to 1.0.
   If sigma is anything else, appropriate weights are calculated.
   */

  apodization_weights.init(hkl_target_eff_res, hkl_target_eff_res.cell()); // Ensure we're functioning at the correct resolution.

  clipper::HKL_data_base::HKL_reference_index ih;

  if (sigma == -1)
  {
    outlog << "No apodization function is applied to the data - apodization weights are unitary " << std::endl;

    for (ih = hkl_target_eff_res.first(); !ih.last(); ih.next())
    {
      apodization_weights[ih].fom() = 1.0; // No Apodization will be applied.
    }
    return; // Break from function.
  }

  // Else, if sigma is not -1, we can calculate appropriate weights.
  double sigma_squared = std::pow(sigma, 2);

  for (ih = hkl_target_eff_res.first(); !ih.last(); ih.next())
  {
    double s_squared = ih.invresolsq();
    double weight = exp(-s_squared / (2 * sigma_squared));
    apodization_weights[ih].fom() = weight;
  }
}

float ipa_worker::calculate_resolution_for_apodization_step(int current_apo_step, int nth_height)
{
  /*

   Calculates an effective resolution limit taking into account apodization of the data.
   The effective resolution limit is defined as the larger of:
   - the maximum resolution limit for the experiment.
   - the resolution corresponding to the half width of the Gaussian apodization function at 1/nth of maximum height.
   If no apodization is being applied, the effective resolution limit is axiomatically the maximum resolution limit for the experiment

    */
  if (current_apo_step < sigma_apodize.size()) // Safety, should always be true.
  {
    if (sigma_apodize[current_apo_step] == -1) // No apodization functiion applied to the data
    {
      outlog << "No apodization function is applied to the data " << std::endl;
      return maximum_resolution; // effective resolution is maximum resolution.
    }
    else
    {
      double half_width_scalar = util::calculate_scale_factor_for_half_width_at_nth_max(nth_height);
      outlog << "A Gaussian apodization function is applied to the data " << std::endl;
      outlog << "Resolution at 1/" << std::setprecision(0) << settings.resolution_gaussian_width_cutoff << "th of maximum height (Å): " << std::setprecision(2) << (float)(1.0 / (half_width_scalar * sigma_apodize[current_apo_step])) << std::endl;
      return std::max((float)(1.0 / (half_width_scalar * sigma_apodize[current_apo_step])), maximum_resolution); // Calculate the effective resolution
    }
  }
  outlog << "ERROR, resolution requested for non-existent apodization step!" << current_apo_step << "." << std::endl;
  throw std::invalid_argument("EXCEPTION: Apodization step for resolution calculation does not exist! Apodization scheme needs to be initialized!");
  return maximum_resolution; // This should also throw an exception, because we should never get here.
}

bool ipa_worker::apodize_data()
{
  // TO BE MADE REDUNDANT.
  // Apodize the data if requested, print out information regarding apo step.

  if (apodize_on_this_step && sigma_apodize[apodization_step] != -1) // We do this because even though apodization may not be requested, we still need to copy data from measured_etc into the working data objects labelled apodized anyway. Also we check whether sigma is -1. as this is the "final" apo step where no apo occurs... but ONLY if no final apo was requested, in which case its a + 1.
  {
    if (apodization_step < apodization_iterate_look_up_table.size()) // Safety check, redundant for now.
    {
      outlog << "\n*--------------------------*" << std::endl;
      outlog << "Apodization Step: " << apodization_step << std::endl;
      outlog << "*--------------------------*\n"
             << std::endl;

      current_resolution_iter = 1 / (2.4477 * sigma_apodize[apodization_step]); // Set a flag.

      outlog << "Apodization function - sigma (1/Å): " << sigma_apodize[apodization_step] << std::endl;
      outlog << "Apodization function - resolution at 1/20th of maximum height (Å): " << 1.0 / (2.4477 * sigma_apodize[apodization_step]) << std::endl;
    }
    else
    {
      outlog << "An apodization request for iterate #" << total_iterate_counter << " on apodization step #" << apodization_step << ", was invalid - this should never happen." << std::endl;
      return false;
    }

    double sigma_squared = std::pow(sigma_apodize[apodization_step], 2);

    clipper::HKL_data_base::HKL_reference_index ih;

    for (ih = hkl_target_eff_res.first(); !ih.last(); ih.next())
    {
      double s_squared = ih.invresolsq();
      double weight = exp(-s_squared / (2 * sigma_squared));

      apodization_weights[ih].fom() = weight;

      // All Data
      // apodized_f_all[ih].f()         = weight*measured_f_all[ih].f();
      // apodized_f_all[ih].sigf()      = weight*measured_f_all[ih].sigf();

      // Resolution Restricted Data
      // apodized_f[ih].f()             = weight*measured_f[ih].f();
      // apodized_f[ih].sigf()          = weight*measured_f[ih].sigf();

      // Working and Test Sets.
      // apodized_f_work_set[ih].f()    = weight * measured_f_work_set[ih].f();
      // apodized_f_work_set[ih].sigf() = weight * measured_f_work_set[ih].sigf();
      // apodized_f_test_set[ih].f()    = weight * measured_f_test_set[ih].f();
      // apodized_f_test_set[ih].sigf() = weight * measured_f_test_set[ih].sigf();
    }

    // Apodize the mean_i_wilson data, so we have a mean_I for those data that are missing from og dataset. (different hkl_object)
    // The model-based estimates for the mean transformed intensities. As intensities these are attenuated by the square of the apodization function
    // weight_sqrd = weight * weight; // pre-calc out of loop.
    // for ( ih = hkl_target_eff_res.first(); !ih.last(); ih.next() )
    //{
    //  mean_i_wilson_apodized[ih].I() = weight_sqrd * mean_i_wilson[ih].I();
    //  mean_i_wilson_apodized[ih].sigI() = weight_sqrd * mean_i_wilson[ih].sigI(); // Nothing is currently carried in SigI, but if we ever change that, this is safe behavior
    //}
  }

  else

  { // If no apodization is to be done, but we're here because its the 0th iterate (or its the final pseudo apo step i.e. sigma = -1) , then do this - which should only happen once.
    clipper::HKL_data_base::HKL_reference_index ih;
    for (ih = measured_f.first(); !ih.last(); ih.next())
    {
      apodization_weights[ih].fom() = 1.0;

      // apodized_f_all[ih].f()         = measured_f_all[ih].f();
      // apodized_f_all[ih].sigf()      = measured_f_all[ih].sigf();

      // apodized_f[ih].f()             = measured_f[ih].f();
      // apodized_f[ih].sigf()          = measured_f[ih].sigf();

      // apodized_f_work_set[ih].f()    = measured_f_work_set[ih].f();
      // apodized_f_work_set[ih].sigf() = measured_f_work_set[ih].sigf();
      // apodized_f_test_set[ih].f()    = measured_f_test_set[ih].f();
      // apodized_f_test_set[ih].sigf() = measured_f_test_set[ih].sigf();
    }

    // deal with mean wilsons...
    // for ( ih = mean_i_wilson.first(); !ih.last(); ih.next() )
    //{
    //  mean_i_wilson_apodized[i].I() = mean_i_wilson[i].I();
    //  mean_i_wilson_apodized[ih].sigI() = mean_i_wilson[ih].sigI(); // Nothing is currently carried in SigI, but if we ever change that, this is safe behavior
    //}

    current_resolution_iter = -1; // We just set this to something so stuff gets calculated in the setup.
  }

  if (apodize_on_this_step || total_iterate_counter == 0)
  {
    // Update the resolution-dependent intensity statistics for the apodized data.
    // Get rid of 0's if its a problem.

    // Clearing is not working?
    mean_I_acentric_apodized_data.clear();
    mean_I_centric_apodized_data.clear();
    /*
     outlog << " DOING THE THINGS WITH MEAN I" << std::endl;
     for (int i = 0 ; i < mean_I_acentric_apodized_data.size() ; i ++ )
     {
     outlog << "Mean_I_acentric|centric: i[" << i << "]: " << mean_I_acentric_apodized_data[i] << "|" << mean_I_centric_apodized_data[i] << std::endl;
     }*/

    mean_I_acentric_apodized_data.resize(settings.res_bins + 1, 0.0);
    mean_I_centric_apodized_data.resize(settings.res_bins + 1, 0.0);
    observation_counts_acentric_apodized_data.resize(settings.res_bins + 1, 0.0);
    observation_counts_centric_apodized_data.resize(settings.res_bins + 1, 0.0);

    // Manual Reset to 0.0, cause nothing else is working properly.
    for (int i = 0; i < mean_I_centric_apodized_data.size(); i++)
    {
      mean_I_centric_apodized_data[i] = 0.0;
      mean_I_acentric_apodized_data[i] = 0.0;
      observation_counts_acentric_apodized_data[i] = 0.0;
      observation_counts_centric_apodized_data[i] = 0.0;
    }

    if (apodize_on_this_step)
    {
      outlog << "\nIntensity Statistics for the target data with apodization applied\n"
             << std::endl;
    }
    else
    {
      outlog << "\nIntensity Statistics for the target data\n"
             << std::endl;
    }

    ipa_functions::compute_intensity_statistics(working_f,
                                                settings.res_bins,
                                                mean_I_acentric_apodized_data,
                                                mean_I_centric_apodized_data,
                                                observation_counts_acentric_apodized_data,
                                                observation_counts_centric_apodized_data,
                                                outlog, true); // Force reporting of statistics here
  }
  return true;
}

bool ipa_worker::initialise_starting_map()
{
  /*
   This function initializes the starting map, xn, and is called as the final step before starting the iterates, it must occur after the working_f and its work and test sets have been initialised. As this effectively sets the F[0,0,0] Fourier term, it is only called once. Additional calls would inflate the mean density overall.

   Relys on: work_set_f, hkl_target_eff_res, and starting_phase_set being valid.

   */
  outlog << "Effective Resolution limit"
         << " " << rlimit << std::endl;
  clipper::Resolution effective_resolution(rlimit);
  outlog << "Shannon Rate: " << settings.shannon_rate << std::endl;
  clipper::Grid_sampling grid_target(hkl_target_eff_res.spacegroup(), hkl_target_eff_res.cell(), effective_resolution, settings.shannon_rate); // create grid object. Grid spacing will be Resolution/(2*Shannon_rate). So a linear Shannon_rate of 3/2 -> Grid Spacing of Resolution/3

  // outlog << "\n TARGET STRUCTURE -- MAP INFO:\n\n";
  outlog << " Grid Dimensions: nu " << grid_target.nu() << " nv " << grid_target.nv() << " nw " << grid_target.nw() << std::endl;
  outlog << " Total Size of the Grid array: " << grid_target.size() << std::endl;

  // initialize all map objects used during the IPA with the new grid object ...
  xn.init(hkl_target_eff_res.spacegroup(), hkl_target_eff_res.cell(), grid_target); // The current iterate

  outlog << "\nGenerating the initial map (xn) from the starting phase set, and the input Structure Factor amplitudes." << std::endl;

  clipper::HKL_data<clipper::data64::F_sigF> temp_fsigf(hkl_target_eff_res);
  clipper::HKL_data<clipper::data64::F_phi> temp_fp(hkl_target_eff_res);

  if (settings.estimate_missing_amplitudes_on_first_iterate) // Replace any missing Fourier amplitudes with the expected value (based on the Wilson model).

  // Since the I's follow a Gamma Distribution, the |F|'s will follow a Nakagami distribution.
  // E(|F|) is not equal to √E(I), and corrective multipliers are needed to interconvert these expectation values
  // See notes in function substitute_Fourier_amplitudes for details

  {
    outlog << "\nThe following Fourier amplitudes are missing, and will be replaced with the expected value at the appropriate resolution" << std::endl;
    outlog << "   h   k   l     |F|expected" << std::endl;

    clipper::HKL_data_base::HKL_reference_index ih;
    for (ih = work_set_f.first(); !ih.last(); ih.next())
    {
      if (work_set_f[ih].missing() && !((ih.hkl().h() == 0) && (ih.hkl().k() == 0) && (ih.hkl().l() == 0))) // if term is missing and not F(000)
      {
        if (ih.hkl_class().centric()) // observation is centric
        {
          temp_fsigf[ih].f() = (sqrt(2.0) / sqrt(clipper::Util::pi())) * sqrt(working_mean_i_wilson[ih.hkl()].I());
          temp_fsigf[ih].sigf() = 0.0;
        }
        else // observation is acentric
        {
          temp_fsigf[ih].f() = (sqrt(clipper::Util::pi()) / 2.0) * sqrt(working_mean_i_wilson[ih.hkl()].I());
          temp_fsigf[ih].sigf() = 0.0;
        }

        outlog << std::setw(4) << ih.hkl().h() << " " << std::setw(4) << ih.hkl().k() << " " << std::setw(4) << ih.hkl().l() << " " << temp_fsigf[ih].f() << " " << std::endl;
      }
      else
      {
        temp_fsigf[ih].f() = work_set_f[ih].f();
        temp_fsigf[ih].sigf() = work_set_f[ih].sigf();
      }
    }

    temp_fp.compute(temp_fsigf, starting_phase_set, clipper::data64::Compute_fphi_from_fsigf_phifom());
  }
  else // Don't do anything about the missing Fourier terms.
  {
    temp_fp.compute(work_set_f, starting_phase_set, clipper::data64::Compute_fphi_from_fsigf_phifom());
  }

  // util::phases_to_log(temp_fp, outlog,"test", 150);

  xn.fft_from(temp_fp); // Fourier transform to generate the starting map (from the starting work_set_f)

  // Calculate the overall mean density
  expected_mean_density_overall = settings.expected_mean_density_solvent * (settings.fraction_solvent * solvent_fraction_multiplier_look_up_table[0]) + settings.expected_mean_density_protein * (1.0 - (settings.fraction_solvent * solvent_fraction_multiplier_look_up_table[0]));
  outlog << "\nxn will be modified to have the expected mean density for the crystal, which is: " << expected_mean_density_overall << std::endl;

  // Put the starting map on an approximately absolute scale
  for (clipper::Xmap_base::Map_reference_index ix = xn.first(); !ix.last(); ix.next())
  {
    xn[ix] += expected_mean_density_overall;
  }

  return true;
}

bool ipa_worker::regrid_working_maps()
{
  // If this gets called, all the target map objects needed for the IPA will be reinitialized using a new grid object, appropriate for the current resolution
  // The information in the old map objects is transferred to the new map objects via interpolation
  // The exception is the statistical weights which are regenerated by lookup.

  // This function enables gridding of the maps to change along with the apodization scheme
  // This is a pain to handle, but speeds things up enormously - at low resolution we need a much coarser grid, making everything quicker.

  outlog << "\nGridding / Regridding all map objects associated with the target structure \n"
         << std::endl;

  /*
   // Test functionality ... Export all the map objects before regridding
   if ( !(total_iterate_counter == 0) )
   {
   clipper::CCP4MAPfile map_out;
   map_out.open_write("temp1.ccp4");
   map_out.export_xmap(xn);
   map_out.close_write();
   map_out.open_write("temp2.ccp4");
   map_out.export_xmap(xn_a);
   map_out.close_write();
   map_out.open_write("temp3.ccp4");
   map_out.export_xmap(xn_b);
   map_out.close_write();
   map_out.open_write("temp4.ccp4");
   map_out.export_xmap(current_mask);
   map_out.close_write();
   map_out.open_write("temp5.ccp4");
   map_out.export_xmap(known_mask);
   map_out.close_write();
   }
   */

  clipper::Xmap_base::Map_reference_index ix;

  // temporary map objects for storing existing data
  clipper::Xmap<float> xn_old;
  clipper::Xmap<float> xn_a_old;
  clipper::Xmap<float> xn_b_old;
  clipper::Xmap<int> current_mask_old;

  if (!(total_iterate_counter == 0)) // on the very first iteration, the map objects have not yet been initialized and carry no information. So nothing to store
  {
    xn_old.init(xn.spacegroup(), xn.cell(), xn.grid_sampling());
    xn_a_old.init(xn.spacegroup(), xn.cell(), xn.grid_sampling());
    xn_b_old.init(xn.spacegroup(), xn.cell(), xn.grid_sampling());
    current_mask_old.init(xn.spacegroup(), xn.cell(), xn.grid_sampling());

    for (ix = xn.first(); !ix.last(); ix.next()) // store xn, xn_a, xn_b, current_mask
    {
      xn_old[ix] = xn[ix];
      xn_a_old[ix] = xn_a[ix];
      xn_b_old[ix] = xn_b[ix];
      current_mask_old[ix] = current_mask[ix];
    }
  }

  // Create a new grid object and use it to reinitalize the map objects (nb we are initializing the map objects for the very first time, if total_iterate_counter == 0)

  // calculate_resolution_limit();  //calculate the effective resolution limit IS NOW SET ELSEWHERE, trust rlimit.

  outlog << "Effective resolution limit: " << rlimit << std::endl;
  outlog << "Shannon Rate: " << settings.shannon_rate << std::endl;
  clipper::Resolution effective_resolution(rlimit);

  // grid_target.init(hkl_target_eff_res.spacegroup(), hkl_target_eff_res.cell(), effective_resolution, 1.5  ); // reinitialize grid object. Linear Shannon rate of 1.5 -> Grid Spacing of Resolution/3
  clipper::Grid_sampling grid_target(hkl_target_eff_res.spacegroup(), hkl_target_eff_res.cell(), effective_resolution, settings.shannon_rate); // create grid object. Linear Shannon rate of 1.5 -> Grid Spacing of Resolution/3

  // outlog << "\n TARGET STRUCTURE -- MAP INFO:\n\n";
  outlog << " Grid Dimensions: nu " << grid_target.nu() << " nv " << grid_target.nv() << " nw " << grid_target.nw() << std::endl;
  outlog << " Total Size of the Grid array: " << grid_target.size() << std::endl;

  // initialize all map objects used during the IPA with the new grid object ...
  xn.init(hkl_target_eff_res.spacegroup(), hkl_target_eff_res.cell(), grid_target);                  // The current iterate
  xn_a.init(hkl_target_eff_res.spacegroup(), hkl_target_eff_res.cell(), grid_target);                // Solution estimate A
  xn_b.init(hkl_target_eff_res.spacegroup(), hkl_target_eff_res.cell(), grid_target);                // Solution estimate B
  current_mask.init(hkl_target_eff_res.spacegroup(), hkl_target_eff_res.cell(), grid_target);        // The current mask
  known_mask.init(hkl_target_eff_res.spacegroup(), hkl_target_eff_res.cell(), grid_target);          // The known mask (test cases only)
  known_mask_inverted.init(hkl_target_eff_res.spacegroup(), hkl_target_eff_res.cell(), grid_target); // The inversion of the known mask (test cases only)
  target_weights.init(hkl_target_eff_res.spacegroup(), hkl_target_eff_res.cell(), grid_target);      // The real space statistical weights associated with the target map

  // One-off generation of the initial iterate xn
  // This should not really be placed in here, but positioning is forced by current program structure
  /*
  if (total_iterate_counter == 0)
  {

    outlog << "\n This is the very first iterate" << std::endl;
    outlog << "Generating xn from the starting phase set, and Structure Factor amplitudes" << std::endl;
    // Compute Fourier map coefficients from the Structure Factor amplitudes, plus Phases and Weights
    clipper::HKL_data<clipper::data64::F_phi> work_fp ( hkl_target_eff_res );
    work_fp.compute( apodized_f_work_set, starting_phase_set, clipper::data64::Compute_fphi_from_fsigf_phifom() );
    xn.fft_from( work_fp );  // Fourier transform to generate the starting map

    // Calculate the overall mean density of the map.
    expected_mean_density_overall = settings.expected_mean_density_solvent*(settings.fraction_solvent*solvent_fraction_multiplier_look_up_table[0]) + settings.expected_mean_density_protein*(1.0 - (settings.fraction_solvent*solvent_fraction_multiplier_look_up_table[0]));
    outlog << "\nxn will be modified to have the expected mean density, which is: " <<  expected_mean_density_overall << std::endl;
    // Put the starting map on an approximately absolute scale
    for (clipper::Xmap_base::Map_reference_index ix = xn.first(); !ix.last(); ix.next() )
    {
      xn[ix] += expected_mean_density_overall;
    }

  }
  */

  // Now restore  *all* of these map objects.

  // Check if we are fixing the current mask - if so we will restore from the input working mask, otherwise we restore from the "old" current mask
  // Strictly speaking we shouldn't need to check for input of the mask (settings.input_working_mask) but we do it here for safety
  bool fix_mask_now = ((settings.n_iterations_fix_mask > 0) && (total_iterate_counter < settings.n_iterations_fix_mask) && settings.input_working_mask);

  if (fix_mask_now)
  {
    outlog << "The mask is currently fixed and will be interpolated from the original input mask" << std::endl;
  }

  for (ix = xn.first(); !ix.last(); ix.next()) // loop over all points in the newly gridded map
  {
    clipper::Coord_frac cf = ix.coord().coord_frac(xn.grid_sampling()); // get the fractional coordinates of the current grid point

    // The iterate and the solutions. Interpolate from old grid to new grid
    //------------------------------------------------------------------------------

    if (!(total_iterate_counter == 0)) // so long as it is not the very first iteration
    {
      // linear interpolation ...
      // xn[ix] = xn_old.interp<clipper::Interp_linear>( cf );
      // xn_a[ix] = xn_a_old.interp<clipper::Interp_linear>( cf );
      // xn_b[ix] = xn_b_old.interp<clipper::Interp_linear>( cf );

      // cubic interpolation  ... use here for increased  accuracy
      // For some discussion, see Appendix B of Afonine PV, Poon BK, Read RJ, Sobolev OV, Terwilliger TC, Urzhumtsev AG, et al. Real-space refinement in PHENIX for cryo-EM and crystallography. Acta Crystallogr D Struct Biol. 2018 Jun 1;74(Pt 6):531–44.
      xn[ix] = xn_old.interp<clipper::Interp_cubic>(cf);
      xn_a[ix] = xn_a_old.interp<clipper::Interp_cubic>(cf);
      xn_b[ix] = xn_b_old.interp<clipper::Interp_cubic>(cf);
    }

    // The working mask
    // ----------------

    if (fix_mask_now) // restore from the working mask, as originally input
    {
      current_mask[ix] = lroundf(input_target_mask_working.interp<clipper::Interp_linear>(cf)); // nb input_target_mask_working is real valued, current_mask is integer valued
    }
    else // restore from current mask, if it exists
    {
      if (!(total_iterate_counter == 0)) // so long as it is not the very first iteration
      {
        current_mask[ix] = lroundf(current_mask_old.interp<clipper::Interp_linear>(cf)); // nb current_mask_old is integer valued, current_mask is integer valued
      }
    }

    // The known mask.
    // --------------

    // If we have one, interpolate from the known mask as originally input
    if (settings.input_known_mask)
    {
      known_mask[ix] = lroundf(input_target_mask_known.interp<clipper::Interp_linear>(cf)); // nb input_target_mask_known is real valued, known_mask is integer valued
    }

    // The real-space statistical weighting factors
    // --------------------------------------------
    // Lookup is expensive, so we do this once and store the result

    target_weights[ix] = 1.0 / xn.multiplicity(ix.coord());
  }

  // If we have a known mask in an achiral space group, invert and store the result
  if ((settings.input_known_mask) && (!space_group_is_chiral))
  {
    known_mask_inverted = ipa_functions::invert_mask(known_mask);
  }

  /*
   // Test functionality ... Export all the map objects after regridding
   if ( !(total_iterate_counter == 0) )
   {
   clipper::CCP4MAPfile map_out;
   map_out.open_write("temp1n.ccp4");
   map_out.export_xmap(xn);
   map_out.close_write();
   map_out.open_write("temp2n.ccp4");
   map_out.export_xmap(xn_a);
   map_out.close_write();
   map_out.open_write("temp3n.ccp4");
   map_out.export_xmap(xn_b);
   map_out.close_write();
   map_out.open_write("temp4n.ccp4");
   map_out.export_xmap(current_mask);
   map_out.close_write();
   map_out.open_write("temp5n.ccp4");
   map_out.export_xmap(known_mask);
   map_out.close_write();
   }
   */

  // All done with regridding

  return true;
}

bool ipa_worker::calculate_reference_histograms()
{
  /*
   If we are apodizing the reference data together with the target data, we need to recalculate the reference histograms every time the apodization scheme changes
   If we are not apodizing the reference data, we just need to calculate the reference histograms from the unmodified map on the very first apodization step, and then leave it alone
   */

  // Generate the reference map used to calculate the electron density and gradient histograms

  // First get the weights required for the apodization, if any
  clipper::HKL_data_base::HKL_reference_index ih;
  for (ih = measured_fp_reference.first(); !ih.last(); ih.next())
  {
    if ((settings.apodize_reference_data) && (apodize_on_this_step))
    {
      double s_squared = ih.invresolsq();
      double sigma_squared = std::pow(sigma_apodize[apodization_step], 2);
      double weight = exp(-s_squared / (2 * sigma_squared));
      apodized_fp_reference[ih].f() = weight * measured_fp_reference[ih].f();
    }
    else
    {
      apodized_fp_reference[ih].f() = measured_fp_reference[ih].f();
    }
  }
  // Now calculate the map
  reference_map.fft_from(apodized_fp_reference);

  // Generate electron density histogram, if required
  if (settings.histogram_matching)
  {
    // Calculate the normalized 1D density histogram

    outlog << std::endl;
    outlog << "Calculating electron density histogram for the reference structure" << std::endl;

    bool calculate_higher_moments = true;
    ipa_functions::calculate_summary_statistics(
        reference_map,
        reference_mask,
        reference_weights,
        protein_region_flag,
        calculate_higher_moments,
        protein_density_stats_reference.mean,
        protein_density_stats_reference.variance,
        protein_density_stats_reference.skewness,
        protein_density_stats_reference.kurtosis,
        protein_density_stats_reference.minimum,
        protein_density_stats_reference.maximum);

    density_histogram_reference_raw.minimum = protein_density_stats_reference.minimum;
    density_histogram_reference_raw.maximum = protein_density_stats_reference.maximum;

    /*
     outlog << " TESTING INPUTS FOR HISTOGRAM OUTPUT:" << std::endl;
     clipper::Xmap_base::Map_reference_index ix;
     int n_counter = 0;
     for ( ix = reference_map.first(); !ix.last(); ix.next() )
     {
     outlog << reference_map[ix] << " " << reference_mask[ix] << " " << reference_weights[ix] << " " << std::endl;

     n_counter++;
     if (n_counter > 100) {break;}
     }
     */

    ipa_functions::calculate_histogram(reference_map, reference_mask, reference_weights, protein_region_flag, density_histogram_reference_raw, outlog);

    // Smooth the raw histogram, using locally weighted regression. Could be tidied up, if there's some spare time.

    // The smoothed histogram has the same bin structure as the raw histogram
    density_histogram_reference.minimum = density_histogram_reference_raw.minimum;
    density_histogram_reference.maximum = density_histogram_reference_raw.maximum;

    outlog << std::endl;
    outlog << "Smoothing the Raw histogram using local regression" << std::endl;
    outlog << std::endl;

    if (settings.verbose)
    {
      outlog << "   Bin #  Lower_Boundary   Raw_Histogram  Smoothed_Histogram Residuals" << std::endl;
    }

    float *loess_bin_boundaries;                                                       // Declare a pointer to an array of float data
    loess_bin_boundaries = new float[density_histogram_reference.number_of_intervals]; // Use the new operator to dynamically allocate space for the array
    float *loess_residuals;                                                            // Declare a pointer to an array of float data
    loess_residuals = new float[density_histogram_reference.number_of_intervals];      // Use the new operator to dynamically allocate space for the array
    float *loess_rw;                                                                   // Declare a pointer to an array of float data
    loess_rw = new float[density_histogram_reference.number_of_intervals];             // Use the new operator to dynamically allocate space for the array

    float histogram_scale_factor = (((float)(density_histogram_reference.number_of_intervals)) - 1.0) / (density_histogram_reference.maximum - density_histogram_reference.minimum);

    for (int j = 0; j < density_histogram_reference.number_of_intervals; j++)
    {
      loess_bin_boundaries[j] = density_histogram_reference.minimum + (float)(j) / histogram_scale_factor;
    }

    // TODO: adjust the fraction of points based on both the number of intervals in the histogram , and the number of points in the map.
    // "Typically" 0.050 works well when there are about 250 bins
    // Something to investigate and fix on a rainy day.

    int nsteps = 2;
    float delta = 0.0;
    float fraction_of_points = 0.050;

    lowess(&loess_bin_boundaries[0], &density_histogram_reference_raw.interval[0], &density_histogram_reference.number_of_intervals, &fraction_of_points, &nsteps, &delta, &density_histogram_reference.interval[0], &loess_rw[0], &loess_residuals[0]);

    float new_sum = 0;

    for (int j = 0; j < density_histogram_reference.number_of_intervals; j++)
    {
      new_sum += density_histogram_reference.interval[j];
      if (settings.verbose)
      {
        outlog << std::fixed << std::setw(5) << j + 1 << " " << std::setw(15) << loess_bin_boundaries[j] << " " << std::setw(15) << density_histogram_reference_raw.interval[j] << " " << std::setw(15) << density_histogram_reference.interval[j] << " " << std::setw(15) << loess_residuals[j] << std::endl;
      }
    }

    // Renormalize the histogram after smoothing
    for (int j = 0; j < density_histogram_reference.number_of_intervals; j++)
    {
      density_histogram_reference.interval[j] /= new_sum;
    }

    delete[] loess_bin_boundaries;
    delete[] loess_residuals;
    delete[] loess_rw;
  }

  // Generate gradient magnitude histogram, if required
  if (settings.gradient_matching)
  {

    outlog << std::endl;
    outlog << "Calculating the gradient magnitude histogram for the reference structure" << std::endl;

    clipper::Xmap<float> reference_grad_x(hkl_reference.spacegroup(), hkl_reference.cell(), grid_reference);    // create map object to store gradient component along X
    clipper::Xmap<float> reference_grad_y(hkl_reference.spacegroup(), hkl_reference.cell(), grid_reference);    // create map object to store gradient component along Y
    clipper::Xmap<float> reference_grad_z(hkl_reference.spacegroup(), hkl_reference.cell(), grid_reference);    // create map object to store gradient component along Z
    clipper::Xmap<float> reference_grad_magn(hkl_reference.spacegroup(), hkl_reference.cell(), grid_reference); // create map object to store gradient magnitude

    ipa_functions::calculate_gradient(apodized_fp_reference, reference_grad_x, reference_grad_y, reference_grad_z, reference_grad_magn, outlog);

    bool calculate_higher_moments = true;
    ipa_functions::calculate_summary_statistics(reference_grad_magn,
                                                reference_mask,
                                                reference_weights,
                                                protein_region_flag,
                                                calculate_higher_moments,
                                                protein_gradient_stats_reference.mean,
                                                protein_gradient_stats_reference.variance,
                                                protein_gradient_stats_reference.skewness,
                                                protein_gradient_stats_reference.kurtosis,
                                                protein_gradient_stats_reference.minimum,
                                                protein_gradient_stats_reference.maximum);

    gradient_histogram_reference.minimum = protein_gradient_stats_reference.minimum;
    gradient_histogram_reference.maximum = protein_gradient_stats_reference.maximum;

    ipa_functions::calculate_histogram(reference_grad_magn, reference_mask, reference_weights, protein_region_flag, gradient_histogram_reference, outlog);

    // Don't smooth the gradient magnitude histogram right now,  as there may be problems with it becoming negative at zero
    /*
     gradient_histogram_reference_raw.minimum = protein_gradient_stats_reference.minimum;
     gradient_histogram_reference_raw.maximum = protein_gradient_stats_reference.maximum;

     ipa_functions::calculate_histogram(reference_grad_magn, reference_mask, reference_weights, protein_region_flag, gradient_histogram_reference_raw);

     // Smooth the raw histogram, using locally weighted regression. Could be tidied up, if there's some spare time.
     // The smoothed histogram has the same bin structure as the raw histogram

     gradient_histogram_reference.minimum  = gradient_histogram_reference_raw.minimum;
     gradient_histogram_reference.maximum  = gradient_histogram_reference_raw.maximum;

     outlog <<  std::endl;
     outlog << "Smoothing the Raw gradient magnitude histogram using local regression" << std::endl;
     outlog <<  std::endl;
     outlog << "   Bin #  Lower_Boundary   Raw_Histogram  Smoothed_Histogram Residuals" << std::endl;

     float *loess_bin_boundaries; // Declare a pointer to an array of float data
     loess_bin_boundaries = new float[gradient_histogram_reference.number_of_intervals]; // Use the new operator to dynamically allocate space for the array
     float *loess_residuals; // Declare a pointer to an array of float data
     loess_residuals = new float[gradient_histogram_reference.number_of_intervals]; // Use the new operator to dynamically allocate space for the array
     float *loess_rw; // Declare a pointer to an array of float data
     loess_rw= new float[gradient_histogram_reference.number_of_intervals]; // Use the new operator to dynamically allocate space for the array


     float histogram_scale_factor = (((float)(gradient_histogram_reference.number_of_intervals))-1.0)/(gradient_histogram_reference.maximum - gradient_histogram_reference.minimum);
     for ( int j= 0; j < gradient_histogram_reference.number_of_intervals; j++ )
     {
     loess_bin_boundaries[j] = gradient_histogram_reference.minimum + (float)(j)/histogram_scale_factor;
     }
     int nsteps = 2;
     float delta = 0.0;
     float fraction_of_points = 0.050;

     lowess(&loess_bin_boundaries[0], &gradient_histogram_reference_raw.interval[0], &gradient_histogram_reference.number_of_intervals, &fraction_of_points, &nsteps, &delta, &gradient_histogram_reference.interval[0], &loess_rw[0], &loess_residuals[0]);

     float new_sum = 0;
     for ( int j= 0; j < gradient_histogram_reference.number_of_intervals; j++ )
     {
     new_sum += gradient_histogram_reference.interval[j];
     outlog << std::fixed <<
     std::setw(5) << j+1 << " " <<
     std::setw(15) << loess_bin_boundaries[j]  << " " <<
     std::setw(15) << gradient_histogram_reference_raw.interval[j] << " " <<
     std::setw(15) << gradient_histogram_reference.interval[j] << " " <<
     std::setw(15) << loess_residuals[j] << std::endl;
     }

     // Renormalize the histogram after smoothing

     for ( int j= 0; j < gradient_histogram_reference.number_of_intervals; j++ )
     {
     gradient_histogram_reference.interval[j] /= new_sum;
     }

     delete [] loess_bin_boundaries;
     delete [] loess_residuals;
     delete [] loess_rw;
     delete [] gradient_histogram_reference_raw.interval;
     */
  }

  return true;
}

void ipa_worker::check_if_histogram_or_gradient_matching()
{
  /*
   Resolves which operation we are doing on the current iterate - currently based on the ap_step_iterate_counter (the number of iterates that have occured SINCE the last change in apodization )
   TODO: - switch so that the decision is made based on the total iterate counter - done. Now this needs to be called every iterate.
   */
  histogram_matching_on_current_iterate = false;
  gradient_matching_on_current_iterate = false;

  if (settings.histogram_matching && (((total_iterate_counter + 1) % settings.histogram_matching_interval) == 0))
  {
    histogram_matching_on_current_iterate = true;
  }

  if (settings.gradient_matching && (((total_iterate_counter + 1) % settings.gradient_matching_interval) == 0))
  {
    gradient_matching_on_current_iterate = true;
    histogram_matching_on_current_iterate = false;
  }
}

bool ipa_worker::generate_a_mask(const double current_filter_radius)
{
  if (!((settings.n_iterations_fix_mask > 0) && (total_iterate_counter < settings.n_iterations_fix_mask))) // The mask is being generated from the current density estimate, though possibly not on each iteration
  {
    if ((total_iterate_counter == 0) || ((total_iterate_counter % settings.compute_mask_interval) == 0))
    // Mask generation is scheduled on this iteration (or it is the very first iteration)
    {
      outlog << "\nBeginning Mask generation" << std::endl;
      // Calculate the current fraction solvent using the appropriate multiplier
      float fraction_solvent_current = settings.fraction_solvent * solvent_fraction_multiplier_look_up_table[total_iterate_counter];
      // Record the current fraction solvent used.
      stat_xn_fraction_solvent = fraction_solvent_current;

      outlog << "Base value for the solvent fraction: " << settings.fraction_solvent << "." << std::endl;
      outlog << "Current value for the solvent fraction multiplier: " << solvent_fraction_multiplier_look_up_table[total_iterate_counter] << "." << std::endl;
      outlog << "Current value for the solvent fraction: " << fraction_solvent_current << "." << std::endl;
      outlog << "Filter radius for mask generation: " << current_filter_radius << "." << std::endl;
      if (total_iterate_counter == 0) // very first iteration only, use xn for the calculation.
      {
        outlog << "Generating a mask from xn (as this is the first iteration)" << std::endl;
        ipa_functions::calculate_mask(settings.compute_mask_from_local_mean, settings.compute_mask_from_local_variance, settings.compute_mask_from_local_mean_and_variance,
                                      current_filter_radius, settings.mask_weighting_factor_mean,
                                      xn, fraction_solvent_current, current_mask, outlog);
      }
      else // ... otherwise use the solution estimate xn_b for the calculation
      {
        outlog << "Generating a mask from solution estimate xn_b." << std::endl;
        ipa_functions::calculate_mask(settings.compute_mask_from_local_mean, settings.compute_mask_from_local_variance, settings.compute_mask_from_local_mean_and_variance,
                                      current_filter_radius, settings.mask_weighting_factor_mean,
                                      xn_b, fraction_solvent_current, current_mask, outlog);
      }

      //  connectivity reinforcement on the target mask, if requested. We erase small "islands" and fill small "voids"

      bool write_on_current_cycle = false;
      bool erase_on_current_cycle = false;
      bool fill_on_current_cycle = false;

      // TODO: ... The output_mask setting is currently ignored. If we reenable it, should probably tie it to the total iterate counter rather than ap_step_iterate counter
      if ((settings.output_mask && (((ap_step_iterate_counter + 1) % settings.mask_output_interval) == 0)) || (ap_step_iterate_counter == n_iterations_ap_step - 1))
        write_on_current_cycle = true;

      // TODO: ... These mask editing settings should probably be tied to the total iterate counter rather than ap_step_iterate counter
      if (settings.erase_islands && (((ap_step_iterate_counter + 1) % settings.erase_islands_interval) == 0))
        erase_on_current_cycle = true;
      if (settings.fill_voids && (((ap_step_iterate_counter + 1) % settings.fill_voids_interval) == 0))
        fill_on_current_cycle = true;

      if (rlimit <= settings.envelope_editing_limit)
      {
        erase_on_current_cycle = false;
        fill_on_current_cycle = false;
      }

      // If mask output is requested, and the mask is to be edited, output it before any editing is performed
      // Currently disabled
      /*
       if (write_on_current_cycle  && (erase_on_current_cycle || fill_on_current_cycle ) )
       {
       clipper::CCP4MAPfile map_out;
       snprintf(io_map_filename, sizeof(io_map_filename), "%smask_pre_edit_run%d_step%d_iteration%d.ccp4", settings.working_dir.c_str(), my_run_id, apodization_step, total_iterate_counter);
       map_out.open_write( io_map_filename );
       map_out.export_xmap(current_mask); // write mask
       map_out.close_write();
       }*/

      int n_connected_sets_eliminated_protein;
      int n_connected_sets_retained_protein;
      int n_connected_sets_eliminated_solvent;
      int n_connected_sets_retained_solvent;

      if (erase_on_current_cycle || fill_on_current_cycle)
      {
        outlog << "Erasing/Filling the mask" << std::endl;
        ipa_functions::enforce_mask_connectivity(erase_on_current_cycle, fill_on_current_cycle, settings.erase_islands_threshold, settings.fill_voids_threshold, current_mask,
                                                 n_connected_sets_eliminated_protein, n_connected_sets_retained_protein, n_connected_sets_eliminated_solvent, n_connected_sets_retained_solvent, outlog);
      }
      /*

       // If mask output is requested, output the final mask
       // Currently disabled

       if (write_on_current_cycle)
       {
       clipper::CCP4MAPfile map_out;
       snprintf(io_map_filename, sizeof(io_map_filename), "%smask_run%d_step%d_iteration%d.ccp4",settings.working_dir.c_str(), my_run_id, apodization_step, ap_step_iterate_counter);
       map_out.open_write( io_map_filename );
       map_out.export_xmap(current_mask); // write mask
       map_out.close_write();
       }
       */
      outlog << "Ending Mask generation" << std::endl;
    }
    else
    {
      outlog << "\nMask will not be updated on this iteration\n"
             << std::endl;
    }
  }
  else // The mask is fixed to its input state
  {
    outlog << "\nThe Mask is currently fixed to its input state, and not being iteratively updated\n"
           << std::endl;
  }
  return true;
}

// ------------------------------------------------------- //
// --------------------- OUTPUTS ------------------------- //
// ------------------------------------------------------- //

void ipa_worker::summarise_iterate_statistics(std::ostream &outstream)
{
  /*
   This function is called on any iterate where we want summary statistics to be sent to the log file for the associated run.

   The calculation of summary statistics costs time, and it's generally overkill to compute these statistics on every iteration.

   member variables xn_a and xn_b contain the solution estimates from the IPA.
   member variables xn_b_fp and xn_a_fp contain the correspondent Fourier amplitudes and phases

   xn_a and xn_b should have been calculated prior to entering the function.
   xn_b_fp and xn_a_fp are not guaranteed to have been calculated for every IPA ... we therefore check the relevant flags  xn_a_fp_is_valid & xn_b_fp_is_valid to see if they need to be generated.

   The function computes:

   A measure of algorithm convergence (aka the agreement between the two solution estimates)
   Agreement with the real space constraints. This is evaluated using the solution estimate PB(…) (=xn_b)  where the Fourier space projection is performed last (& calculation of statistics will use xn_b directly)
   Agreement with the Fourier space constraints. This is evaluated using the solution estimate PA(…) (=xn_a)  where the Real Space projection is performed last (& calculation of statistics will use xn_a_fp )
   Phase agreement statistics, if an known phase set has been input.
   Mask agreement statistics, if an known mask has been input.

   Once the relevant statistics are evaluated a line is printed to the output log for this run (the ostream object)

   */

  // xn_b_fp and xn_a_fp may not be valid on entry, depending on the IPA employed, so generate them if needed

  if (!xn_a_fp_is_valid && xn_a_is_valid)
  {
    xn_a.fft_to(xn_a_fp);
    xn_a_fp_is_valid = true;
    outlog << "Generated the FFT of Solution A for statistical use" << std::endl;
  }
  else if (!xn_a_fp_is_valid)
  {
    outlog << "Error when calculating statistics, couldn't generate xn_a_fp from xn_a. xn_a invalid." << std::endl;
  }

  if (!xn_b_fp_is_valid && xn_b_is_valid)
  {
    xn_b.fft_to(xn_b_fp);
    xn_b_fp_is_valid = true;
    outlog << "Generated the FFT of Solution B for statistical use" << std::endl;
  }
  else if (!xn_b_fp_is_valid)
  {
    outlog << "Error when calculating statistics, couldn't generate xn_b_fp from xn_b. xn_b invalid." << std::endl;
  }

  // xn_a and xn_b *should* be valid on entry (this function should only be called at the conclusion of an IPA iterate) but we check the associated flags for safety

  if (!xn_a_is_valid && xn_a_fp_is_valid)
  {
    xn_a.fft_from(xn_a_fp);
    xn_a_is_valid = true;
  }
  else if (!xn_a_is_valid)
  {
    outlog << "Error when calculating statistics, couldn't generate xn_a from xn_a_fp. xn_a_fp invalid." << std::endl;
  }

  if (!xn_b_is_valid && xn_b_fp_is_valid)
  {
    xn_b.fft_from(xn_b_fp);
    xn_b_is_valid = true;
  }
  else if (!xn_b_is_valid)
  {
    outlog << "Error when calculating statistics, couldn't generate xn_b from xn_b_fp. xn_b_fp invalid." << std::endl;
  }

  // Now the two solution estimates and their Fourier transforms are valid for sure,  so we are good to go.

  // Calculate the RMSD between the two solution estimates, which is a measure of IPA convergence.

  stat_delta_n = ipa_functions::rmsd_between_maps(xn_a, xn_b, target_weights);
  outlog << std::endl
         << "RMSD between the two solution estimates (Delta_n): " << stat_delta_n << std::endl;

  // Calculate the real space agreement statistics (with xn_b)
  calculate_real_space_agreement_statistics();

  // Calculate the Fourier space agreement statistics (with xn_a_fp)
  calculate_fourier_space_agreement_statistics();

  // Calculate phase agreement statistics if we have known phases (test cases only). Again use solution A (xn_a_fp) for this
  // Note that the map correlation coefficient is now calculated here, in addition to the mean absolute phase difference
  if (settings.input_known_phases)
  {
    calculate_agreement_with_known_phase_set();
  }

  // Calculate some mask agreement statistics if we have a known mask (test cases only)

  if (settings.input_known_mask)
  {
    calculate_agreement_with_known_mask();
  }

  // Quickly cheat and check if a solution was found!
  if (!solution_found && settings.input_known_phases)
  {
    if (stat_xn_map_cc > 0.5 || stat_xn_map_cc_inverse > 0.5) // Check both solutions
    {
      solution_found = true;                          // Set to true, a solution can not be found twice.
      solution_found_iterate = total_iterate_counter; // Set the iterate the solution was found.
    }
  }

  write_formatted_summary_to_log(outstream);

  if (settings.verbose)
  {
    outlog << "Data written to file..." << std::endl;
  }

  if (settings.output_phases)
  // Output phase to an mtz if requested.
  // TODO ... This should be moved to its own function if we want to be able to write phases independent of writing summary statistics, which is the way the settings were originally constructed
  {
    char io_mtz_filename[256];
    clipper::CCP4MTZfile mtzout;

    //           clipper::String output_col_fp         = "/*/*/PB_PA_xn";
    clipper::String output_col_fp = "/*/*/[F,PHIC]";
    snprintf(io_mtz_filename, sizeof(io_mtz_filename), "%sxn_run%d_step%d_iteration%d.mtz", settings.working_dir.c_str(), my_run_id, apodization_step, ap_step_iterate_counter);
    mtzout.set_column_label_mode(clipper::CCP4MTZfile::Legacy);
    mtzout.open_write(io_mtz_filename);
    mtzout.export_hkl_info(hkl_target_eff_res);

    if (!xn_b_fp_is_valid && xn_b_is_valid)
    {
      outlog << "Calculating xn_b_fp..." << std::endl;
      xn_b.fft_to(xn_b_fp); // Quickly calculate this if we havent already...
      xn_b_fp_is_valid = true;
    }

    mtzout.export_hkl_data(xn_b_fp, output_col_fp); // Pretty sure xn_b_fp is also appropriate here, as even for Error Reduction, xn() is just copied from xn_b, and xn_b will always be the most recent fourier projection.
    mtzout.close_write();
  }
}

void ipa_worker::write_formatted_summary_to_log(std::ostream &outstream)
{
  // Write Header to filestream on first iterate.
  if (total_iterate_counter == 0)
  {
    outstream
        << "       Run#                    " // 1
        << "       Apodization_Step        "
        << "       Iterate                 "
        << "  Phase    Real    Cusum       "
        << "      xn_solvent_mean          "
        << "      xn_solvent_variance      "
        << "      xn_protein_mean          "
        << "      xn_protein_variance      "
        << "      xn_histogram_agreement   "
        << "      xn_correlation_coeff     "
        << "      xn_r-factor              "
        << "      xn_correlation_coeff_free"
        << "      xn_r-factor_free         "
        << "       delta_n                 "
        << "       beta                    "
        << "       radius                  "
        << "       fraction_solvent        ";
    if (settings.input_known_phases)
    {
      outstream
          << "  xn_mean_phase_diff_all       "
          << "  xn_mean_phase_diff_cent      "
          << "  xn_mean_phase_diff_acent     "
          << "  xn_map_cc                    " // 23
          << "  xn_mean_phase_diff_all_inv   " // 24
          << "  xn_mean_phase_diff_cent_inv  " // 25
          << "  xn_mean_phase_diff_acent_inv " // 26
          << "  xn_map_cc_inv                "
          << "  solution_iter                ";
    }
    if (settings.input_known_mask)
    {
      outstream
          << "       mask_cc                 "
          << "       smc                     "
          << "       Cohen_kappa             "
          << "       Matthews_cc             "
          << "       mask_cc_inverse         "
          << "       smc_inverse             "
          << "       Cohen_kappa_inverse     "
          << "       Matthews_cc_inverse     ";
    }
    outstream << std::endl;
  }

  // Write Data to stream
  char buffer[512];
  std::sprintf(buffer, "%10i                       %10i                       %10i             %10.1f %10.3f %10.1f                 ",
               my_run_id,
               apodization_step,
               total_iterate_counter,
               stat_phase_consistency,
               stat_real_space_consistency,
               stat_phase_cusum);
  outstream << buffer;

  std::sprintf(buffer, "%10.8f                     %10.8f                     %10.8f                     %10.8f                     ",
               stat_xn_solvent_mean,
               stat_xn_solvent_variance,
               stat_xn_protein_mean,
               stat_xn_protein_variance);
  outstream << buffer;

  std::sprintf(buffer, "%10.8f                     ", stat_xn_histogram_agreement);
  outstream << buffer;

  std::sprintf(buffer, "%10.8f                     %10.6f                     ",
               stat_xn_correlation_coeff_work,
               stat_xn_r_factor_work);
  outstream << buffer;

  std::sprintf(buffer, "%10.8f                     %10.6f                     ",
               stat_xn_correlation_coeff_free,
               stat_xn_r_factor_free);
  outstream << buffer;

  std::sprintf(buffer, "%10.8f                     ", stat_delta_n);
  outstream << buffer;
  std::sprintf(buffer, "%10.8f                     ", beta);
  outstream << buffer;
  std::sprintf(buffer, "%10.8f                     ", stat_xn_radius);
  outstream << buffer;
  std::sprintf(buffer, "%10.8f                     ", stat_xn_fraction_solvent);
  outstream << buffer;
  if (settings.input_known_phases)
  {
    std::sprintf(buffer, "%10.4f                     %10.4f                     %10.4f                     %10.4f                     ",
                 stat_xn_mean_phase_diff_all,
                 stat_xn_mean_phase_diff_centric,
                 stat_xn_mean_phase_diff_acentric,
                 stat_xn_map_cc);
    outstream << buffer;

    std::sprintf(buffer, "%10.4f                     %10.4f                     %10.4f                     %10.4f                     %10i                     ",
                 stat_xn_mean_phase_diff_all_inverse,
                 stat_xn_mean_phase_diff_centric_inverse,
                 stat_xn_mean_phase_diff_acentric_inverse,
                 stat_xn_map_cc_inverse,
                 solution_found_iterate);
    outstream << buffer;
  }
  if (settings.input_known_mask)
  {
    std::sprintf(buffer, "%10.7f                     %10.7f                     %10.7f                     %10.7f                     ",
                 stat_mask_cc,
                 stat_smc,
                 stat_cohen_kappa,
                 stat_mcc);
    outstream << buffer;

    std::sprintf(buffer, "%10.7f                     %10.7f                     %10.7f                     %10.7f                     ",
                 stat_mask_cc_inverse,
                 stat_smc_inverse,
                 stat_cohen_kappa_inverse,
                 stat_mcc_inverse);
    outstream << buffer;
  }
  outstream << std::endl;
}

void ipa_worker::calculate_real_space_agreement_statistics()
{
  outlog << "\nCalculating agreement between Solution Estimate B, and the real space constraints" << std::endl;

  // Validity of xn_b has been checked in calling routine, so no need to do that again here

  map_region_summary_statistics solvent_density_stats;
  map_region_summary_statistics protein_density_stats;

  bool calculate_higher_moments = true;
  // Compute map statistics for the solvent region
  ipa_functions::calculate_summary_statistics(xn_b,
                                              current_mask,
                                              target_weights,
                                              solvent_region_flag,
                                              calculate_higher_moments,
                                              solvent_density_stats.mean,
                                              solvent_density_stats.variance,
                                              solvent_density_stats.skewness,
                                              solvent_density_stats.kurtosis,
                                              solvent_density_stats.minimum,
                                              solvent_density_stats.maximum);
  stat_xn_solvent_mean = (float)solvent_density_stats.mean;
  stat_xn_solvent_variance = (float)solvent_density_stats.variance;
  stat_xn_solvent_minimum = (float)solvent_density_stats.minimum;
  stat_xn_solvent_maximum = (float)solvent_density_stats.maximum;

  outlog << "Solution Estimate B - real space representation: Solvent region variance: " << solvent_density_stats.variance << std::endl;

  // Compute map statistics for the protein region
  ipa_functions::calculate_summary_statistics(xn_b,
                                              current_mask,
                                              target_weights,
                                              protein_region_flag,
                                              calculate_higher_moments,
                                              protein_density_stats.mean,
                                              protein_density_stats.variance,
                                              protein_density_stats.skewness,
                                              protein_density_stats.kurtosis,
                                              protein_density_stats.minimum,
                                              protein_density_stats.maximum);
  stat_xn_protein_mean = (float)protein_density_stats.mean;
  stat_xn_protein_variance = (float)protein_density_stats.variance;
  stat_xn_protein_minimum = (float)protein_density_stats.minimum;
  stat_xn_protein_maximum = (float)protein_density_stats.maximum;

  // Don't currently compute or store overall statistics. Should we ?
  // double mean;
  // double variance;
  // double skewness;
  // double kurtosis;
  // double minimum;
  // double maximum;

  // ipa_functions::calculate_summary_statistics(xn_b, current_mask, target_weights, all_regions_flag, calculate_higher_moments, mean, variance, skewness, kurtosis, minimum, maximum);
  /*
   outlog << "Solution Estimate B - real space representation: Overall Mean: " << mean << std::endl;
   outlog << "Solution Estimate B - real space representation: Overall Variance: " << variance << std::endl;

   outlog << "Solution Estimate B - real space representation: Mean in the solvent region: "     << xn_solvent_mean[ap_step_iterate_counter]     << std::endl;
   outlog << "Solution Estimate B - real space representation: Variance in the solvent region: " << xn_solvent_variance[ap_step_iterate_counter] << std::endl;
   outlog << "Solution Estimate B - real space representation: Minimum in the solvent region: "  << xn_solvent_minimum[ap_step_iterate_counter]  << std::endl;
   outlog << "Solution Estimate B - real space representation: Maximum in the solvent region: "  << xn_solvent_maximum[ap_step_iterate_counter]  << std::endl;

   outlog << "Solution Estimate B - real space representation: Mean in the protein region: "     << xn_protein_mean[ap_step_iterate_counter]     << std::endl;
   outlog << "Solution Estimate B - real space representation: Variance in the protein region: " << xn_protein_variance[ap_step_iterate_counter] << std::endl;
   outlog << "Solution Estimate B - real space representation: Minimum in the protein region: "  << xn_protein_minimum[ap_step_iterate_counter]  << std::endl;
   outlog << "Solution Estimate B - real space representation: Maximum in the protein region: "  << xn_protein_maximum[ap_step_iterate_counter]  << std::endl;
   */

  if (settings.histogram_matching) // Calculate histogram distance metric if we are histogram matching
  {
    if (settings.verbose)
    {
      outlog << "Histogram for solution estimate B:" << std::endl;
    }
    density_histogram.minimum = protein_density_stats.minimum;
    density_histogram.maximum = protein_density_stats.maximum;

    ipa_functions::calculate_histogram(xn_b, current_mask, target_weights, protein_region_flag, density_histogram, outlog, settings.verbose);

    if (settings.verbose)
    {
      outlog << "Will rescale the reference density so it has the mean and variance of the solution estimate " << std::endl;
    }

    double a;
    double b;

    // Suppose we have reference density with mean E(X) and Variance Var(X)
    // We need a linear transformation of the reference density ( Y=aX+b ), such that we end up with the required mean E(Y) and variance Var(Y).
    // Under this linear transformation E(Y) = aE(X)+ b and Var(Y) = a^2 * Var(X)
    // Hence the parameters of the linear transfornation are  a = Sqrt( Var(Y)/Var(X) ), and b = E(Y)- aE(X)

    // In  this instance
    // E(X) = mean of the reference density in the protein region
    // Var(X) = variance of the reference density in the protein region
    // E(Y) = mean of xn_b in the protein region
    // Var(Y) = variance of xn_b in the protein region

    a = sqrt(protein_density_stats.variance / protein_density_stats_reference.variance); //  a = Sqrt( Var(Y)/Var(X) )
    b = protein_density_stats.mean - (a * protein_density_stats_reference.mean);         // b = E(Y)- aE(X)

    if (settings.verbose)
    {
      outlog << std::endl;
      outlog << "Linear rescaling factors for the reference density" << std::endl;
      outlog << " slope (a):      " << a << std::endl;
      outlog << " intercept (b):  " << b << std::endl;
    }

    clipper::Xmap_base::Map_reference_index ix;
    for (ix = reference_map_rescaled.first(); !ix.last(); ix.next())
    {
      reference_map_rescaled[ix] = a * reference_map[ix] + b;
    }

    // calculate a reference histogram for the rescaled map, using the *same* bin intervals as for solution estimate B
    if (settings.verbose)
    {
      outlog << "Histogram of the rescaled Reference density, using the *same* bin intervals:" << std::endl;
    }
    density_histogram_reference_rescaled.minimum = density_histogram.minimum;
    density_histogram_reference_rescaled.maximum = density_histogram.maximum;
    ipa_functions::calculate_histogram(reference_map_rescaled, reference_mask, reference_weights, protein_region_flag, density_histogram_reference_rescaled, outlog, settings.verbose);

    stat_xn_histogram_agreement = ipa_functions::earth_movers_distance(density_histogram, density_histogram_reference_rescaled);

    outlog << "Solution Estimate B - real space representation: Wasserstein distance from reference histogram: " << stat_xn_histogram_agreement << std::endl;
  }
}

void ipa_worker::calculate_fourier_space_agreement_statistics()
{
  outlog << "\nCalculating agreement between Solution Estimate A, and the Fourier space constraints" << std::endl;

  // Validity of xn_a_fp has been checked in calling routine, so no need to to that again here

  // Calculate the correlation coefficient and the R-factor between xn_a_fp and the experimental amplitudes (appropriately apodized of course)

  double esd_threshold;
  double scale_factor;
  double cc;
  double rf;
  double cc_free;
  double rf_free;

  ipa_functions::simple_linear_scale(work_set_f, test_set_f, apodization_weights, xn_a_fp, settings.use_test_set, esd_threshold, settings.res_bins, scale_factor, cc, rf, cc_free, rf_free, outlog, settings.verbose); // Compute linear scale Factor that will put Fourier amplitudes on the same scale as the experimental measurements, Calculate correlation coefficient and R-factor

  stat_xn_correlation_coeff_work = (float)cc;
  stat_xn_r_factor_work = (float)rf;

  outlog << "Solution estimate A - Fourier space representation : Correlation Coefficient with work set: " << stat_xn_correlation_coeff_work << std::endl;
  outlog << "Solution estimate A - Fourier space representation : R-Factor with work set: " << stat_xn_r_factor_work << std::endl;

  if (settings.use_test_set)
  {
    stat_xn_correlation_coeff_free = (float)cc_free;
    stat_xn_r_factor_free = (float)rf_free;
    outlog << "Solution estimate A - Fourier space representation : Correlation Coefficient with free set: " << stat_xn_correlation_coeff_free << std::endl;
    outlog << "Solution estimate A - Fourier space representation : R-Factor with free set: " << stat_xn_r_factor_free << std::endl;
  }
  else
  {
    stat_xn_correlation_coeff_free = 0.0;
    stat_xn_r_factor_free = 0.0;
  }
}

void ipa_worker::calculate_agreement_with_known_mask()
{
  settings.compute_real_space_agreement_statistics = true;
  clipper::Coord_frac origin_shift;
  outlog << "\nCalculating Agreement between Current Working Mask and known mask  \n"
         << std::endl;

  ipa_functions::compute_mask_agreement(current_mask, known_mask, origin_shift, settings.compute_real_space_agreement_statistics, stat_mask_cc, stat_smc, stat_cohen_kappa, stat_mcc, outlog, settings.verbose);

  if (!space_group_is_chiral)
  {
    outlog << "\nCalculating Agreement between Current Working Mask and Inversion of the known mask  \n"
           << std::endl;

    ipa_functions::compute_mask_agreement(current_mask, known_mask_inverted, origin_shift, settings.compute_real_space_agreement_statistics, stat_mask_cc_inverse, stat_smc_inverse, stat_cohen_kappa_inverse, stat_mcc_inverse, outlog, settings.verbose);
  }
  else
  { // Just set these to 0 if they are not being calculated...
    stat_mask_cc_inverse = 0.0;
    stat_smc_inverse = 0.0;
    stat_cohen_kappa_inverse = 0.0;
    stat_mcc_inverse = 0.0;
  }
}

void ipa_worker::calculate_agreement_with_known_phase_set()
{
  outlog << "\nCalculating agreement between phases of Solution Estimate A and the known phase set, allowing for alternate origin choice:" << std::endl;

  // No need to check the validity of xn_a_fp - this is done in the calling routine
  // We use solution A here - all Fourier space agreement measures are calculated using this

  // Get measured amplitudes, appropriately apodized, and calculated phases into a common object

  clipper::HKL_data<clipper::data64::F_phi> apodized_f_calculated_phase(hkl_target_eff_res);

  apodized_f_calculated_phase.compute(working_f, known_phase_set, clipper::data64::Compute_fphi_from_fsigf_phifom());

  bool allow_origin_shifts = true;
  clipper::Coord_frac origin_shift(0.0, 0.0, 0.0);

  if (!space_group_is_chiral)
  {
    if (settings.verbose)
    {
      outlog << "\nWithout Inversion: \n"
             << std::endl;
    }
  }

  bool invert = false;
  ipa_functions::compute_phase_agreement(xn_a_fp, apodized_f_calculated_phase, apodization_weights,
                                         settings.res_bins, allow_origin_shifts, invert, stat_xn_mean_phase_diff_all, stat_xn_mean_phase_diff_centric, stat_xn_mean_phase_diff_acentric, stat_xn_map_cc, origin_shift, outlog, settings.verbose);

  if (!space_group_is_chiral)
  {
    if (settings.verbose)
    {
      outlog << "\nWith Inversion: \n"
             << std::endl;
    }
    invert = true;
    ipa_functions::compute_phase_agreement(xn_a_fp, apodized_f_calculated_phase, apodization_weights,
                                           settings.res_bins, allow_origin_shifts, invert, stat_xn_mean_phase_diff_all_inverse, stat_xn_mean_phase_diff_centric_inverse, stat_xn_mean_phase_diff_acentric_inverse, stat_xn_map_cc_inverse, origin_shift, outlog, settings.verbose);
  }
  else
  {
    stat_xn_mean_phase_diff_all_inverse = 0.0;
    stat_xn_mean_phase_diff_centric_inverse = 0.0;
    stat_xn_mean_phase_diff_acentric_inverse = 0.0;
    stat_xn_map_cc_inverse = 0.0;
  }
}

void ipa_worker::copy_F_sigF_data(clipper::HKL_data<clipper::data64::F_sigF> &copy,
                                  clipper::HKL_data<clipper::data64::F_sigF> &paste)
{
  /*
   Helper function that copies data from one object to another, but in a way that is always valid via hkl lookup. Importantly, this probably only works in a high to low resolution direction. And might throw an exception if you try copy low res stuff into a high res data object.

   i.e. we iterate over the hkl indices of the pasted object, and lookup the appropriate hkl() in the copied object. They have different hkl_info associations.
   */
  clipper::HKL_data_base::HKL_reference_index ih;
  for (ih = paste.first(); !ih.last(); ih.next())
  {
    paste[ih] = copy[ih.hkl()];
  }
}

void ipa_worker::report_hkl_info_and_data(clipper::HKL_info &info, clipper::HKL_data_base &data64, std::ostream &outlog)
{
  // Reports some basic information about crystallographic data  to the log file.
  // Note we just accept a hkl_data_base object reference as an input here.

  outlog << "\n CRYSTALLOGRAPHIC INFO:\n\n";

  outlog << " Space group " << info.spacegroup().symbol_xhm() << std::endl;

  bool is_chiral = util::space_group_is_chiral(info.spacegroup().spacegroup_number());

  if (is_chiral)
  {
    outlog << " This space group is chiral" << std::endl;
  }
  else
  {
    outlog << " This space group is achiral" << std::endl;
  }

  outlog << " Cell Dimensions " << info.cell().format() << std::endl;
  outlog << " Cell Volume (Angstroms^3) " << info.cell().volume() << std::endl;

  outlog << "\n DATA INFO:\n\n";

  outlog << " Resolution range " << 1.0 / sqrt(data64.invresolsq_range().min()) << " - " << 1.0 / sqrt(data64.invresolsq_range().max()) << std::endl;

  outlog << " Number of data - possible " << info.num_reflections() << " - observed " << data64.num_obs() << std::endl;
}

// ------------------------------------------------------- //
// -------------------- DEVELOPER ------------------------ //
// ------------------------------------------------------- //

// The following simulate X for run are helper functions only for trouble shooting the generation of the LookUpTables without having to perform an entire run. They could be switched off / verbose only/  or gotten rid of as we get to the stage where we no longer require them.
void ipa_worker::simulate_beta_for_run()
{
  std::ofstream filewrite("beta_simulation.txt", std::ofstream::out);

  if (filewrite.is_open())
  {
    outlog << "Writing beta values to file..." << std::endl;

    int total_iterates = beta_look_up_table.size();

    for (int i = 0; i < total_iterates; i++)
    {
      filewrite << i << " " << beta_look_up_table[i] << "\n";
    }

    filewrite.close();

    outlog << "File written! Success." << std::endl;
  }
}

void ipa_worker::simulate_filter_radius_for_run()
{
  std::ofstream filewrite("filter_radius_simulation.txt", std::ofstream::out);

  if (filewrite.is_open())
  {

    outlog << "Writing filter radius values to file..." << std::endl;

    int total_iterates = filter_radius_look_up_table.size();

    for (int i = 0; i < total_iterates; i++)
    {
      filewrite << i << " " << filter_radius_look_up_table[i] << "\n";
    }

    filewrite.close();

    outlog << "File written! Success." << std::endl;
  }
}

void ipa_worker::simulate_solvent_fraction_multiplier_for_run()
{
  std::ofstream filewrite("solvent_fraction_multiplier_simulation.txt", std::ofstream::out);

  if (filewrite.is_open())
  {

    outlog << "Writing solvent fraction multiplier values to file..." << std::endl;

    int total_iterates = solvent_fraction_multiplier_look_up_table.size();

    for (int i = 0; i < total_iterates; i++)
    {
      filewrite << i << " " << solvent_fraction_multiplier_look_up_table[i] << "\n";
    }

    filewrite.close();

    outlog << "File written! Success." << std::endl;
  }
}

bool ipa_worker::copy_centric_or_acentric_phases_only(clipper::HKL_data<clipper::data64::Phi_fom> source_phases, clipper::HKL_data<clipper::data64::F_phi> target_fp, bool copy_centric, bool copy_acentric)
{
  /*
   Copies the centric and/or acentric data from one dataset to another - this includes both phi and fom.

   This is curremtly just used for testing whether solution of structures in projection might help with the general probem of phase retrieval

   Determination of phases of the centric reflections (always 0 or π) is significantly easier than the acentric case.

   it's assumed both the the source and the target share the same hkl list which is a bit unsafe.

   */
  for (clipper::HKL_data_base::HKL_reference_index ih = hkl_target_eff_res.first(); !ih.last(); ih.next())
  {
    if (ih.hkl_class().centric() == copy_centric)
    { // Copy centric reflections if copy centric is true;
      target_fp[ih].phi() = source_phases[ih].phi();
    }
    else if (ih.hkl_class().centric() != copy_acentric)
    { // copy acentric reflections if copy_acentric isn't false.
      target_fp[ih].phi() = source_phases[ih].phi();
    }
  }

  // target_fp.compute( apodized_f, calculated_phase_set, clipper::data64::Compute_fphi_from_fsigf_phifom() );

  return true;
}

void ipa_worker::write_fingerprint_to_log(clipper::HKL_data<clipper::data64::F_phi> &xn_fp, clipper::HKL_data<clipper::data64::F_phi> &xn_a_fp, clipper::HKL_data<clipper::data64::F_phi> &xn_b_fp, clipper::HKL_data<clipper::data64::F_sigF> &apodized_f_work_set, clipper::Xmap<float> &xn, clipper::Xmap<float> &xn_a, clipper::Xmap<float> &xn_b, std::ofstream &outlog)
{
  /*
   Must be called during iterate to make sure things are declared,
   This is a DEBUG and developer only function.
   */
  outlog << "###############################################################" << std::endl;
  outlog << "##################### DEBUG FINGERPRINT #######################" << std::endl;
  outlog << "###############################################################" << std::endl;

  outlog << "\n Real space representation (xn, xn_a, xn_b) - First 10 data points:" << std::endl;
  outlog << "    n       (xn)        (xn_a)        (xn_b)  " << std::endl;
  int i = 0;
  for (clipper::Xmap_base::Map_reference_index ix = xn.first(); i <= 10; ix.next())
  {

    outlog << std::setprecision(8) << std::setw(5) << ix.index() << " " << std::setw(10) << xn[ix] << " " << std::setw(10) << xn_a[ix] << " " << std::setw(10) << xn_b[ix] << " " << std::fixed << std::endl;

    i++;
  }

  outlog << "\n Fourier space representation (xn_fp||work_fp) - First 20 Observations:" << std::endl;
  outlog << "    n     h     k     l             |F|          Phi             |Fobs|(after apodization)     Starting Phi" << std::endl;
  for (clipper::HKL_data_base::HKL_reference_index ih = xn_fp.first(); ih.index() <= 20; ih.next())
  {
    outlog << std::setprecision(8) << std::setw(5) << ih.index() << " " << std::setw(5) << ih.hkl().h() << " " << std::setw(5) << ih.hkl().k() << " " << std::setw(5) << ih.hkl().l() << " " << std::setw(20) << xn_fp[ih].f() << " " << std::setw(10) << xn_fp[ih].phi() << " " << std::setw(20) << apodized_f_work_set[ih].f() << " " <<
        // std::setw(20) << starting_test_set[ih].phi() << " " <<
        // std::setw(20) <<  test_set[ih].flag() <<
        std::fixed << std::endl;
  }

  outlog << "\n Fourier space representation (xn_A_fp) - First 20 Observations:" << std::endl;
  outlog << "    n     h     k     l             |F|          Phi" << std::endl;
  for (clipper::HKL_data_base::HKL_reference_index ih = xn_a_fp.first(); ih.index() <= 20; ih.next())
  {
    outlog << std::setprecision(8) << std::setw(5) << ih.index() << " " << std::setw(5) << ih.hkl().h() << " " << std::setw(5) << ih.hkl().k() << " " << std::setw(5) << ih.hkl().l() << " " << std::setw(20) << xn_a_fp[ih].f() << " " << std::setw(10) << xn_a_fp[ih].phi() << " " << std::fixed << std::endl;
  }

  outlog << "\n Fourier space representation (xn_B_fp) - First 20 Observations:" << std::endl;
  outlog << "    n     h     k     l             |F|          Phi" << std::endl;
  for (clipper::HKL_data_base::HKL_reference_index ih = xn_b_fp.first(); ih.index() <= 20; ih.next())
  {
    outlog << std::setprecision(8) << std::setw(5) << ih.index() << " " << std::setw(5) << ih.hkl().h() << " " << std::setw(5) << ih.hkl().k() << " " << std::setw(5) << ih.hkl().l() << " " << std::setw(20) << xn_b_fp[ih].f() << " " << std::setw(10) << xn_b_fp[ih].phi() << " " << std::fixed << std::endl;
  }

  outlog << "###############################################################" << std::endl;
  outlog << "###############################################################" << std::endl;
  outlog << "###############################################################" << std::endl;
}

void ipa_worker::substitute_trusted_phases(clipper::HKL_data<clipper::data64::F_phi> &fp)
{
  // Updates the phases in the data object fp with phases present in the data object trusted_phase_set
  // trusted_phase_set is defined during initialisation and contains any phases considered "reliable", currently just used for testing and experimentation

  clipper::HKL_data_base::HKL_reference_index ih;
  for (ih = fp.first(); ih.index() <= 20; ih.next())
  {
    if (!trusted_phase_set[ih.hkl()].missing()) // ensure we only copy data that is present.
    {
      fp[ih].phi() = trusted_phase_set[ih.hkl()].phi(); // set the phase to the trusted value.
    }
  }
}

float ipa_worker::update_phase_consistency(clipper::HKL_data<clipper::data64::F_phi> &current_soln, clipper::HKL_data<clipper::data64::Phi_fom> &current_apodization_weights)
{
  /*
  This experimental function takes in the Fourier coefficients of the current solution estimate, and calculates some real and Fourier space consistency measures 
  with the solution estimate from n iterations prior. 
  
  Hence it is calculating some measures of the distance moved over n iterations.
  
  By making a CUSUM plot of some consistency measure, and looking at its gradient, we can get an indicator of the successful location of a solution
   
  */
  // Set consistency measures to -1, implying that they have not yet been calculated. This is just to enable clean reporting in the log file
  stat_phase_consistency = -1; // 
  stat_mask_consistency = -1;  // 
  stat_phase_cusum = -1;
  stat_real_space_consistency = -1; 

  float mean_distance;   // The mean distance travelled in the complex plane - currently unused.
  float phi_mean_diff;   // The mean absolute phase difference.
  double map_cc;         // The real space map correlation coefficient
  float new_cusum_value;


  if (previous_data_objects_are_valid)
  {

    clipper::Coord_frac orig;
    float mask_cc, smc, ccoef, mcc;

    // ipa_functions::compute_mask_agreement(current_mask, previous_mask_tn, orig, false, mask_cc, smc, ccoef, mcc, outlog, false);
    stat_mask_consistency = -1; // TURNED OFF BINARY MASK CHECK FOR NOW - THIS HAS POTENTIAL TO CAUSE SEGFAULT DURING PHASE DETERMINATION DUE TO CHANGES IN APODIZATION

    ipa_functions::compute_phase_total_difference(current_soln, current_apodization_weights, previous_soln_tn, previous_apodization_weights_tn, mean_distance, map_cc, phi_mean_diff, false, outlog, false);
    
    // stat_phase_consistency = ipa_functions::calculate_circular_correlation_coeff(current_soln, previous_soln_tn); // The circular correlation coefficient provides an alternative to the mean absolute phase difference - unused for now ...
    
    stat_phase_consistency = phi_mean_diff;
    
    stat_real_space_consistency = map_cc;

    // This constrols what statistic we are using for the CUSUM plot
    new_cusum_value = phi_mean_diff;
    phi_difference_plot.push_back(new_cusum_value);

    // Now that we have a new value for the consistency, we can update the CUSUM plot...
    // Calculate new mean, based off all previous mean phase difference values
    double new_mean = 0.0;
    for (int i = 0; i < phi_difference_plot.size(); i++)
    {
      new_mean += phi_difference_plot[i];
    }
    new_mean /= phi_difference_plot.size();

    float old_cusum = 0.0;
    if (phi_CUSUM_plot.size() > 0)
    {
      old_cusum = phi_CUSUM_plot.back();
    }

    float new_cusum = old_cusum + (new_mean - new_cusum_value);
    phi_CUSUM_plot.push_back(new_cusum);
    stat_phase_cusum = new_cusum;
  }

  if (phi_CUSUM_plot.size() > 10)
  {
    // Calculate a gradient.
    // store in a stat
  }

  if (!stat_predicted_success)
  {
    // Use gradient information to predict success
    // stat_predicted_success = ipa_functions::predict_success(list_of_phases, stat_successful_iterate, 5);
  }

  // Store the current solution estimate, apodization weights and envelope, so they are ready for the next time this function is called.

  clipper::HKL_info::HKL_reference_index ih;

  previous_soln_tn.init(current_soln.hkl_info(), current_soln.hkl_info().cell());
  previous_apodization_weights_tn.init(current_apodization_weights.hkl_info(), current_apodization_weights.hkl_info().cell());

  // store the Fourier coefficients of the solution estimate
  for (ih = current_soln.hkl_info().first(); !ih.last(); ih.next())
  {
    previous_soln_tn[ih].f() = current_soln[ih].f();
    previous_soln_tn[ih].phi() = current_soln[ih].phi();
  }

  // store the apodization weights
  for (ih = current_apodization_weights.hkl_info().first(); !ih.last(); ih.next())
  {
    previous_apodization_weights_tn[ih].phi() = 0; // The phase has no utility here, this object is simply used to store the apodization weights
    previous_apodization_weights_tn[ih].fom() = current_apodization_weights[ih].fom();
  }

  // Store the molecular envelope.
  previous_mask_tn.init(current_mask.spacegroup(), current_mask.cell(), current_mask.grid_sampling());
  clipper::Xmap_base::Map_reference_index ix;
  for (ix = current_mask.first(); !ix.last(); ix.next())
  {
    previous_mask_tn[ix] = current_mask[ix];
  }

  previous_data_objects_are_valid = true;

  return stat_phase_consistency;
}
