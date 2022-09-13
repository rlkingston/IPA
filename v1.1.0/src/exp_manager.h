// IPA: Iterative Projection Algorithms for protein crystallography
// Header for the experiment manager class

#pragma once

#include <iostream>
#include <fstream> // file streams
#include <iomanip> // std::setw & std::setprecison & std::boolalpha
#include <algorithm> // std::sort & std::min
#include <cmath> // lround, lroundf, std::abs, std::pow
#include <time.h> // clock
#include <random> // random number generation
#include <thread> // Good for sleeping thread.
#include <sys/stat.h>

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

#include "ipa_manager.h"
#include "util-lib.h"
#include "ipa-lib.h"
#include "consensus_worker.h"
#include "svn-version.h"

#ifndef EXP_MANAGER
#define EXP_MANAGER

class exp_manager
{
public:
  
  // ------------------------------------------------------- //
  // ---------------- PUBLIC Methods ----------------------- //
  // ------------------------------------------------------- //
  
  bool run(std::string experimentparams); // All we need.
  
  inline b_factor_data estimate_B_from_diffraction_data(clipper::HKL_data<clipper::data64::F_sigF>& measured_f, std::ofstream& bfactorlog);
  
  // ------------------------------------------------------- //
  // ---------------- PUBLIC Variables --------------------- //
  // ------------------------------------------------------- //
  ipa_settings envelope_settings; // Struct containing settings for envelope determinatiom
  ipa_settings phasing_settings; // Struct containing settings for phase detrmination
  consensus_mask_settings envelope_consensus_mask_settings; // Struct containing settings for consensus envelope determination
  consensus_phase_settings phase_consensus_phase_settings; // Struct containing settings for consensus phase determination 
  
  // Major directories for job specific locations (i.e. set to be within /job_name/)
  std::string dir_envelope_determination;
  std::string dir_phase_determination;
  std::string dir_envelope_consensus;
  std::string dir_phase_consensus;
  std::string dir_final_outputs;
  
  // Unique Job identifier for this experiment.
  std::string job_name;
  
  // Input mtz file path & name - also stored in ipa_settings.
  std::string mtz_input_file;
  
  int n_threads;
  
private:
  
  // ------------------------------------------------------- //
  // ----------------------- Methods ----------------------- //
  // ------------------------------------------------------- //
  
  // JOB methods:
  bool run_initialisation();
  bool run_envelope_determination_experiment();
  bool run_envelope_consensus();
  bool run_phase_determination_experiment();
  bool run_phase_consensus();

  bool validate_settings(ipa_settings& setting);
  
  int import_settings_from_file(ipa_settings& env_settings,
                                ipa_settings& pha_settings,
                                std::string filename,
                                std::ostream& outlog);
  
  // Helper functions:
  void output_settings_summary(std::ostream& outlog);
  
  bool unpack_settings_as_log(std::ostream& outlog,
                                  ipa_settings& settings);
  bool unpack_settings_as_log(std::ostream& outlog,
                                  consensus_mask_settings& settings);
  bool unpack_settings_as_log(std::ostream& outlog,
                                  consensus_phase_settings& settings);
                                                  
  std::string return_algorithm_name_from_enum(IPA_updaterule_command myalg);
  
  void lower_string(std::string& string);
  
  void remove_comments_from_string(std::string& line);
  
  void parse_string_into_string_array(std::vector<std::string>& stringArray,
                                      const std::string& input,
                                      const std::vector<char>& delims);
  
  bool read_string_setting_into_struct(ipa_settings& SettingsToUpdate,
                                       std::vector<std::string>& info,
                                       std::ostream& outlog);
  
  bool read_mask_consensus_setting_into_struct(std::vector<std::string>& info,
                                               std::ostream& outlog);
  
  bool read_phase_consensus_setting_into_struct(std::vector<std::string>& info,
                                                std::ostream& outlog);
  
  bool read_string_algorithm_info_into_struct(ipa_settings& settingstruct,
                                              std::vector<std::string>& info,
                                              std::ostream&outlog);
  
  bool read_string_with_core_label(std::vector<std::string>& lines,
                                    std::ostream& outlog);
  
  bool read_job_list_from_file(std::vector<std::string>& line, std::ostream& outlog);
  
  bool convert_string_into_bool(std::string input);
  
  // Parsing sequential parameter files within experiment.params:
  bool read_string_beta_params_info_into_struct(std::vector<beta_param>& writeparams,
                                                std::vector<std::string>& lines,
                                                std::ostream& outlog);
  bool read_string_solvent_params_info_into_struct(std::vector<solvent_fraction_param>& writeparams,
                                                   std::vector<std::string>& lines,
                                                   std::ostream& outlog);
  bool read_string_filter_params_info_into_struct(std::vector<filter_param>& writeparams,
                                                  std::vector<std::string>& lines,
                                                  std::ostream& outlog);
  
  // Helper functions that will generate default examples for all input parameter files.
  void generate_default_experiment_params_file(std::string paramname);
  
  bool validate_input_mtz(); // Simply read the input mtz.
  
  void initialise_hkl_info_and_data_objects_for_settings(ipa_settings& settings); // Copy the member variable data into the ipa_settings struct.
  
  void print_summary_information_page(std::ostream& custom_outlog);
  
  std::ostream& summarise_cluster_data(cluster_data& data, std::ostream& outlog);

  // ------------------------------------------------------- //
  // ----------------------- Variables --------------------- //
  // ------------------------------------------------------- //
  
  std::ofstream summary_log; // Log file for writing the summary of the whole job.
  std::string summary_log_filename; // Saved string on log file create, thus allowing easy opening and closing between C++ and fortran.

  int cl = 78; // Custom Length of printed outputs.
  
  // Related to calculating the exact atomic composition of the crystal.
  const int n_elements = 98; // Any change to this number must be mirrored in the declaration section of fortran module atom_count
  
  char protein_sequence_filename[128];
  char dna_sequence_filename[128];
  char rna_sequence_filename[128];
  
  int protein_sequence_ncopies_in_asu = 1;
  int dna_sequence_ncopies_in_asu = 1;
  int rna_sequence_ncopies_in_asu = 1;
  
  int protein_sequence_nchains = 1;
  int dna_sequence_nchains = 1;
  int rna_sequence_nchains = 1;
  
  bool protein_sequence_is_provided = false;
  bool rna_sequence_is_provided = false;
  bool dna_sequence_is_provided = false;
  
  // Related to debugging settings
  int n_bad_inputs_count;
  std::string my_param_name;
  
  // Job TODO flags:
  bool JOB_initialisation;
  bool JOB_envelope_determination;
  bool JOB_envelope_consensus;
  bool JOB_phase_determination;
  bool JOB_phase_consensus;
  bool JOB_custom;
  
  // Flags to enable tracking of essential program input
  bool found_job_list;
  bool found_beta_params;
  bool found_filter_params;
  bool found_solvent_radius_params;
  bool core_b_factor_provided;
  bool core_solvent_fraction_provided;
  bool custom_job_id;
  bool found_input_file;
  bool found_input_file_columns;
  
  // Container to hold all information related to output of consensus solutions from consensus workers.
  std::vector<cluster_data> list_of_solutions;

  std::vector<double> envelope_consensus_epsilon; // List of epsilon values to use during envelope consensus.

  // Special settings
  std::vector<std::string> list_of_working_masks; // Container to hold the masks to be input into phase determination.
  int user_input_mask_count = 0; // Counter indicating how many of the masks were input via the parameter file.
  int phase_determination_counter = 0; // Counts the number of phase Determination steps performed, for summary output.
  bool infinity_mode = false; // Endlessly calculate phases until solution is found.
  bool no_display = false;
  double high_resolution_truncation = -1.0; // High resolution cut-off to use on data import. If this is non-negative, the cutoff will be applied ...
  clipper::HKL_info hkl_target; // Hkl info objects to handle the data below.
  clipper::HKL_info hkl_wilson;
  clipper::HKL_data <clipper::data64::F_sigF> measured_f_all; // Used to store the user supplied data 
  clipper::HKL_data <clipper::data64::I_sigI>  mean_i_wilson; // Used to store estimates for the mean intensity as a function of hkl : <I(hkl)>. This needs to be copied to the ipa_manager and worker. 
  
  bool mean_i_wilson_is_valid = false; // Defaults to false, to ensure it is calculated at least once.
  
  bool free_lunch = false; // Are we extending the data beyond the experimental limit ?

  // 32 = space, 9 = Horizontal Tab. See https://en.cppreference.com/w/cpp/language/ascii for details.
  std::vector<char> deliminators = {32,9}; 

  b_factor_data target_b_factor_data;
    
};

#endif
