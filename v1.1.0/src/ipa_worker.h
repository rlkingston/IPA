// IPA: Iterative Projection Algorithms for protein crystallography
// Header for the IPA worker class

#pragma once

#include <iostream>
#include <fstream>   // file streams
#include <iomanip>   // std::setw & std::setprecison & std::boolalpha
#include <algorithm> // std::sort & std::min
#include <cmath>     // lround, lroundf, std::abs, std::pow
#include <time.h>    // clock
#include <random>    // random number generation
#include <atomic>    // For parsing between threads.

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

#include "ipa-lib.h"
#include "svn-version.h"

// ------------------------------------------------------- //
// ----------------------- FORTRAN ----------------------- //
// ------------------------------------------------------- //

// Fortran 2008 subroutines required for computing apodization functions
extern "C"
{
  void get_gaussian_width(double *smax, double *target_integral, double *sigma);
}
extern "C"
{
  double gaussian_integral(double *sigma, double *smax);
}

// Fortran 90 subroutines required for local regression ... used for smoothing histograms.
extern "C"
{
  void lowess(float *x, float *y, int *n, float *f, int *nsteps, float *delta, float *ys, float *rw, float *res);
}

// Fortran 90 function required for calculation of the Incomplete Gamma Integral
extern "C"
{
  double gammad(double *x, double *p, int *ifault);
}

// ------------------------------------------------------- //
// ----------------------- STRUCTS ----------------------- //
// ------------------------------------------------------- //

struct pdb_data
{
  clipper::String pdb_id;
  float mean_bfactor; // TODO ... change this globally to "overall_bfactor" reflecting that the overall B factor will be systematically higher than the mean B factor when there is a distribution of B values in the cell.
  float resolution_limit;
};

// Statistics are organised with structs based on use
// Makes for ease of data handling with std::vector
// member variables to be easily accessed within object.
// Are updated throughout each run.
struct real_space_statistic
{
  float xn_solvent_variance;
  float xn_solvent_mean;
  float xn_solvent_minimum;
  float xn_solvent_maximum;
  float xn_protein_variance;
  float xn_protein_mean;
  float xn_protein_minimum;
  float xn_protein_maximum;
  float xn_histogram_agreement;
};

struct fourier_space_statistic
{
  float xn_correlation_coeff_work;
  float xn_r_factor_work;
  float xn_correlation_coeff_free;
  float xn_r_factor_free;
};

struct test_case_phase_statistic
{
  float xn_mean_phase_diff_all;
  float xn_mean_phase_diff_centric;
  float xn_mean_phase_diff_acentric;
  float xn_mean_phase_diff_all_inverse;
  float xn_mean_phase_diff_centric_inverse;
  float xn_mean_phase_diff_acentric_inverse;
};

// TODO test case envelope agreement statistics also belong here  ??

#ifndef IPA_WORKER
#define IPA_WORKER

class ipa_worker
{

public:
  // ------------------------------------------------------- //
  // ----------------------- Methods ----------------------- //
  // ------------------------------------------------------- //

  ipa_worker(ipa_settings new_settings); // Constructor, worth adding a settings argument in the future?

  // This is the only thing that should be called externally to this object.
  bool run_experiment(std::atomic<float> &ThreadProgress,
                      std::atomic<int> &IterateCounter,
                      std::atomic<float> &CUSUM,
                      const int run_id);

  // Debug 'simulation' functions, can also be called externally when the correct settings have been applied.
  void simulate_beta_for_run();

  void simulate_filter_radius_for_run();

  void simulate_solvent_fraction_multiplier_for_run();

  // ------------------------------------------------------- //
  // ----------------------- Variables --------------------- //
  // ------------------------------------------------------- //

private:
  // ------------------------------------------------------- //
  // ----------------------- Methods ----------------------- //
  // ------------------------------------------------------- //

  bool initialise_run();

  bool run_loop(std::atomic<float> &ThreadProgress, const IPA_updaterule_command my_update_rule, const IPA_constraint_command my_constraints);

  bool finish_run();

  void populate_reference_data();

  void calculate_filter_radius_for_iterate(const int n_iterates, const bool lerp);

  void calculate_ipa_beta_for_iterate(int n_iterates, bool lerp);

  void calculate_solvent_fraction_multiplier_for_iterate(const int n_iterates,
                                                         const bool found_solvent_fraction_multiplier_params,
                                                         const bool lerp);

  void calculate_apodization_look_up_table();

  void check_for_apodization(int iterate);

  void perform_fourier_space_projection(const clipper::Xmap<float> &map,
                                        clipper::HKL_data<clipper::data64::F_phi> &map_fp,
                                        clipper::Xmap<float> &map_out);

  void perform_real_space_projection(const clipper::Xmap<float> &map_in,
                                     clipper::Xmap<float> &map_out,
                                     const IPA_constraint_command my_constraints);

  void difference_map_algorithm(const IPA_constraint_command &my_constraints);

  void relax_relax_reflect_algorithm(const IPA_constraint_command &my_constraints);

  void reversed_relax_relax_reflect_algorithm(const IPA_constraint_command &my_constraints);

  void error_reduction_algorithm(const IPA_constraint_command &my_constraints);

  void calculate_how_many_iterations_for_this_apodization_step();

  void remove_low_resolution_data(clipper::HKL_data<clipper::data64::F_sigF> &data,
                                  clipper::HKL_data<clipper::data64::Flag> &datacutoff);

  bool safely_import_mask(clipper::Xmap<float> &mask,
                          const std::string &mask_file,
                          const clipper::HKL_info &hkl_check);

  bool select_reference_structure();

  bool initialise_data_objects();
  bool read_reference_mtz();
  bool define_working_and_test_sets();

  bool define_reference_density_histograms();

  bool initialise_apodization_scheme();

  bool generate_a_starting_phase_set();

  bool apodize_data();

  bool regrid_working_maps();

  bool calculate_reference_histograms();

  void check_if_histogram_or_gradient_matching();

  bool generate_a_mask(const double current_filter_radius);

  void summarise_iterate_statistics(std::ostream &outstream);

  void calculate_real_space_agreement_statistics();
  void calculate_fourier_space_agreement_statistics();

  void calculate_protein_solvent_histogram_density_statistics();

  void calculate_agreement_with_known_mask();

  void calculate_agreement_with_known_phase_set();

  inline void substitute_Fourier_amplitudes(clipper::HKL_data<clipper::data64::F_phi> &fp,
                                            const clipper::HKL_data<clipper::data64::F_sigF> &measured_f,
                                            const double &p_threshold,
                                            const bool &replace_largeI_with_square_of_expectedF,
                                            const double &intensity_multiplier,
                                            const bool &remove_low_res_data,
                                            const double low_res_cutoff,
                                            const int &res_bins,
                                            std::ostream &outlog,
                                            const bool verbose = false);

  void copy_F_sigF_data(clipper::HKL_data<clipper::data64::F_sigF> &copy,
                        clipper::HKL_data<clipper::data64::F_sigF> &paste);

  void report_hkl_info_and_data(clipper::HKL_info &info,
                                clipper::HKL_data_base &data64,
                                std::ostream &outlog);

  void write_fingerprint_to_log(clipper::HKL_data<clipper::data64::F_phi> &xn_fp,
                                clipper::HKL_data<clipper::data64::F_phi> &xn_a_fp,
                                clipper::HKL_data<clipper::data64::F_phi> &xn_b_fp,
                                clipper::HKL_data<clipper::data64::F_sigF> &apodized_f_work_set,
                                clipper::Xmap<float> &xn,
                                clipper::Xmap<float> &xn_a,
                                clipper::Xmap<float> &xn_b,
                                std::ofstream &outlog);

  bool set_starting_phase_set_from_phi(const clipper::HKL_data<clipper::data64::F_phi> input_phi);

  bool copy_centric_or_acentric_phases_only(clipper::HKL_data<clipper::data64::Phi_fom> source,
                                            clipper::HKL_data<clipper::data64::F_phi> target,
                                            bool copy_centric = true,
                                            bool copy_acentric = false); // defaults to copying centric only.

  void substitute_trusted_phases(clipper::HKL_data<clipper::data64::F_phi> &map_fp);

  void rebuild_hkl_target_objects();

  void prepare_data_objects(const int &current_apo_step);

  void define_working_resolution(float low_resolution_cutoff, float high_resolution_cutoff);

  float calculate_resolution_for_apodization_step(int apodization_step, int nth_height);

  void define_apodization_weights(const double sigma);

  bool initialise_starting_map();

  float update_phase_consistency(clipper::HKL_data<clipper::data64::F_phi> &current_phi, clipper::HKL_data<clipper::data64::Phi_fom> &current_apodization_weights);

  void write_formatted_summary_to_log(std::ostream &outstream);

  // --------------------------------------------------------- //
  // ----------------------- Variables ----------------------- //
  // --------------------------------------------------------- //

  // These two are for testing purposes only, as they are only populated if a known phase set has been passed.
  bool solution_found;             // Whether or not a solution has been found.
  int solution_found_iterate = -1; // The iterate a solution was found

  ipa_settings settings; // Struct containing all pertinent settings

  int my_run_id; // Instantly set when run_experiment is called, can be referred to for the unique run identifier.

  double expected_mean_density_overall = 0.0;

  std::ofstream outlog;       // Output stream for logging.
  std::ofstream summary_file; // Output condensed summary file for consensus steps.

  // Filename storage
  std::string ref_msk_filename;
  std::string ref_mtz_filename;
  std::string src_code_location;

  // Current actual member variables not in the struct:

  // Algorithm Look Up Tables (LUT):
  std::vector<double> beta_look_up_table;                        // stores the value of beta for every iterate
  std::vector<double> solvent_fraction_multiplier_look_up_table; // stores the value of the solvent fraction multiplier for every iterate
  std::vector<double> filter_radius_look_up_table;               // stores the value of the filter radius for every iterate
  std::vector<int> apodization_iterate_look_up_table;            // stores the iterates on which the apodization function will change.

  bool apodize_on_this_step;                    // Whether or not the apodization function changes on the current iterate.
  bool fix_trusted_phases_this_iterate = false; // Whether or not we fix phases from a trusted source, for the current iterate
  int apodization_step;                         // The current apodization step. (Updates)
  bool space_group_is_chiral;

  // These two values store the effective resolution on the last iterate, and the current iterate, if they are the same, we don't need to regrid maps or adjust reflection lists
  double prev_resolution_iter = 0; // Must be initialised, else uh oh. (set second AFTER test if same as current)
  double current_resolution_iter;  // This is set FIRST based on apo LUT.
  // TODO: The below will superceed the above values in short order.
  bool updated_resolution = true; // When set to true, the effective resolution limit has changed and all data objects are updated to reflect the new limit.
  double prev_rlimit = -1;        // The previous effective resolution limit.

  float target_min_resn;    // Minimum resolution of the input data
  float target_max_resn;    // Maximum resolution of the input data
  float wilson_max_resn;    // Maximum resolution of the model-based <I/ε> estimates. This should always be the same as the maximumum resolution of the image reconstruction. Needed ?
  float rlimit;             // The effective resolution limit at the current apodization step,
  float maximum_resolution; // The maximum resolution of the image reconstruction (which may be greater than the maximum resolution of the input data)

  int current_test_set; // Integer flag for which test set we are currently using.

  // Data stored with the object as member variables, we may still wish to pass them to particular functions (always by reference)
  clipper::HKL_info hkl_target_eff_res; // HKL_info object associated with working hkl list that changes effective resolution over a run.
  clipper::HKL_info hkl_reference;      // HKL_info object associated with the reference dataset.
  clipper::HKL_info hkl_target_max_res; // HKL_info object associated with those data objects at maximum resolution, such as mean I Wilson.
  clipper::HKL_info hkl_target_exp_res; // HKL_info object associated with data objects at the experimental input data resolution.

  clipper::HKL_data<clipper::data64::F_sigF> measured_f;   // The original data
  clipper::HKL_data<clipper::data64::Flag> test_set_flags; // Integer flags used to construct the test and working sets, if requested

  clipper::HKL_data<clipper::data64::F_sigF> working_f;             // The resolution truncated and apodized data
  clipper::HKL_data<clipper::data64::F_sigF> test_set_f;            // The resolution truncated and apodized_data - test set only
  clipper::HKL_data<clipper::data64::F_sigF> work_set_f;            // The resolution truncated and apodized data - work set only
  clipper::HKL_data<clipper::data64::Phi_fom> apodization_weights;  // Stores the weights associated with the Gaussian apodization function, used to apodize the data.
  clipper::HKL_data<clipper::data64::Flag> resolution_cutoff_flags; // Stores the flags used to apply a low-resolution cutoff to the data.

  // Major Fourier space data objects
  clipper::HKL_data<clipper::data64::Phi_fom> starting_phase_set;
  clipper::HKL_data<clipper::data64::Phi_fom> known_phase_set;
  clipper::HKL_data<clipper::data64::Phi_fom> trusted_phase_set; // These trusted phases can come from a variety of places, and can be used as additional Fourier space constraints within the IPA

  // These data objects are used when checking that the reconstructed Fourier amplitudes of any "missing data" have not evolved to physically unrealistic values
  clipper::HKL_data<clipper::data64::I_sigI> mean_i_wilson;         // Carries estimates for <I/ε>, based on Wilson statistics without apodization applied. Estimated once from the data, and then left unchanged
  clipper::HKL_data<clipper::data64::I_sigI> working_mean_i_wilson; // Carries estimates for <I/ε>, based on Wilson statistics, at the current apodization step.

  // Custom data flag, useful in custom debugging.
  clipper::HKL_data<clipper::data64::Flag> custom_flag_01;
  clipper::HKL_data<clipper::data64::Flag> custom_flag_02;
  clipper::HKL_data<clipper::data64::Flag> custom_flag_03;

  clipper::HKL_data<clipper::data64::F_phi> xn_fp;
  clipper::HKL_data<clipper::data64::F_phi> xn_a_fp;
  clipper::HKL_data<clipper::data64::F_phi> xn_b_fp;
  // clipper::HKL_data <clipper::data64::Flag>    resolution_cutoff_flags; // Stores the low resolution cutoff flags for use to apply a resolution cutoff to data.

  // Boolean variables which indicate if a valid real space or Fourier space solution estimate is currently stored in the associated objects. As an extra safety check.
  bool xn_a_is_valid;
  bool xn_b_is_valid;
  bool xn_a_fp_is_valid;
  bool xn_b_fp_is_valid;
  bool xn_fp_is_valid;

  // Reference Structure Objects:
  clipper::HKL_data<clipper::data64::F_phi> measured_fp_reference;
  clipper::HKL_data<clipper::data64::F_phi> apodized_fp_reference;

  // Reference Grids, Reference Maps, and Masks.
  clipper::Grid_sampling grid_reference;

  clipper::Xmap<float> input_reference_mask;

  clipper::Xmap<float> reference_map;
  clipper::Xmap<int> reference_mask;
  clipper::Xmap<float> reference_weights;
  clipper::Xmap<float> reference_map_rescaled;

  // Operating maps:
  clipper::Xmap<float> xn;                // map object to store the real space iterate
  clipper::Xmap<float> xn_a;              // store the first solution estimate (Real space constraints applied last)
  clipper::Xmap<float> xn_b;              // store the second solution estimate (Fourier space constraints applied last)
  clipper::Xmap<int> current_mask;        // integer valued map object to hold the working protein mask
  clipper::Xmap<int> known_mask;          // integer valued map object to hold the known protein mask (test cases only)
  clipper::Xmap<int> known_mask_inverted; // integer valued map object to hold the inversion of the known protein mask (test cases only)
  clipper::Xmap<float> target_weights;    // map object to store the statistical weights for each point in the ASU.

  clipper::Xmap<float> input_target_mask_known;   // the known protein mask as originally input (test cases only)
  clipper::Xmap<float> input_target_mask_working; // a working mask as originally input

  double beta; // This changes over time and is for quick reference.. could scope it...

  // It seems like we must initialize the iterate counters here to avoid issues ... there may be a better place to do this
  int total_iterate_counter = 0;   // This tracks the number of IPA iterations currently executed during a run
  int ap_step_iterate_counter = 0; // This tracks the number of IPA iterations currently executed at a given apodization step
  int n_iterations_ap_step = 0;    // This is the total number of IPA iterations that will be executed at any given apodization step.

  bool histogram_matching_on_current_iterate;
  bool gradient_matching_on_current_iterate;

  histogram density_histogram;
  histogram gradient_histogram;

  histogram density_histogram_reference_raw;
  histogram density_histogram_reference;
  histogram density_histogram_reference_rescaled;

  histogram gradient_histogram_reference_raw;
  histogram gradient_histogram_reference;

  // Holds the summary statistics for the protein region of the reference map, changes depending on what was last called.
  map_region_summary_statistics protein_density_stats_reference;
  // Holds the summary statistics for the protein region of the gradient magnitude reference map, changes depending on what was last called.
  map_region_summary_statistics protein_gradient_stats_reference;

  // The convention for defining the protein and solvent regions in the integer-valued mask: Solvent = 0, Protein >=1.
  // The following flags are used in calls to the density modification routines
  const int solvent_region_flag = 0;
  const int protein_region_flag = 1;
  const int all_regions_flag = -1;

  // Random Number generation
  std::default_random_engine MyGenerator; // Reference to generator we're using.
  std::uniform_real_distribution<double> uniform_distribution;

  // Reference Structures array
  int n_reference_structures;
  int ref_id;
  std::vector<pdb_data> reference_structures;

  // Apodisation variables
  std::vector<double> sigma_apodize; // To replace the [] array thats dynamically allocated at start of a run. This dictates the sigma to use at each apodization step.

  // Variables for storing transformed intensity statistics <I/ε> for the centric and acentric data - unchanging
  std::vector<double> mean_I_acentric_measured_data;
  std::vector<double> mean_I_centric_measured_data;
  std::vector<int> observation_counts_acentric_measured_data;
  std::vector<int> observation_counts_centric_measured_data;

  // Variables for storing transformed intensity statistics <I/ε> for the centric and acentric data with apodization applied
  std::vector<double> mean_I_acentric_apodized_data;
  std::vector<double> mean_I_centric_apodized_data;
  std::vector<int> observation_counts_acentric_apodized_data;
  std::vector<int> observation_counts_centric_apodized_data;

  // Variables for determining internal success metrics:
  std::vector<float> phi_difference_plot; // Values of every phase difference calculated, used for calculating a mean for CUSUM.
  std::vector<float> phi_CUSUM_plot;      // Values of the CUSUM plot, based on the mean phase difference between current and previous iterates.
  bool stat_predicted_success = false;    // Whether or not we have predicted if we have found a solution or not.
  int stat_successful_iterate = -1;       // The predicted iterate at which a solution was found.

  // Calculations of phase and envelope consistency
  bool previous_data_objects_are_valid = false;                                // Indicates if structure factor amplitudes, envelope etc have been stored for the consistency checks.
  clipper::HKL_data<clipper::data64::F_phi> previous_soln_tn;                  // Stores previous structure factor estimates
  clipper::HKL_data<clipper::data64::Phi_fom> previous_apodization_weights_tn; // Stores previous weights associated with the Gaussian apodization function
  clipper::Xmap<int> previous_mask_tn;                                         // Stores previous mask

  // Stats for the current iterate are stored individually, and dumped to the output text file at the end of each iterate (if requested).
  bool write_stats_this_iterate;

  float stat_xn_solvent_mean;
  float stat_xn_solvent_variance;
  float stat_xn_solvent_minimum;
  float stat_xn_solvent_maximum;
  float stat_xn_protein_mean;
  float stat_xn_protein_variance;
  float stat_xn_protein_minimum;
  float stat_xn_protein_maximum;
  float stat_xn_histogram_agreement;
  float stat_xn_correlation_coeff_work;
  float stat_xn_r_factor_work;
  float stat_xn_correlation_coeff_free;
  float stat_xn_r_factor_free;
  float stat_delta_n;
  float stat_xn_radius;
  float stat_xn_fraction_solvent;
  float stat_xn_mean_phase_diff_all;
  float stat_xn_mean_phase_diff_centric;
  float stat_xn_mean_phase_diff_acentric;
  float stat_xn_map_cc;
  float stat_xn_mean_phase_diff_all_inverse;
  float stat_xn_mean_phase_diff_centric_inverse;
  float stat_xn_mean_phase_diff_acentric_inverse;
  float stat_xn_map_cc_inverse;
  float stat_mask_cc;
  float stat_smc;
  float stat_cohen_kappa;
  float stat_mcc;
  float stat_mask_cc_inverse;
  float stat_smc_inverse;
  float stat_cohen_kappa_inverse;
  float stat_mcc_inverse;
  float stat_xn_beta;

  // These are updated as predictive measures of success.
  float stat_phase_consistency;      // The mean absolute phase difference between the current and previous iterate (based on settings.phase_consistency_interval as a period).
  float stat_mask_consistency;       // The correlation coefficient for the binary mask.
  float stat_real_space_consistency; // The real-space correlation coefficient.
  float stat_phase_cusum;            // The cumulative sumation plot for success.
};

#endif
