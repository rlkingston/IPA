// IPA: Iterative Projection Algorithms for protein crystallography
// Header for the core function library 

#pragma once

#include "util-lib.h"

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>

#include <Eigen/Eigen>

#include "string"
#include "fstream"

// The following are from consensus_mask
#include <random> // random number generation
#include <sstream>
#include <clipper/clipper-contrib.h>
#include <iomanip>

#ifndef IPA_LIBRARY
#define IPA_LIBRARY

// Forward Declarations:
struct ipa_setting_algorithm_container;
struct beta_param;
struct solvent_fraction_param;
struct filter_param;

struct map_region_summary_statistics
{
  double mean;
  double variance;
  double skewness;
  double kurtosis;
  double minimum;
  double maximum;
};

struct histogram
{
  int    number_of_intervals;
  double minimum;
  double maximum;
  float *interval;
};
// nb the "number of intervals" in the structure histogram is redundant, since it is implict in the dimensioning of "interval"
// a better data structure might be ...
//
// struct histogram
//  {
//   double minimum;
//   double maximum;
//  float *bin_values;
//  };
//
// TODO: ... rework to eliminate this redundancy at some point ???

struct b_factor_data
{
  // Struct variables
  double SCALEK;              // Absolute Scale Factor
  std::vector<double> U_CIF;  // Symmetric tensor U_CIF describing atomic displacement - for definition see see Grosse-Kunstleve & Adams J Appl Cryst. 2002;35(4):477–80.
  std::vector<double> U_STAR; // Symmetric tensor U_STAR describing atomic displacement - for definition see see Grosse-Kunstleve & Adams J Appl Cryst. 2002;35(4):477–80.
  std::vector<double> B_STAR; //  Symmetric tensor B describing atomic displacement B = 2π^2 * U_STAR
  double B_ISO;               // Equivalent isotropic displacement parameter B
  double U_ISO;               // Equivalent isotropic displacement parameter U

  // Null Constructor
  b_factor_data() 
  {
    U_STAR.resize(9);
    U_CIF.resize(9);
    B_STAR.resize(9);
  }

  // Copy Constructor
  b_factor_data(const b_factor_data &copy) : SCALEK(copy.SCALEK),
                                             B_ISO(copy.B_ISO),
                                             U_ISO(copy.U_ISO)
  {
    U_STAR.resize(9);
    U_CIF.resize(9);
    B_STAR.resize(9);

    for (int i = 0; i < 9; i++)
    {
      U_STAR[i] = copy.U_STAR[i];
      U_CIF[i] = copy.U_CIF[i];
      B_STAR[i] = copy.B_STAR[i];
    }
  }

  // Constructor from std::vectors
  b_factor_data(const double &u_iso,
                const double &b_iso,
                const double &scale_k,
                const std::vector<double> &u_star,
                const std::vector<double> &u_cif,
                const std::vector<double> &b_star) : U_ISO(u_iso),
                                                     B_ISO(b_iso),
                                                     SCALEK(scale_k)
  {
    U_STAR.resize(9);
    U_CIF.resize(9);
    B_STAR.resize(9);

    for (int i = 0; i < 9; i++)
    {
      U_STAR[i] = u_star[i];
      U_CIF[i] = u_cif[i];
      B_STAR[i] = b_star[i];
    }
  }

  clipper::U_aniso_orth get_clipper_U_aniso_orth()
  {
    // Note the matrices we store are 3x3, however there are only six unique elements, as the matrix is diagonally symmetric
    // b11, b22, b33, b12, b13, b23
    // or:
    // B11 B12 B13
    //     B22 B23
    //         B33
    return clipper::U_aniso_orth(U_STAR[0], U_STAR[4], U_STAR[8], U_STAR[1], U_STAR[2], U_STAR[5]);
  }
};

struct ipa_settings
{
  // container for holding all settings information required for an IPA run.
  
  // Special CORE data objects, copies made from experiment_manager... Need to check these copies are legit when being passed this way, as the setting struct is copied from experiment manager to ipa manager to worker.
  
  clipper::HKL_info hkl_target; // HKL_info for the data input by the user
  clipper::HKL_info hkl_wilson; // Same as target HKL_info *except* that the resolution may extend beyond the measured limit. 
  clipper::HKL_data <clipper::data64::F_sigF>  measured_f_all; // Stores the unmodified original data
  clipper::HKL_data <clipper::data64::I_sigI>  mean_i_wilson;  // Stores model-based estimates of the mean intensity <I>

  // Pass by const reference to any required functions so that all settings can be accessed in simple manner without passing multiple different values.
  
  bool job_is_envelope_determination = false;
  bool job_is_phase_determination  = false; 
  bool check_for_consensus_phases = false;
  
  std::vector<ipa_setting_algorithm_container> UpdateRuleRegime; // Details the IP algorithms that will be used during the experiment
  
  // Custom parameters to modulate beta / the solvent fraction / the filter radius over the course of a run
  // These are used to calculate the values for beta, the solvent fraction multiplier and the filter radius, as a function of iterate
  std::vector<beta_param> beta_params_look_up_table;
  std::vector<solvent_fraction_param> solvent_fraction_multiplier_params_look_up_table;
  std::vector<filter_param> filter_radius_params_look_up_table;
  
  int n_total_iterates = 0; // Is set to total number of iterates, by exp_manager.
  
  clipper::String input_filename_target_data = "NONE";
  clipper::String target_col_fo              = "NONE";
  clipper::String target_col_pw              = "NONE";
  
  clipper::String target_col_pw_calc         = "NONE"; // To be made redundant, new names are below...
  
  // These are the new centralised way to grab known phases from an existing mtz. (not appended to the input mtz.)
  clipper::String input_filename_known_phase_set = "NONE";
  clipper::String known_phase_set_fp_column_labels = "NONE";

  
  // Note "final" outputs should be put into envelope_determination/outputs/ and labelled with job_name, so that its easy to use multiple jobs in consensus
  std::string job_name = "default_job"; // Helps identify where to save log files etc
  std::string working_dir = "ERRORDIR"; // This is either "envelope_determination/job_name/" OR "phasing_determination/job_name/", and is useful for ensuring all the runs standard outputs happen in one easy to reach place.

  bool input_known_phases = false; // Do we have a known phase set for testing and debugging purposes ?

  bool output_envelopes_and_phases = false; // If true, the program outputs both envelopes *and* phase sets at the end of a run. By default only envelopes are output during envelope determination, and phase sets are output during phase determination

  // These two control the output of phase sets *during* a run ... currently disabled   
  bool output_phases = false;
  int phase_output_interval = 1;
  
  // Settings associated with manipulation of phases for the centric observations. Experimental. 
  bool keep_known_centric_phases = false; // This simply copies centric phases to the random starting phase set at the start of a run.
  float centric_random_percent = 1.0; // Float representing the % of centric phases to randomly convert.
  bool control_for_phase_copying = false; // Activate a control for acentric phase changes only.
  int n_fix_centric_phases = 0; // Fix the centric phases for this number of iterates
  
  int fix_trusted_phases = -1; // Any value greater than 0, will fix the trusted phases for that number of iterates (they will be used as Fourier space constraints within the IPA)
  
  // mask handling
  bool input_working_mask = false;
  clipper::String input_filename_target_mask_working = "NONE"; // This is the current target mask we're using.
  
  bool input_known_mask = false;
  clipper::String input_filename_known_mask = "NONE";
  bool fix_mask = false; // This is deprecated I think, or should be. This meant "fix mask throughout". It does nothing that cannot be achieved with n_iterations_fix_mask > 0
  int n_iterations_fix_mask = 0;
  bool compute_mask_from_local_mean = false;
  bool compute_mask_from_local_variance = true;
  bool compute_mask_from_local_mean_and_variance = false;
  int compute_mask_interval = 1;
  float mask_weighting_factor_mean = 0.25;
  double filter_radius = 8.0; 
  double percentage_radius_expansion = 0.0; // This is obsolete ... does nothing that can't be done with the parameters controlling the the filter radius in experiment.params


  // These two control the output of masks *during* a run ... currently disabled 
  bool output_mask = false;
  int mask_output_interval = 100;
  
  // These indicate if parameters controlling variation of beta / the filter radius / the solvent fraction have been provided.
  bool found_filter_params = false;
  bool found_beta_params = false;
  bool found_solvent_fraction_multiplier_params = false;
  
  // TODO ... Shouldn't beta, the solvent_fraction_multiplier, and the filter radius *all* be defined in the settings, and given appropriate default values ? (0.75, 1.0, 8.0?)
  // Then if there are no parameters given which control their behavior as a function of iterate, we have some constant values which would be used in place ...
  
  // Mask editing
  
  bool erase_islands = false;
  int erase_islands_interval = 10;
  float erase_islands_threshold = 0.10;
  
  bool fill_voids = false;
  int fill_voids_interval = 10;
  float fill_voids_threshold = 0.05;
  
  int envelope_editing_limit = 0;
  
  // Real space density modification
  
  // TODO: All of these should be gone soon - they have been made redundant by the new way of specifying update rules and constraints. 
  bool solvent_flattening = true;
  bool histogram_matching = true;
  int histogram_matching_interval = 1;
  bool gradient_matching  = false;
  int gradient_matching_interval  = 5;
  
  bool apodize_reference_data = true; // Determines if we apodize the reference data to match what's been done to the target data.
  
  
  bool input_starting_phases = false; // Do we have some starting phase estimates ?

  // Controls the evaluation of "Free" statistics (the use of test and work sets on each iteration)
  bool use_test_set = false;
  int n_test_sets = 50;
  
  // Extension beyond the experimental limit of the data ... aka a free lunch. Is there any such thing? 
  bool free_lunch = false;
  float free_lunch_scale_factor = 1.0; // Multiplier by which the volume of the sphere enclosing the data in Fourier space will be increased 
  
  // Target crystal properties

  b_factor_data target_b_factor_data; // Holds the information about scale and B-factor for the target dataset

  double fraction_solvent; // No default value - this is obligatory user input (or the program needs to make an informed guess, which it doesn't at present)
  double expected_mean_density_protein = 0.44;
  double expected_mean_density_solvent = 0.38;
  
  // Settings for hemihedrally twinned data ... this is work in progress
  bool hemihedral_twinning = false; // Indicates if the target data is hemihedrally twinned
  int hemihedral_twin_case = 0; // The possible scenarios for hemihedral twinning will be enumerated in an internal table. This indicates which case we are dealing with. Case zero is untwinned for safety.
  double hemihedral_twin_fraction = 0.5; // The hemihedral twin fraction
  
  // Data apodization and resolution variation
  
  bool apodize = false;
  int n_apodization_steps = 0;         // How many Apodization steps are there ?
  int n_apodization_first_iterate = 0; // Data Apodization starts changing at this iterate (0 is always apodized)
  int n_apodization_last_iterate = 0;  // Data Apodization stops changing at this iterate
  
  double initial_apodization_sigma = 1.83842542; // Gaussian falling to 1/20th of the maximum height at 4.5 Å resolution
  double final_apodization_sigma = 0.40853898; // Gaussian falling to 1/20th of the maximum height at 1.0 Å resolution
  bool upper_apodization_limit = false;
 
  int resolution_gaussian_width_cutoff = 120; // Used in setting the effective resolution of the data set (and hence the generation of hkl lists) 
  double shannon_rate = 1.32; // Controls the gridding of the maps at any effective resolution
  bool estimate_missing_amplitudes_on_first_iterate = false; // If true, then the missing Fourier amplitude data will be replaced with the expected value of F, when generating the initial iterate xn
  
  // core IPA parameters
  bool linear_freq_change = false;
  bool fix_beta = true;  // deprecated I think. If you want constant beta, do it with the parameter settings, plus fixed beta should be the default behavior if there are no beta parameter settings
  int  beta_function = 1;
  bool compute_real_space_agreement_statistics = true; // controls whether we calculate real space agreement statistics when evaluating mask agreement. Probably kill this and always do it 
  
  int report_statistics_interval = 20;

  int n_runs=1; // Number of unique runs.

  int first_run_id_offset = 0; // This can be set to offset over-ridding log files etc.
  
  int n_iterations_shrink = 250;
  
  // phase generation
  bool really_random = true;
  bool fourier_space_phase_generation = true;
  int n_sphere_phase_generation = 1;
  double rad_sphere_phase_generation = 12.0;
  int max_attempts = 20;
  
  // Fourier space projection  .. low resolution cutoff and treatment of missing terms.
  bool remove_low_res_data = false;
  double low_res_cutoff = 25.0;
  double p_threshold = 0.995;
  bool replace_largeI_with_square_of_expectedF = false;
  double intensity_multiplier = 0.6;

  
  int res_bins = 10;
  bool verbose = false;
  
  // Dimensioning of the histograms for density and gradient matching
  int n_density_bins  = 220;
  int n_gradient_bins = 220;

  int phase_consistency_interval = 30; // The interval between iterates by which the phase velocity is calculated.
};


struct consensus_mask_settings
{   // These settings are used to pass to the calculate consensus mask job in exp_manager.
    // They contain all necessary inputs to the function, and will be redundantly set based on the exp manager experiment.params file, utilising either the -envelope keyword, or a custom -envelope_consensus keyword for example...
  std::vector<clipper::String> file_list; // No more saving to a file, we pass the names of files directly.
  unsigned int    min_pts = 5; // DBSCAN parameter MinPoints
  bool            epsilon_is_absolute = false;
  double          epsilon = 4.0; // DBSCAN parameter epsilon specified as a percentile of the k-distance distribution
  int             distance_measure = 2;  // Specifies the correlation-based distance measure to be used during clustering
  bool            erase_islands = true;
  bool            fill_voids = true;
  float           erase_islands_threshold = 0.05;
  float           fill_voids_threshold = 0.05;
  bool            input_known_mask = false;
  clipper::String input_filename_known_mask = "NONE";
  bool            input_Fourier_amplitudes = false;
  clipper::String input_filename_target_data = "NONE";
  clipper::String target_col_fo         = "NONE";
  clipper::String mask_list_file = "DEFAULT";
  bool            apodization = false;
  double          apodization_limit = 4.5;
  bool            shift_to_known_origin = false; // If true, consensus masks will be shifted to the same origin as a known mask, and written out (useful for test cases only). 
  bool            allow_inversion = false; // If true, both the mask and its inverse will be tested when calculating distances for clustering (if the space group is achiral). During averaging, input masks will be inverted as required.  
  double          fraction_solvent_target = -1.0; 
  bool            verbose = false;
};

struct consensus_phase_settings
{
  std::vector<clipper::String> file_list; // No more saving to a file, we pass the names of files directly.
  // note filelist above will soon be obsolete due to consensus worker.
  std::string working_dir;
  
  unsigned int min_pts = 2; // Need at least two highly consistent phase sets to form a cluster.
  double epsilon_input;
  clipper::String col_fp    = "/*/*/[F,PHIC]"; // This is the default output for our phasing routines.
  bool apodization = false;  
  float effective_resolution_limit = 5.0;
  int res_bins = 25;
  bool threshold_specified = false;
  bool input_known_phases = false;
  clipper::String input_filename_known_phase_set = "NONE";
  clipper::String known_phase_set_fp_column_labels    = "/*/*/[FC,PHIC]"; // This is the defaults for our sigmaa.mtz files.
  
  clipper::String phase_list_file = "DEFAULT";
  bool allow_origin_shifts = true; // by default will search for alternate origins when calculating map agreement.
  bool shift_to_known_origin = true; // If true, consensus phase sets will be shifted to the same origin as a known phase set, and written out (useful for test cases only).
  bool  verbose = false;
};


// Augment MapFilterFn_base from the clipper library with some new map filter functions

//! Tricube-function radial map filter
/*! \ingroup g_mapf */
class MapFilterFn_tricube : public clipper::MapFilterFn_base {
  public:
    //! constructor: takes radius for step function cutoff
    MapFilterFn_tricube( const clipper::ftype& radius ) : radius_( radius ) {}
    //! destructor
    ~MapFilterFn_tricube() {}
    //! evaluate radial step function: (1-abs(r/r0)^3)^3 if inside or 0.0 outside
    clipper::ftype operator() ( const clipper::ftype& radius ) const
      { return (radius<radius_)?pow(1.0-pow(abs(radius/radius_),3),3):0.0; }
  private:
    clipper::ftype radius_;
};

//! Triweight-function radial map filter
/*! \ingroup g_mapf */
class MapFilterFn_triweight : public clipper::MapFilterFn_base {
public:
  //! constructor: takes radius for step function cutoff
  MapFilterFn_triweight( const clipper::ftype& radius ) : radius_( radius ) {}
  //! destructor
  ~MapFilterFn_triweight() {}
  //! evaluate radial step function: (1-(r/r0)^2)^3 if inside or 0.0 outside
  clipper::ftype operator() ( const clipper::ftype& radius ) const
    { return (radius<radius_)?pow(1.0-pow(radius/radius_,2),3):0.0; }
private:
  clipper::ftype radius_;
};

class ipa_functions {
public:
  
  static bool compute_consensus_mask(consensus_mask_settings settings,
                                     std::ostream& outlog,
                                     std::string out_directory,
                                     std::vector<std::string>& concensus_list);
  
  // a simple scaling function
  static void simple_linear_scale (const clipper::HKL_data<clipper::data64::F_sigF> &measured_fourier_amplitudes_work_set,
                                   const clipper::HKL_data<clipper::data64::F_sigF> &measured_fourier_amplitudes_test_set,
                                   const clipper::HKL_data<clipper::data64::Phi_fom> &apodization_weights,
                                   const clipper::HKL_data<clipper::data64::F_phi>  &working_amplitude_and_phase,
                                   const bool &use_test_set,
                                   const double &threshold,
                                   const int &res_bins,
                                   double &scale_factor,
                                   double &correlation_coefficient,
                                   double &r_factor,
                                   double &correlation_coefficient_test_set,
                                   double &r_factor_test_set,
                                   std::ostream &outlog,
                                   bool verbose = false);
  
  // Functions for calculating Fourier cofficients that allow derivates of the electron density to be obtained by FFT
  static void  compute_fourier_coeffs_first_order_deriv(const clipper::HKL_data<clipper::data64::F_phi> &fp,
                                                        clipper::HKL_data<clipper::data64::F_phi> &fp_dpdx,
                                                        clipper::HKL_data<clipper::data64::F_phi> &fp_dpdy,
                                                        clipper::HKL_data<clipper::data64::F_phi> &fp_dpdz);
  
  static void compute_fourier_coeffs_second_order_deriv(const clipper::HKL_data<clipper::data64::F_phi> &fp,
                                                        clipper::HKL_data<clipper::data64::F_phi> &fp_lapl);
  
  static void calculate_gradient(const clipper::HKL_data<clipper::data64::F_phi> &fp,
                                 clipper::Xmap<float>  &gradient_x,
                                 clipper::Xmap<float>  &gradient_y,
                                 clipper::Xmap<float>  &gradient_z,
                                 clipper::Xmap<float>  &gradient_magnitude,
                                 std::ostream &outlog,
                                 bool verbose = false);
  
  // Functions for calculating mask connectivity
  static void mask_connections (const int &mask_id,
                                const bool &apply_translational_symm,
                                const bool &check_diagonals,
                                const float &fractional_threshold,
                                clipper::Xmap<int> &mask,
                                int &n_connected_sets_eliminated,
                                int &n_connected_sets_retained,
                                std::ostream &outlog,
                                bool verbose = false);
  
  static void update_connection_matrix (int &current_anchor_set,
                                        const int &working_anchor_set,
                                        Eigen::MatrixXi &connection_matrix);

  // Function for computing alternate origins
  static void compute_alternate_origins (const int &space_group_number,
                                         const clipper::Grid_sampling &grid,
                                         std::vector<clipper::Coord_frac> &origins);
  
  // Function for computing phased translation function
  static void phased_translation_function (const clipper::HKL_data<clipper::data64::F_phi> &fa,
                                           const clipper::HKL_data<clipper::data64::F_phi> &fb,
                                           const bool &allow_origin_shifts,
                                           clipper::Coord_frac &x,
                                           float &cc);

  // Function for computing intensity statistics
  static void compute_intensity_statistics (const clipper::HKL_data<clipper::data64::F_sigF> &f,
                                            const int &res_bins,
                                            std::vector<double> &mean_I_acentric,
                                            std::vector<double> &mean_I_centric,
                                            std::vector<int> &observation_counts_acentric,
                                            std::vector<int> &observation_counts_centric,
                                            std::ostream &outlog,
                                            bool verbose = false);
  
  // Function for computing weighted agreement between two sets of phases - phases taken from fp1 and fp2, weights from pw
  static void compute_phase_agreement (const clipper::HKL_data<clipper::data64::F_phi>& fp1,
                                       clipper::HKL_data<clipper::data64::F_phi> fp2,
                                       const clipper::HKL_data<clipper::data64::Phi_fom>& pw,
                                       const int &res_bins,
                                       const bool &allow_origin_shifts,
                                       const bool &invert,
                                       float &diff_all,
                                       float &diff_centric,
                                       float &diff_acentric,
                                       float &map_cc,                                       
                                       clipper::Coord_frac &origin_shift,
                                       std::ostream &outlog,
                                       bool verbose = false);
  
  // Function for calculating a mask from an electron density map, based on the locally computed mean or variance
  static void calculate_mask (const bool &compute_mask_from_local_mean,
                              const bool &compute_mask_from_local_variance,
                              const bool &compute_mask_from_local_mean_and_variance,
                              const double &filter_radius,
                              const float &mask_weighting_factor_mean,
                              const clipper::Xmap<float> &map,
                              const float &fracion_solvent,
                              clipper::Xmap<int> &mask,
                              std::ostream &outlog,
                              bool verbose = false);
  
  static void enforce_mask_connectivity (const bool &erase_islands,
                                         const bool &fill_voids,
                                         const float &erase_islands_threshold,
                                         const float &fill_voids_threshold,
                                         clipper::Xmap<int> &mask,
                                         int &n_connected_sets_eliminated_protein,
                                         int &n_connected_sets_retained_protein,
                                         int &n_connected_sets_eliminated_solvent,
                                         int &n_connected_sets_retained_solvent,
                                         std::ostream &outlog,
                                         bool verbose = false);
  
  static void compute_mask_agreement (const clipper::Xmap<int> &mask_a,
                                      const clipper::Xmap<int> &mask_b,
                                      clipper::Coord_frac &origin_shift,
                                      const bool &compute_agreement_by_voxel,
                                      float &mask_cc,
                                      float &simple_matching_coefficient,
                                      float &cohen_coefficient,
                                      float &mcc, std::ostream &outlog,
                                      bool verbose = false);
  
  static clipper::Xmap<int> complement_mask (const clipper::Xmap<int> &input_mask);
  
  static clipper::Xmap<int> shift_mask(const clipper::Xmap<int> &input_mask,
                                       const clipper::Coord_frac &origin_shift);
  
  static void shift_phases(const clipper::HKL_data<clipper::data64::F_phi> &input_fp,
                           clipper::HKL_data<clipper::data64::F_phi> &output_fp,
                           const clipper::Coord_frac &origin_shift);
  
  static clipper::Xmap<int> invert_mask(const clipper::Xmap<int> &input_mask);
  
  static float correlation_based_distance_measure(const float &cc,
                                                  const int &distance_measure);
  
  
  static void perform_ncs_density_averaging (clipper::Xmap<float> &map,
                                             const clipper::Xmap<int> &ncs_avg_mask,
                                             const std::vector< std::vector< clipper::RTop_frac > > &ncs_operations);
  
  static void calculate_summary_statistics(const clipper::Xmap<float> &map,
                                           const clipper::Xmap<int> &mask,
                                           const clipper::Xmap<float> &weights,
                                           const int &region_flag,
                                           const bool &calculate_higher_moments,
                                           double &mean,
                                           double &variance,
                                           double &skewness,
                                           double &kurtosis,
                                           double &minimum,
                                           double &maximum);
  
  static void calculate_histogram(const clipper::Xmap<float> &map,
                                  const clipper::Xmap<int> &mask,
                                  const clipper::Xmap<float> &weights,
                                  const int &region_flag,
                                  histogram &hg,
                                  std::ostream &outlog,
                                  bool verbose = false);
  
  static void  flatten_map_region(clipper::Xmap<float> &map,
                                  const clipper::Xmap<int> &mask,
                                  const int &region_flag,
                                  const float &value);
  
  static void    shift_map_region(clipper::Xmap<float> &map,
                                  const clipper::Xmap<int> &mask,
                                  const int &region_flag,
                                  const float &value);
  
  static void histogram_match_map_region(clipper::Xmap<float> &map,
                                         const clipper::Xmap<int> &mask,
                                         const int &region_flag,
                                         const histogram &observed_histogram,
                                         const histogram &reference_histogram,
                                         const double &a,
                                         const double &b,
                                         const bool &integrate_from_top);
  
  static double earth_movers_distance(const histogram &hg1,
                                      const histogram &hg2);
  
  static void reconstruct_image_from_gradient(const clipper::HKL_info &hkl_info,
                                              const clipper::Xmap<float>  &gradient_x,
                                              const clipper::Xmap<float>  gradient_y,
                                              const clipper::Xmap<float> &gradient_z,
                                              clipper::Xmap<float> output_map);
  
  static void density_modification (const bool &solvent_flattening,
                                    const bool &histogram_matching,
                                    const bool &gradient_matching,
                                    const double  &expected_mean_density_protein,
                                    const double &expected_mean_density_solvent,
                                    histogram &density_histogram,
                                    const histogram &density_histogram_reference,
                                    const map_region_summary_statistics &protein_density_stats_reference,
                                    histogram &gradient_histogram,
                                    const histogram &gradient_histogram_reference,
                                    const map_region_summary_statistics &protein_gradient_stats_reference,
                                    const clipper::Xmap<float> &input_map,
                                    const clipper::Xmap<int> &mask,
                                    const clipper::Xmap<float> &weights,
                                    const clipper::HKL_info &hkl_target,
                                    clipper::Xmap<float> &output_map,
                                    std::ostream &outlog,
                                    bool verbose = false);
  
  static double rmsd_between_maps(const clipper::Xmap<float> &map1,
                                  const clipper::Xmap<float> &map2,
                                  const clipper::Xmap<float> &weights);
  
  static double flattened_sinusoid(const double squareness,
                                   const double phase,
                                   const double integral);
  
  static double pulse_wave(const double pulse_duration,
                           const double phase,
                           const double integral);
    
  static void compute_phase_total_difference(const clipper::HKL_data<clipper::data64::F_phi>& fp1, 
                                            const clipper::HKL_data <clipper::data64::Phi_fom>& apodization_weights1,
                                            const clipper::HKL_data<clipper::data64::F_phi>& fp2, 
                                            const clipper::HKL_data <clipper::data64::Phi_fom>& apodization_weights2,
                                            float& distance_travelled, 
                                            double& mean_cc , 
                                            float& mean_phi_difference, 
                                            bool centrics_only, 
                                            std::ostream& outlog, 
                                            bool verbose);

static double calculate_circular_correlation_coeff(const clipper::HKL_data<clipper::data64::F_phi> &fp1, const clipper::HKL_data<clipper::data64::F_phi> &fp2);

static void skeletonize_map(clipper::Xmap<float>& xn_input, const double& solvent_density, const double& density_cutoff, const int& sqrrad, const double& beta);

static void smooth_filter_map_region(clipper::Xmap<float>& map_to_smooth, // Input map to be smoothed.
                              const clipper::Xmap<int>& mask, // 1 is the effected region flag by default.
                              const double& smoothfilter_radius, 
                              double& smoothed_map_mean,
                              const int region_flag = 1,
                              bool ignore_unaffected_region = true); // Might want to change this to true, but its more expensive.

static void constant_bias_scale_map_region(clipper::Xmap<float>& map, // Input map 
                                    const clipper::Xmap<int>& mask, // 1 is the effected region flag by default.
                                    double& new_map_mean, // The resultant mean, after the operation.
                                    const double& constant, // Constant bias offset to be applied to map
                                    const double& scale,  // scaling to be applied.
                                    const int region_flag = 1);
                                    
static void power_multiply_map_region(clipper::Xmap<float>& map, // Input map 
                              const clipper::Xmap<int>& mask, // 1 is the effected region flag by default.
                              double& new_map_mean,
                              const double& power,
                              const double& bias, 
                              const int region_flag = 1);

static void atomize_map_region(clipper::Xmap<float>& map, const clipper::Xmap<int>& mask, const clipper::Atom_list& list_of_atoms, const clipper::HKL_info& hkl, const double threshold = 0.1, const int region_flag = 1, const double nearest_distance = 2.5);

private:
  
};

#endif 

