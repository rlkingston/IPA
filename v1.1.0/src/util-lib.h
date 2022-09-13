// IPA: Iterative Projection Algorithms for protein crystallography
// Header for the utility function library 

#pragma once

#include "ipa-lib.h"

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>

#include <Eigen/Dense>

#include "string"
#include "fstream"

// The following are from consensus_mask
#include <random> // random number generation
#include <sstream>
#include <clipper/clipper-contrib.h>
#include <iomanip>

// Common Colours and formatting ANSI hacks.
#define cBlack "\u001b[30m"
#define cRed "\u001b[31m"
#define cGreen "\u001b[32m"
#define cYellow "\u001b[33m"
#define cBlue "\u001b[34m"
#define cMagenta "\u001b[35m"
#define cCyan "\u001b[36m"
#define cWhite "\u001b[37m"
#define cReset "\u001b[37m"
#define cBold "\u001b[1m"

#ifndef UTILITY_LIBRARY
#define UTILITY_LIBRARY

struct ipa_settings; // Forward Declaration.

enum IPA_updaterule_command 
{
  DifMap, // Difference Map Algorithm
  ErrRed, // Error Reduction Algorithm
  ReReRe, // Relax Relax Reflect Algorithm
  RevRRR  // Reverse Relax Reflect Relax Algorithm
};

// These are intended to be added to, and updated to be acted upon within ipa_worker::perform_real_space_projection()
// A lot of these are currently highly experimental ...
enum IPA_constraint_command
{
  pHistogramMatch_sFlatten,            // Density value histogram matching of the protein region, flattening of solvent region (The default)
  pGradientMatch_sFlatten,             // Gradient Restoration of the protein region, flattening of solvent region 
  psSkeletonize,                       // Skeletonization of the protein and solvent regions.
  pSkeletonize_sFlatten,               // Skeletonization of the protein region, flattening of solvent region 
  pWeightedSkele_sFlatten,             // Skeletonization of the protein region, merged with weighting to pre-skeleton map, flattening of solvent region 
  pAtomize_sFlatten,                   // Atomization of the protein region, flattening of solvent region
  pGlobicize_sFlatten,                 // Globicization of the protein region, flattening of solvent region
  pHistogramMatch_sSmooth,             // Density value histogram matching of the protein region, smoothing of solvent region
  pGradientMatch_sSmooth,              // Gradient Restoration of the protein region, smoothing of solvent region  
  pGlobicize_pHistogramMatch_sFlatten  // Globicization followed by Histogram matching of the protein region, flattening of solvent region
};

struct beta_param // populated in the experiment parameter file
{
  int     iterate;            // Iterate at which the settings apply
  double  midline;            // Midline value of oscillations
  double  amplitude;          // Amplitude of oscillations
  float   phase;              // Phase offset, specified as a fraction of 2π
  float   period;             // Period of oscillations (specified in iterations)
  float   functionparam_a;    // a function specific parameter
  float   functionparam_b;    // a function specific parameter
};

struct filter_param // populated in the experiment parameter file  
{
  int     iterate;            // Iterate at which the settings apply
  double  radius;             // Filter radius for computation of the local mean or local variance (Angstroms)
};

struct solvent_fraction_param // populated in the experiment parameter file 
{
  int     iterate;            // Iterate at which the settings apply
  double  multiplier;         // Multiplier to be applied to the target solvent fraction
};

struct ipa_atom
{
  /*
  Its all very well and good using Clippers inbuilt Atom class, 
  however it has no basic constructor, except for from another "atom like" class.
  This structure should do the trick.
  */
  clipper::String element_;
  clipper::U_aniso_orth u_aniso_;
  float occ_;
  clipper::Coord_orth xyz_;

  // Simple constructor
  ipa_atom(clipper::String element, clipper::U_aniso_orth u_aniso, float occ, clipper::Coord_orth xyz) : element_(element), u_aniso_(u_aniso), occ_(occ), xyz_(xyz) {}; 
  
  // These methods allow this object to be used as an clipper::Atom constructor object.
  clipper::String& element() {return element_;}
  clipper::U_aniso_orth& u_aniso_orth() {return u_aniso_;}
  float& occupancy() {return occ_;}
  clipper::Coord_orth& coord_orth() {return xyz_;}
};

struct ipa_setting_algorithm_container
{
  // Total number of iterations in the regime
  int n_iterates;

  // List of update rules and constrants to be used, in the order performed. These are repeated until we reach the total number of iterations
  
  std::vector<IPA_updaterule_command> commands; // The update rule to apply
  std::vector<int> n_commands; // The number of times the update rule is applied

  std::vector<IPA_constraint_command> constraints; // The constraints to apply
  std::vector<int> n_constraints; //  The number of times the constraints are applied 
};

class util 
{
    public:
    static std::string return_name_from_enum(IPA_updaterule_command myalg);
    static std::string return_name_from_enum(IPA_constraint_command myalg); 

    // Helper functions to quickly retrieve or store a Sigma value from a clipper:String as is used in a mtz title.
    static bool extract_sigma_from_mtz_title(const clipper::String title,  double& sigma);
    static clipper::String encode_sigma_into_mtz_title(clipper::String old_title, const double sigma);

    static bool predict_success(const std::vector<float>& phase_velocity, int& iterate_of_success, int edge_padding = 5);

    static float calculate_new_resolution_limit_from_volume_scale_factor(const float old_resolution, const float scale_factor);
    static double calculate_scale_factor_for_half_width_at_nth_max(const int divisor);
    static void print_sigma_resolution_look_up_table(std::ostream& outlog, double start_res = 1.5, double end_res = 5.0, double increment = 0.1, double nth_height = 20);

    static void phases_to_log(clipper::HKL_data<clipper::data64::F_phi>& data,
                        std::ostream& outlog,
                        std::string heading = "Summary of phases: ",
                        int extent = 50);
    
    // These are for the Command Line only really... Parsing them to a log file will likely output the fancy ANSI hacks too... maybe.
    static std::string   parse_line(int cl, std::string line, int ansi_hacks = 0, const std::string& Contentcolour = "", const std::string& Linecolour = "");                                     // Change text to a formated Section output
    static std::ostream& print_section_line(const int size, std::ostream& outlog, const std::string& Linecolour = "");                 // Print a solid divider Section line
    static std::ostream& print_empty_line(const int size, std::ostream& outlog, const std::string& Linecolour = "");                   // Print an empty Section line.
    static std::ostream& print_progress_bar(const int size, float progress, std::ostream& outlog, const std::string& Contentcolour = "", const std::string& Linecolour = ""); // Print a progress bar, of given size and progress, formated as part of a section.
    
    // Helper function fo calculating a Debye-Waller factor for a given hkl, given an Anisotropic ADP (U*)
    static double calculate_debye_waller_factor(int h, int k, int l, std::vector<double> U_Star);
    
    static void print_all_settings_to_command_line(std::ostream& outlog, ipa_settings& settings);

    static char crystal_family(const int space_group_number);
    
    static int lattice_multiplicity(const std::string &H_M_symbol);
    
    static int enantiomorphic_space_group(const int space_group_number);
    
    static bool space_group_is_chiral(const int space_group_number);

    static void db_scan(const std::vector< std::vector<float> > &distance,
                        const double &epsilon,
                        const double &minPoints,
                        int &n_clusters,
                        std::vector<int> &cluster_id,
                        std::ostream &outlog);

  // Function for printing a comparison between experimental intensities (derived from F-SigF) and model intensities based on Wilson statistics.
  static void print_intensity_statistics_data_vs_model(const clipper::HKL_data<clipper::data64::F_sigF>& f, 
                                                        const clipper::HKL_data<clipper::data64::I_sigI>& I, 
                                                        const int res_bins, 
                                                        std::ostream& outlog,
                                                        const double res_bin_start = -1,
                                                        const double res_bin_end = -1);

  static std::vector<std::string> convert_fasta_to_globule_list(const std::string input, 
                                                      bool sorted = true);
  
  static void version_fingerprint(std::ostream& outlog);
};

#endif