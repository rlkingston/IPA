// IPA: Iterative Projection Algorithms for protein crystallography
// Utility function library 

#include <iomanip>

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

#include "util-lib.h"

#include "svn-version.h"

#include <cmath> // std::pow 

#include <iomanip> // std::setw & std::setprecison & std::boolalpha

#include <Eigen/Dense>

std::string util::return_name_from_enum(IPA_updaterule_command myalg)
{
  /*
   Simple helper function which provides
   a string of the full algorithm name
   from the the codes enum for logging
   purposes.
   */
  switch(myalg)
  {
    case IPA_updaterule_command::DifMap : return "Difference Map Algorithm";
    case IPA_updaterule_command::ErrRed : return "Error Reduction Algorithm ";
    case IPA_updaterule_command::ReReRe : return "Relax Reflect Relax Algorithm";
    case IPA_updaterule_command::RevRRR : return "Reverse Relax Reflect Relax Algorithm";
  }
  
  return "Undefined";
}

std::string util::return_name_from_enum(IPA_constraint_command myalg)
{
  /*
   Simple helper function which provides
   a string of the full algorithm name
   from the the codes enum for logging
   purposes.
   */
  switch(myalg)
  {
  case IPA_constraint_command::pHistogramMatch_sFlatten : return "Histogram Matching of the protein region, with solvent flattening";
  case IPA_constraint_command::pGradientMatch_sFlatten : return "Gradient Restoration of the protein region, with solvent flattening";
  case IPA_constraint_command::psSkeletonize : return "Skeletonization of the protein and solvent regions";
  case IPA_constraint_command::pSkeletonize_sFlatten : return "Skeletonization of the protein region, with solvent flattening";
  case IPA_constraint_command::pWeightedSkele_sFlatten : return "Skeletonization of the protein region, merged with weightings to pre-skeleton map, with solvent flattening";
  case IPA_constraint_command::pAtomize_sFlatten : return "Atomization of the protein region, with solvent flattening";
  case IPA_constraint_command::pGlobicize_sFlatten : return "Globicization of the protein region, with solvent flattening";
  case IPA_constraint_command::pHistogramMatch_sSmooth : return "Histogram matching of the protein region, with solvent smoothing";
  case IPA_constraint_command::pGradientMatch_sSmooth : return "Gradient Restoration of the protein region, with solvent smoothing";
  case IPA_constraint_command::pGlobicize_pHistogramMatch_sFlatten : return "Globicization and Histogram Matching of the protein region, with solvent flattening.";
  }
  
  return "Undefined";
}

bool util::predict_success(const std::vector<float>& data, int& iterate_of_success, int edge_padding)
{
  /*
  Expolratory function for monitoring algorithm performance
  This simple function is used to predict whether a solution has been located, based on the change in some statistic with the iterate
  The idea is simple, there is a point within the data, where the values on the left are significantly different to the values on the right.
  We locate this point by looping through the data, and calculating the average on each side of the current point
  */
  bool success = false;
  double left_average = 0;
  double right_average = 0;
  int r_count = 0;
  int l_count = 0;

  double maximum_difference = 0;
  int iterate_of_max_dist = 0;

  if (data.size() < edge_padding*2) return false; // Not enough data yet to determine success.
  if (!(edge_padding > 0)) return false; // Need SOME edge padding else this doesn't work.

  for (int i = edge_padding ; i < (data.size()-edge_padding); i++ )
  {
    for (int l = 0; l < i; l++)
    {
      l_count++;
      left_average += data[l];
    }
    for (int r = i; r < data.size(); r++)
    {
      r_count++;
      right_average += data[r];
    }

    right_average /= r_count;
    left_average /= l_count;

    double distance = left_average - right_average;

    if (distance > maximum_difference)
    {
      iterate_of_max_dist = i;
      maximum_difference = distance;
    }
  }

  // From here, we have deduced the point of maximum difference between the left-mean and the right-mean  
  // We could then fit two straight lines to the data, on either side of our point. 
  // If these lines have difference Y intercept values and similar gradients, that might mean we have found a solution.

  if (maximum_difference  > 5)
  {
    iterate_of_success = iterate_of_max_dist;
    success = true;
  }

  return success;
}

bool util::extract_sigma_from_mtz_title(const clipper::String title,  double& sigma)
{
  sigma = -1;
  int i = title.find("SA-");
  if (i != -1)
  {
    clipper::String sub = title.substr(i+3,title.size() - i+2); // Sigma is always appended to end.
    sigma = sub.f64();
    return true;
  }
  return false; // No SA- keyword exists, therefore there is no sigma to extract.
}

clipper::String util::encode_sigma_into_mtz_title(clipper::String old_title, const double sigma)
{
  clipper::String new_str = old_title;
  new_str.append("SA-" + std::to_string(sigma));
  return new_str;
}

void util::print_sigma_resolution_look_up_table(std::ostream& outlog, double start_res, double end_res, double increment, double nth_height)
{
    outlog << "Sigma versus Half-width of Gaussian at "<< std::setprecision(0) << std::fixed << int(nth_height) << "th the maximum height" << std::endl;
    double half_width_scale = std::sqrt( 2.0 * std::log(nth_height) );
    outlog << "*--------------*---------*" << std::endl;
    outlog << "| Sigma (1/Å)  | Res.(Å) |" << std::endl;
    outlog << "*--------------*---------*" << std::endl;
    for (double r = start_res; r <= end_res; r += increment)
    {
      double sigma = 1.0/(half_width_scale * r);
      outlog << "| " << std::fixed << std::setw(12) << std::setprecision(8) << sigma << " |   " << std::setw(5) << std::setprecision(2) << std::fixed << r << " |" << std::endl;
    }
    outlog << "*--------------*---------*" << std::endl;
    return;
}

double util::calculate_scale_factor_for_half_width_at_nth_max(const int divisor)
{
  /*
   Calculates the multiplier to compute the half width at nth of maximum of a Gaussian function, from the standard deviation 
 
   Used for specification and reporting of the apodization function.
  
   The divisor n indicates the fraction of the maximum at which the Gaussian half width is reported
  
   i.e. with the divisor = 20, the function returns the multiplier for calculating the  "Half width at a Twentieth of Maximum" which is 2.4477
   
   / N.B. For a Gaussian distribution, with standard deviation σ
   // Full Width at Half Maximum (FWHM)      = 2*SQRT(2*ln[ 2])*σ  = 2.3548*σ, Hence Half Width at Half Maximum (HWHM)     = 1.1774*σ ...
   // Full Width at 1/10th of Maximum (FWTM) = 2*SQRT(2*ln[10])*σ  = 4.2919*σ, Hence Half Width at Tenth of Maximum (HWTM) = 2.1460*σ ...
   // Full Width at 1/20th of Maximum        = 2*SQRT(2*ln[20])*σ  = 4.8955*σ, Hence Half Width at Twentieth of Maximum    = 2.4477*σ ...
   // Full Width at 1/100th of Maximum       = 2*SQRT(2*ln[100])*σ = 6.0697*σ, Hence Half Width at Hundredth of Maximum    = 3.0349*σ
   */
  return std::sqrt( 2.0 * std::log(divisor) );
}

float util::calculate_new_resolution_limit_from_volume_scale_factor(const float old_resolution_limit, const float scale_factor)
{
  /*
   Returns the resolution limit in Fourier space, that results from increasing the volume of the sphere enclosing the data by the supplied scale factor 
   Volume of the sphere in Fourier space = 4/3 π s ^ 3 where s is the magnitude of the scattering vector = 1/resolution
   */
  float old_s_limit = 1/old_resolution_limit; 
  float old_volume = (4.0/3.0) * clipper::Util::pi() * std::pow(old_s_limit,3);
  float new_volume = old_volume * scale_factor;
  return 1.0/std::cbrt( 3.0 * ( new_volume / ( 4.0 *  clipper::Util::pi() ) ));
}

void util::phases_to_log(clipper::HKL_data<clipper::data64::F_phi>& data, std::ostream& outlog, std::string heading, int extent)
{
  // DEBUG quick Print out some of the known phases to ensure we have the right ones.
  outlog << heading << std::endl;
  clipper::HKL_data_base::HKL_reference_index ih;
  outlog << " Index    h    k    l          F           Phi" << std::endl;
  for ( ih = data.first(); ih.index() <= extent; ih.next() )
  {
    outlog << std::setprecision(8) <<
    std::setw(5) << ih.index()   << " " <<
    std::setw(5) << ih.hkl().h() << " " <<
    std::setw(5) << ih.hkl().k() << " " <<
    std::setw(5) << ih.hkl().l() << " " <<
    std::setw(20) << data[ih].f() <<  " "  <<
    std::setw(10) << data[ih].phi() << " " <<
    std::fixed << std::endl;
  }
  outlog << "*---------*" << std::endl;
}

std::ostream& util::print_section_line(const int size, std::ostream& outlog, const std::string& Linecolour)
{
  // Helper function to draw an section line for update screen.
  outlog << Linecolour;
  outlog << "◉";
  for (int i = 0; i < size; i++)
  {outlog << "━";}
  outlog << "◉\n";
  if (Linecolour != "") outlog << cReset; // Only resetting when needed, to clean actualy .txt logs.
  return outlog;
}

std::ostream& util::print_empty_line(const int size, std::ostream& outlog, const std::string& Linecolour)
{
  // Helper function to draw an empty output line for update screen.
  outlog << Linecolour;

  outlog << "┃";
  for (int i = 0; i < size; i++)
  {outlog << " ";}
  outlog << "┃\n";

  if (Linecolour != "") outlog << cReset; // Only resetting when needed, to clean actualy .txt logs.
  return outlog;
}

std::ostream& util::print_progress_bar(const int size, float progress, std::ostream& outlog, const std::string& Contentcolour, const std::string& Linecolour)
{
  // Helper function to draw a progress bar with a given progress for update screen.
  std::string char_array[] = {"▏","▎","▍","▌","▋","▊","▉","█"};
  outlog << Linecolour;
  outlog << "┃   ◁";
  outlog << Contentcolour;
  int barsize = size - 12;
  float pos;// = barsize * progress;
  float frac = std::modf((barsize * progress), &pos); // Get the integeral, and fractional part of progress.
  int array_pos = frac * 8; // Which character to get based on fractional progress.
  for (int i = 0; i < barsize; i++)
  {
    if (i < pos) outlog << "█"; // Full Load Bar
    else if (i == pos) outlog << char_array[array_pos]; // Fractional Load Bar
    else outlog << " "; // Empty Load Bar
  }
  
  char buffer[64];
  sprintf(buffer,"%3.0f", progress*100);
  outlog << Linecolour;
  outlog << "▷ " << buffer << "%  ┃\n";
  outlog << cReset;
  return outlog;
}

std::string util::parse_line(int cl, std::string line, int ansi_hacks, const std::string& Contentcolour, const std::string& Linecolour)
{
  // Helper function to format text into a format that will adhere to our output style.
  // This function takes a string of any length, and breaks it at the nearest previous space, and pads the output with the formatting lines.
  // The problem with the ANSI hacks offset, is that it can only be applied to a single line, and there's no a priori way to know which line it'll end up in.
  // So the Ansi hack works with single lines, but not multiple lines
  
  std::string new_output = "";
  
  while (line.length() > (cl-7))
  {
    int cut_point = cl-7; // Define the default cut point
    
    for (int i = cut_point; i > 0; i --) // Find the next best cut point
    {
      if (line[i] == ' ' || line[i] == '/') // Lets also cut at these slashes, as sometimes theres no spaces.
      {
        cut_point = i;
        break;
      }
    }
    
    std::string line_sub = line.substr(0,cut_point);
    line_sub.resize(cl-7+ansi_hacks,' '); // Pad to correct size. and allow ANSI hacks
    new_output += Linecolour + "┃    " + Contentcolour + line_sub + Linecolour + "   ┃\n";
    // Alter line by removing what we've already added to output.
    line = line.substr(cut_point,line.length()); // Truncate what we've added so far.
  }
  
  // Print whats left.
  if (line.length() < cl-7)
  {   // Extend that line-O
    line.resize(cl-7+ansi_hacks,' '); // Pad to correct size, and allow for ANSI hacks
    new_output += Linecolour + "┃    " + Contentcolour + line + Linecolour + "   ┃\n";
  }
  if (Contentcolour != "" || Linecolour != "") new_output += cReset; // Only resetting when needed, to clean actualy .txt logs.
  return new_output;
}

double util::calculate_debye_waller_factor(int h, int k, int l, std::vector<double> U_Star)
{
  /*
   This function takes in the symmetric second-rank tensor U*, and Miller indices hkl, and returns the Debye-Waller Factor, as described in Grosse-Kunstieve and Adams (2002).
   
   Note the U_Star 3x3 matrix is submitted as a 1-dimensional vector of size 9. Indices 0 through to 8.
   
   The function being evaluated is 
  
   T(h) = exp ( -2PI^2 ( ht U* h ) )
  
   where
   T is the Debye-Waller Factor, as a function of h.
   h is the column vector containing indices h, k, and l
   ht is the transpose of h, i.e. a row vector containing indices h, k, and l.
   U* is the dimensionless and symmetric 3x3 tensor describing the atomic displacements
      
   Returns -1 if the calculation fails.
   
   */
  if (U_Star.size() != 9)
  {
    std::cerr << " Input incorrectly sized U_Star mean-sqaured displacement parameter, requires a 9 sized vector, representing a 3x3 matrix of values..." << std::endl;
    throw std::invalid_argument( "Received incorrectly sized mean-squared displacement U* parameter" );
    return -1;
  }

  //if (U_Star[1] != U_Star[3]) {std::cout << "U* Symmetry check failed, 1 and 3." << std::endl;}
  //if (U_Star[2] != U_Star[6]) {std::cout << "U* Symmetry check failed, 2 and 6." << std::endl;}
  //if (U_Star[5] != U_Star[7]) {std::cout << "U* Symmetry check failed, 5 and 7." << std::endl;}
  
  // TODO: Double check the ordering is right here... However whether it is row ordered or column ordered will not matter in this particular context
  // as U* is symmetric around the diagonal
  
  //std::cout << std::endl;
  //std::cout << U_Star[0] << " " << U_Star[1] << " " << U_Star[2] << std::endl;
  //std::cout << U_Star[3] << " " << U_Star[4] << " " << U_Star[5] << std::endl;
  //std::cout << U_Star[6] << " " << U_Star[7] << " " << U_Star[8] << std::endl;
  //std::cout << std::endl;
  
  // First multiply U* and h, to generate the components of the 3x1 matrix [uh uk ul]t
  
  double uh = U_Star[0]*h + U_Star[1]*k + U_Star[2]*l;
  double uk = U_Star[3]*h + U_Star[4]*k + U_Star[5]*l;
  double ul = U_Star[6]*h + U_Star[7]*k + U_Star[8]*l;
  
  double negativetwopisquared = -2.0 * (clipper::Util::pi() * clipper::Util::pi());
  
  return std::exp(negativetwopisquared * ((uh * h) + (uk * k) + (ul * l)));
}

void util::print_all_settings_to_command_line(std::ostream& outlog, ipa_settings& settings)
{
  for (int i = 0; i < settings.beta_params_look_up_table.size(); i ++)
  {
    outlog << "Beta Param [" << i << "]: ";
    outlog << settings.beta_params_look_up_table[i].iterate << ", ";
    outlog << settings.beta_params_look_up_table[i].midline << ", ";
    outlog << settings.beta_params_look_up_table[i].amplitude << ", ";
    outlog << settings.beta_params_look_up_table[i].phase << ", ";
    outlog << settings.beta_params_look_up_table[i].period << ", ";
    outlog << settings.beta_params_look_up_table[i].functionparam_a << ", ";
    outlog << settings.beta_params_look_up_table[i].functionparam_b << std::endl;
  }
  
  for (int i = 0; i < settings.solvent_fraction_multiplier_params_look_up_table.size(); i ++)
  {
    outlog << "solvent multiplier Param [" << i << "]: ";
    outlog << settings.solvent_fraction_multiplier_params_look_up_table[i].iterate << ", ";
    outlog << settings.solvent_fraction_multiplier_params_look_up_table[i].multiplier << std::endl;
  }
  
  for (int i = 0; i < settings.filter_radius_params_look_up_table.size(); i ++)
  {
    outlog << "filter radius Param [" << i << "]: ";
    outlog << settings.filter_radius_params_look_up_table[i].iterate << ", ";
    outlog << settings.filter_radius_params_look_up_table[i].radius << std::endl;
  }
  
  outlog << "n_total_iterates: ";
  outlog << settings.n_total_iterates << std::endl;
  outlog << "input_filename_target_data: ";
  outlog << settings.input_filename_target_data << std::endl;
  outlog << "target_col_fo: ";
  outlog << settings.target_col_fo << std::endl;
  outlog << "target_col_pw: ";
  outlog << settings.target_col_pw << std::endl;
  outlog << "target_col_pw_calc: ";
  outlog << settings.target_col_pw_calc << std::endl;
  outlog << "output_phases: ";
  outlog << settings.output_phases << std::endl;
  outlog << "phase_output_interval: ";
  outlog << settings.phase_output_interval << std::endl;
  outlog << "output_envelopes_and_phases: ";
  outlog << settings.output_envelopes_and_phases << std::endl;
  outlog << "input_known_phases: ";
  outlog << settings.input_known_phases << std::endl;
  outlog << "input_working_mask: ";
  outlog << settings.input_working_mask << std::endl;
  outlog << "input_filename_target_mask_working ";
  outlog << settings.input_filename_target_mask_working << std::endl;
  outlog << "input_known_mask: ";
  outlog << settings.input_known_mask << std::endl;
  outlog << "input_filename_known_mask: ";
  outlog << settings.input_filename_known_mask << std::endl;
  //outlog << "fix_mask: ";
  //outlog << settings.fix_mask << std::endl;
  outlog << "n_iterations_fix_mask: ";
  outlog << settings.n_iterations_fix_mask << std::endl;
  outlog << "compute_mask_from_local_mean: ";
  outlog << settings.compute_mask_from_local_mean << std::endl;
  outlog << "compute_mask_from_local_variance: ";
  outlog << settings.compute_mask_from_local_variance << std::endl;
  outlog << "compute_mask_from_local_mean_and_variance: ";
  outlog << settings.compute_mask_from_local_mean_and_variance << std::endl;
  outlog << "compute_mask_interval: ";
  outlog << settings.compute_mask_interval << std::endl;
  outlog << "mask_weighting_factor_mean: ";
  outlog << settings.mask_weighting_factor_mean << std::endl;
  outlog << "filter_radius: ";
  outlog << settings.filter_radius << std::endl;
  //outlog << "percentage_radius_expansion: ";
  //outlog << settings.percentage_radius_expansion << std::endl;
  outlog << "output_mask: ";
  outlog << settings.output_mask << std::endl;
  outlog << "mask_output_interval: ";
  outlog << settings.mask_output_interval << std::endl;
  outlog << "found_filter_params: ";
  outlog << settings.found_filter_params << std::endl;
  outlog << "found_beta_params: ";
  outlog << settings.found_beta_params << std::endl;
  outlog << "found_solvent_fraction_multiplier_params: ";
  outlog << settings.found_solvent_fraction_multiplier_params << std::endl;
  outlog << "erase_islands: ";
  outlog << settings.erase_islands << std::endl;
  outlog << "erase_islands_interval: ";
  outlog << settings.erase_islands_interval << std::endl;
  outlog << "erase_islands_threshold: ";
  outlog << settings.erase_islands_threshold << std::endl;
  outlog << "fill_voids: ";
  outlog << settings.fill_voids << std::endl;
  outlog << "fill_voids_interval: ";
  outlog << settings.fill_voids_interval << std::endl;
  outlog << "fill_voids_threshold: ";
  outlog << settings.fill_voids_threshold << std::endl;
  outlog << "envelope_editing_limit: ";
  outlog << settings.envelope_editing_limit << std::endl;
  outlog << "solvent_flattening: ";
  outlog << settings.solvent_flattening << std::endl;
  outlog << "histogram_matching: ";
  outlog << settings.histogram_matching << std::endl;
  outlog << "apodize_reference_data: ";
  outlog << settings.apodize_reference_data << std::endl;
  outlog << "histogram_matching_interval: ";
  outlog << settings.histogram_matching_interval << std::endl;
  outlog << "gradient_matching: ";
  outlog << settings.gradient_matching << std::endl;
  outlog << "gradient_matching_interval: ";
  outlog << settings.gradient_matching_interval << std::endl;
  outlog << "input_starting_phases: ";
  outlog << settings.input_starting_phases << std::endl;
  outlog << "use_test_set: ";
  outlog << settings.use_test_set << std::endl;
  outlog << "n_test_sets: ";
  outlog << settings.n_test_sets << std::endl;
  outlog << "target_bfactor: ";
  outlog << settings.target_b_factor_data.B_ISO << std::endl;
  outlog << "fraction_solvent: ";
  outlog << settings.fraction_solvent << std::endl;
  outlog << "expected_mean_density_protein: ";
  outlog << settings.expected_mean_density_protein << std::endl;
  outlog << "expected_mean_density_solvent: ";
  outlog << settings.expected_mean_density_solvent << std::endl;
  outlog << "hemihedral_twinning: ";
  outlog << settings.hemihedral_twinning << std::endl;
  outlog << "hemihedral_twin_case: ";
  outlog << settings.hemihedral_twin_case << std::endl;
  outlog << "hemihedral_twin_fraction: ";
  outlog << settings.hemihedral_twin_fraction << std::endl;
  outlog << "apodize: ";
  outlog << settings.apodize << std::endl;
  outlog << "n_apodization_steps: ";
  outlog << settings.n_apodization_steps << std::endl;
  outlog << "n_apodization_first_iterate: ";
  outlog << settings.n_apodization_first_iterate << std::endl;
  outlog << "n_apodization_last_iterate: ";
  outlog << settings.n_apodization_last_iterate << std::endl;
  outlog << "initial_apodization_sigma: ";
  outlog << settings.initial_apodization_sigma << std::endl;
  outlog << "final_apodization_sigma: ";
  outlog << settings.final_apodization_sigma << std::endl;
  outlog << "upper_apodization_limit: ";
  outlog << settings.upper_apodization_limit << std::endl;
  outlog << "linear_freq_change: ";
  outlog << settings.linear_freq_change << std::endl;
  //outlog << "fix_beta: ";
  //outlog << settings.fix_beta << std::endl;
  outlog << "beta_function: ";
  outlog << settings.beta_function << std::endl;
  outlog << "compute_real_space_agreement_statistics: ";
  outlog << settings.compute_real_space_agreement_statistics << std::endl;
  //outlog << "algorithm_dm: ";
  //outlog << settings.algorithm_dm << std::endl;
  //outlog << "algorithm_rrr: ";
  //outlog << settings.algorithm_rrr << std::endl;
  outlog << "n_runs: ";
  outlog << settings.n_runs << std::endl;
  //outlog << "n_iterations_ipa_first: ";
  //outlog << settings.n_iterations_ipa_first << std::endl;
  //outlog << "n_iterations_ipa_regular: ";
  //outlog << settings.n_iterations_ipa_regular << std::endl;
  //outlog << "n_iterations_ipa_last: ";
  //outlog << settings.n_iterations_ipa_last << std::endl;
  //outlog << "n_iterations_error_reduction_first: ";
  //outlog << settings.n_iterations_error_reduction_first << std::endl;
  //outlog << "n_iterations_error_reduction_regular: ";
  //outlog << settings.n_iterations_error_reduction_regular << std::endl;
  //outlog << "n_iterations_error_reduction_last: ";
  //outlog << settings.n_iterations_error_reduction_last << std::endl;
  //outlog << "n_iterations_shrink: ";
  //outlog << settings.n_iterations_shrink << std::endl;
  outlog << "really_random: ";
  outlog << settings.really_random << std::endl;
  outlog << "fourier_space_phase_generation: ";
  outlog << settings.fourier_space_phase_generation << std::endl;
  outlog << "n_sphere_phase_generation: ";
  outlog << settings.n_sphere_phase_generation << std::endl;
  outlog << "rad_sphere_phase_generation: ";
  outlog << settings.rad_sphere_phase_generation << std::endl;
  outlog << "max_attempts: ";
  outlog << settings.max_attempts << std::endl;
  outlog << "remove_low_res_data: ";
  outlog << settings.remove_low_res_data << std::endl;
  outlog << "low_res_cutoff: ";
  outlog << settings.low_res_cutoff << std::endl;
  outlog << "p_threshold: ";
  outlog << settings.p_threshold << std::endl;
  outlog << "replace_largeI_with_square_of_expectedF: ";
  outlog << settings.replace_largeI_with_square_of_expectedF << std::endl; 
  outlog << "intensity_multiplier: ";
  outlog << settings.intensity_multiplier << std::endl;
  outlog << "res_bins: ";
  outlog << settings.res_bins << std::endl;
  outlog << "verbose: ";
  outlog << settings.verbose << std::endl;
  outlog << "n_density_bins: ";
  outlog << settings.n_density_bins << std::endl;
  outlog << "n_gradient_bins: ";
  outlog << settings.n_gradient_bins << std::endl;  
  outlog << "shannon_rate: ";
  outlog << settings.shannon_rate << std::endl;
  /*
   outlog << ": ";
   outlog << settings. << std::endl;
   */
}

int util::lattice_multiplicity(const std::string &H_M_symbol)
{
  // Simple minded function to return the lattice multiplicity given the extended Hermann-Mauguin symbol for the space group
  // The lattice multiplcity (mL) is important for calculating the enhancement in the average intensity of the observations arising when a centered cell is employed
  // If there is a systematic extinction of a fraction [1 - (1/mL )] of the Fourier amplitudes due to lattice centering, the total scattering must be concentrated in the allowed fraction 1/mL of the Fourier amplitudes. 
  // For discussion see Blessing RH, Guo DY, Langs DA. Intensity Statistics and Normalization. In: Fortier S, editor. Direct Methods for Solving Macromolecular Structures. Kluwer Academic Publishers; 1998. pp. 47–71. 

  int lattice_multiplicity = 1; // Function will return 1 as default, this is the correct multiplicity for a Primitive lattice, also the Rhombohedral setting of the Rhombohedral space groups
   
  if      ( (H_M_symbol.front() == 'C') || (H_M_symbol.front() == 'B') || (H_M_symbol.front() == 'A') || (H_M_symbol.front() == 'I') ) {lattice_multiplicity = 2;}
  else if   (H_M_symbol.front() == 'F')                                                                                                {lattice_multiplicity = 4;}
  else if  ((H_M_symbol.front() == 'R')  && (H_M_symbol.back() == 'H') )                                                               {lattice_multiplicity = 3;} // Applicable to the hexagonal setting of the rhombohedral space groups which will have extended symbols R 3 :H and R 3 2 :H
  
  return lattice_multiplicity;
  
}

char util::crystal_family(const int space_group_number)
  
{
  // a simple function that returns the standard symbol for the crystal family, given the space group number as an argument
  // Symbols for the families are a (anorthic = triclinic), m (monoclinic), o (orthorhombic), t (tetragonal), h (hexagonal) and c (cubic).
  // see Nespolo M, Aroyo MI, Souvignier B. Crystallographic shelves: space‐group hierarchy explained. J Appl Cryst.2018 Oct 1;51(5):1481–91. 
   
  char crystal_family;

  if      ( (space_group_number >= 1) && (space_group_number <= 2) )
  {
    crystal_family = 'a'; // anorthic = triclinic
  }
  else if ( (space_group_number >= 3) && (space_group_number <= 15) )
  {
    crystal_family = 'm'; // monoclinic
  }
  else if ( (space_group_number >= 16) && (space_group_number <= 74) )
  {
    crystal_family = 'o'; // orthorhombic
  }
  else if ( (space_group_number >= 75) && (space_group_number <= 142) )
  {
     crystal_family = 't'; // tetragonal
  }
  else if ( (space_group_number >= 143) && (space_group_number <= 194) )
  {
    crystal_family = 'h'; // hexagonal
  }
  else if ( (space_group_number >= 195) && (space_group_number <= 230) )
  {
    crystal_family = 'c'; // cubic
  }
  else 
  {
    std::cerr << "\n Illegal space group number supplied to function crystal_family \n"  << std::endl; 
  }

  
  return crystal_family;
  
}

int util::enantiomorphic_space_group(const int space_group_number)
{
  // Simple minded function to delivers the enantiomorph of a chiral space group (with space groups identified by number)
  // Only works for the 65 space groups relevant to protein crystallography


  int enantiomorphic_space_group;
    
  switch (space_group_number){

    //---------------------------------
    case 76 : // space group P41
      enantiomorphic_space_group = 78;
      break;
      
    case 78 : // space group P43
      enantiomorphic_space_group = 76;
      break;

    //---------------------------------    
    case 91 : // space group P4122
      enantiomorphic_space_group = 95;
      break;
    
    case 95 : // space group P4322
      enantiomorphic_space_group = 91;
      break;

    //--------------------------------- 
    case 92 : // space group P41212
      enantiomorphic_space_group = 96;
      break;
 
    case 96 : // space group P43212
      enantiomorphic_space_group = 92;
      break;

    //--------------------------------- 
    case 144 : // space group P31
      enantiomorphic_space_group = 145;
      break;
 
    case 145 : // space group P32
      enantiomorphic_space_group = 144;
      break;
 
     //---------------------------------
    case 151 : // space group  P3112
      enantiomorphic_space_group = 153;
      break;
 
    case 153 : // space group  P3212
      enantiomorphic_space_group = 151;
      break;
 
     //---------------------------------
    case 152 : // space group  P3121
      enantiomorphic_space_group = 154;
      break;
 
    case 154 : // space group  P3221
      enantiomorphic_space_group = 152;
      break;

    //--------------------------------- 
    case 169 : // space group P61
      enantiomorphic_space_group = 170;
      break;
 
    case 170 : // space group P65
      enantiomorphic_space_group = 169;
      break;
      
     //---------------------------------
    case 171 : // space group P62
      enantiomorphic_space_group = 172;
      break;
 
    case 172 : // space group P64
      enantiomorphic_space_group = 171;
      break;
 
     //---------------------------------
    case 178 : // space group P6122
      enantiomorphic_space_group = 179;
      break;
 
    case 179 : // space group P6522
      enantiomorphic_space_group = 178;
      break;

    //--------------------------------- 
    case 180 : // space group P6222
      enantiomorphic_space_group = 181;
      break;
 
    case 181 : // space group P6422
      enantiomorphic_space_group = 180;
      break;
      
    //--------------------------------- 
    case 212 : // space group P4332
      enantiomorphic_space_group = 213;
      break;
 
    case 213 : // space group P4132
      enantiomorphic_space_group = 212;
      break;

    //---------------------------------
           
    default : 
      enantiomorphic_space_group = space_group_number;
      
    }
    
  return enantiomorphic_space_group;
  
}

bool util::space_group_is_chiral(const int space_group_number)
{
  // In theory clipper spacegroup().invariant_under_change_of_hand() should take care of this, but it seems to have some issues.
  // So we do this the simple-minded way - the function returns true if space group is chiral and false if it is not, based on space group number
  
  bool space_group_is_chiral;
  
  switch (space_group_number){
      
      // the 43 achiral space groups relevant to protein crystallography
      
    case 1 : // space group P1
    case 3 : // space group P2
    case 4 : // space group P21
    case 5 : // space group C2
    case 16 : // space group P222
    case 17 : // space group P2122 / P2212 / P2221
    case 18 : // space group P22121 / P21221 / P21212
    case 19 : // space group P212121
    case 20 : // space group C2221
    case 21 : // space group C222
    case 22 : // space group F222
    case 23 : // space group I222
    case 24 : // space group I212121
    case 75 : // space group P4
    case 77 : // space group P42
    case 79 : // space group I4
    case 80 : // space group I41
    case 89 : // space group P422
    case 90 : // space group P4212
    case 93 : // space group P4222
    case 94 : // space group P42212
    case 97 : // space group I422
    case 98 : // space group I4122
    case 143 : // space group P3
    case 146 : // space group R3/H3
    case 149 : // space group  P312
    case 150 : // space group  P321
    case 155 : // space group  R32/H32
    case 168 : // space group P6
    case 173 : // space group P63
    case 177 : // space group P622
    case 182 : // space group P6322
    case 195 : // space group P23
    case 198 : // space group P213
    case 196 : // space group F23
    case 197 : // space group I23
    case 199 : // space group I213
    case 207 : // space group P432
    case 208 : // space group P4232
    case 209 : // space group F432
    case 210 : // space group F4132
    case 211 : // space group I432
    case 214 : // space group I4132
      
      space_group_is_chiral = false;
      
      break;
      
      // the 22 chiral space groups arranged by enantiomorphic pair
      
    case 76 : // space group P41
    case 78 : // space group P43
    case 91 : // space group P4122
    case 95 : // space group P4322
    case 92 : // space group P41212
    case 96 : // space group P43212
    case 144 : // space group P31
    case 145 : // space group P32
    case 151 : // space group  P3112
    case 153 : // space group  P3212
    case 152 : // space group  P3121
    case 154 : // space group  P3221
    case 169 : // space group P61
    case 170 : // space group P65
    case 171 : // space group P62
    case 172 : // space group P64
    case 178 : // space group P6122
    case 179 : // space group P6522
    case 180 : // space group P6222
    case 181 : // space group P6422
    case 212 : // space group P4332
    case 213 : // space group P4132
      
      space_group_is_chiral = true;
      
      break;
      
    default :
      std::cerr << "\n Uh-Oh - Looks like the chirality of your space group isn't encoded yet \n"  << std::endl;
      // avert your eyes c++ aficionados
  }
  
  return space_group_is_chiral;
  
}

void util::db_scan(const std::vector< std::vector<float> > &distance, const double &epsilon, const double &minPoints, int &n_clusters, std::vector<int> &cluster_id, std::ostream &outlog)
{
  // Clustering using algorithm DB-SCAN
  // Input is the pairwise distances of the elements to be clustered & DB-SCAN algorithm parameters epsilon and minPoints
  // Output is n_clusters and cluster_id, the latter allowing lookup of the cluster assignments of each element
  
  // Within cluster_id
  // An entry of -1 indicates the element is unclassified
  // An entry of -2 indicates the element is classified as a noise point by the clustering algorithm
  // Entries 0, 1, ... ,n_clusters-1 index the clusters
  
  // See https://christianmaxmike.github.io/mindnotes/ml-clustering-dbscan
  //     https://en.wikipedia.org/wiki/DBSCAN
  
  int n_elements = distance[0].size();
  
  cluster_id.clear();
  cluster_id.resize(n_elements, -1); // initialize vector of size n_elements, with integer entries -1, to allow rapid lookup of the cluster assignment
  
  int cluster_counter = 0;
  
  // Loop over all elements
  for(int i = 0; i < n_elements; i++)
  {
    outlog << "Processing Element: " << i << std::endl;
    
    if ( cluster_id[i] == -1 ) // only execute if this element is currently unclassified
    {
      std::vector<int> neighbors; // create a integer valued vector to hold the neighbors
      for(int j = 0; j < n_elements; j++) // find the neighbors - this list *excludes* the query point istelf
      {
        if ( (i != j) && (distance[i][j] <= epsilon) )
          neighbors.push_back(j); // Add this point to the result, by augmenting the vector containing the neighbors
      }
      
      if ( (neighbors.size() + 1)  < minPoints)  // Fails the core point condition
      {
        cluster_id[i] = -2;                   // Label this element as a noise point
      }
      else                                       // Satisfies the Core point condition, so now we expand the cluster
      {
        
        cluster_id[i] = cluster_counter; // Assign the current element to this cluster
        
        std::vector<int> seed_set;
        seed_set = neighbors;
        
        for( std::vector<int>::size_type CountSeeds = 0, n = seed_set.size(); CountSeeds < n; ++CountSeeds ) // Loop over all seed points - note that the bounds of this loop will be subsequently adjusted if the seed set gets augmented in size ...
        {
          
          if (cluster_id[seed_set[CountSeeds]] == -2) cluster_id[seed_set[CountSeeds]] = cluster_counter; // Change status - Noise point becomes a border point
          
          if (cluster_id[seed_set[CountSeeds]] == -1) // The point is unclassified so deal with it
          {
            cluster_id[seed_set[CountSeeds]] = cluster_counter; // Label the point
            
            neighbors.clear();
            for(int j = 0; j < n_elements; j++) // find the neighbors - this list now *includes* the query point
            {
              if (distance[seed_set[CountSeeds]][j] <= epsilon)
                neighbors.push_back(j); // Add this point to the result, by augmenting the vector containing the neighbors
            }
            
            if (neighbors.size() >= minPoints) // Satisfies the Core point condition
            {
              seed_set.insert( seed_set.end(), neighbors.begin(), neighbors.end() ); // add new neighbors to the seed set
            }
            n = seed_set.size();
            
          }
        }
        cluster_counter += 1; // We have finished processing a core point, so augment the cluster counter
        
      }
    }
  } // end of loop over all elements
  
  n_clusters = cluster_counter;
}

void util::print_intensity_statistics_data_vs_model(const clipper::HKL_data<clipper::data64::F_sigF>& f, const clipper::HKL_data<clipper::data64::I_sigI>& I, const int res_bins, std::ostream& outlog, const double res_bin_start, const double res_bin_end)
{
  // Note that Wilson model now specifies < I > and not < I/ε > ... this function has been modified accordingly
   
  bool resolution_range_override = false;
  if (res_bin_start > 0 && res_bin_end > 0)
  {
    resolution_range_override = true;
  }
  // We initialise to + 2, to store the totals in the final +1 th bin. (length is already +1 for bin ID)
  // Containers for information of the histogram outputs.
  std::vector<double> data_i_bins(res_bins+2);
  std::vector<int> data_i_counts_bins(res_bins+2);
  std::vector<double> model_i_bins(res_bins+2);
  std::vector<int> model_i_counts_bins(res_bins+2);

  // Simply used to calculate and display the completeness of these datasets.
  int low_res_counts_model = 0;
  int low_res_counts_data = 0;
  int high_res_counts_model = 0;
  int high_res_counts_data = 0;

  double max_value; // To scale both Histograms to an equal playing field later...

  clipper::HKL_info hkl_data =  f.base_hkl_info();
  clipper::HKL_info hkl_model =  I.base_hkl_info();
  clipper::Resolution_ordinal resord_data; // Resolution Ordinal - once initiated this returns the approximate fractional ordinal within a dataset (in the range 0...1) for any specified value of 1 / s^2
  clipper::Resolution_ordinal resord_model; // Resolution Ordinal - once initiated this returns the approximate fractional ordinal within a dataset (in the range 0...X) for any specified value of 1 / s^2

  resord_model.init(hkl_model, 1.0);

  // Perform binning.
  clipper::HKL_data_base::HKL_reference_index ih;
  for ( ih = I.first(); !ih.last(); ih.next() )
  {
    int dh = hkl_data.index_of(ih.hkl()); // The index of the same hkl for the f data.
    int bin;
    // Calculate our bins
    bin = std::min ( int(static_cast<double>(res_bins) * resord_model.ordinal(ih.invresolsq())) , res_bins-1 );

    if (resolution_range_override)
    {
      double bin_total = res_bin_start - res_bin_end;
      bin = int(static_cast<double>(res_bins) * ( (res_bin_start - (1.0 / sqrt(ih.invresolsq())) ) / bin_total));
      if (bin > res_bins-1 || bin < 0)
      {
        continue; // if this is an out of range resolution reflection, ignore it.
      }
    }

    // Update the low vs high super bins
    if ((1.0 / sqrt(ih.invresolsq())) > 10.0 && (1.0 / sqrt(ih.invresolsq())) <= 25.0)
    {
      low_res_counts_model++;
      if (dh != -1 && !f[dh].missing())
      {
        low_res_counts_data++;
      }
    }
    else
    {
      high_res_counts_model++;
      if (dh != -1 && !f[dh].missing())
      {
        high_res_counts_data++;
      }
    }

    // We use to only compare when present in both. However we want to ensure the Wilson model observations are ALWAYS counted, as they may extend to lower resolution than the experimental data.
    model_i_bins[bin] += I[ih].I();
    model_i_bins[res_bins] += I[ih].I();
    model_i_counts_bins[bin] += 1;
    model_i_counts_bins[res_bins] += 1;

    if (dh != -1 && !f[dh].missing())
    {
      // Here are our measured_f's.
      data_i_bins[bin] += pow(f[dh].f(), 2);
      data_i_bins[res_bins] += pow(f[dh].f(), 2);
      data_i_counts_bins[bin] += 1;
      data_i_counts_bins[res_bins] += 1;
    }
  }
  
  // Compute the mean values from the raw sums
  
  for ( int bin = 0; bin < (res_bins+1); bin++ )
  {
    if (data_i_counts_bins[bin] > 0)
    {
      data_i_bins[bin] =  data_i_bins[bin] / static_cast<double>(data_i_counts_bins[bin]);
    }
    if (model_i_counts_bins[bin] > 0)
    {
      model_i_bins[bin] = model_i_bins[bin] / static_cast<double>(model_i_counts_bins[bin]);
    }
  }
  
  // Calculate the super bins completeness.
  double low_res_completeness = (double)low_res_counts_data / (double)low_res_counts_model;
  double high_res_completeness = (double)high_res_counts_data / (double)high_res_counts_model;

  double lower_limit;
  double upper_limit;
  double reso_min = 1.0/sqrt(f.invresolsq_range().min());


  resord_model.invert();

  outlog << "●━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━●" << std::endl;
  outlog << "┃              INTENSITY STATISTICS <I> : INPUT DATA vs WILSON MODEL                   ┃" << std::endl;
  outlog << "●━━━━━┯━━━━━━━━━━━━━━┯━━━━━━━━┯━━━━━━━━━━━━━━━┯━━━━━━━━┯━━━━━━━━━━━━━━━┯━━━━━━━━━━━━━━━●" << std::endl;

  outlog << "┃ Bin │ Resolution Å │      # │          Data │      # │         Model │  Data / Model ┃" << std::endl;
  outlog << "┠─────┼──────────────┼────────┼───────────────┼────────┼───────────────┼───────────────┨" << std::endl;

  double bincrement = (res_bin_start - res_bin_end)/(double)res_bins;
  for (int bin = 0 ; bin < (res_bins + 1) ; bin++ )
  {
  if (bin == 0)
  {
    if (!resolution_range_override)
    {lower_limit = std::numeric_limits<double>::infinity();}
    else 
    { lower_limit = res_bin_start;} //float().infinity //25.0; // Hard coded, because infinity.
    upper_limit = 1.0/sqrt(resord_model.ordinal( static_cast<double>(bin+1)/static_cast<double>(res_bins) ) );
    if (resolution_range_override)
    {
      lower_limit = res_bin_start;
      upper_limit = res_bin_start-(bincrement * (bin+1));
    }
    
  }
  else if ((bin > 0) && (bin < res_bins))
  {
    lower_limit = 1.0/sqrt(resord_model.ordinal( static_cast<double>(bin)/static_cast<double>(res_bins) ) );
    upper_limit = 1.0/sqrt(resord_model.ordinal( static_cast<double>(bin+1)/static_cast<double>(res_bins) ) );

    if (resolution_range_override)
    {
      lower_limit = res_bin_start-(bincrement * bin);
      upper_limit = res_bin_start-(bincrement * (bin+1));
    }
  }
  else
  {
    lower_limit = reso_min;//res_bin_start-(((float)(res_bin_start - res_bin_end)/(float)res_bins) * bin);;
    upper_limit = 1.0/sqrt(resord_model.ordinal( 1.0 ));

    if (resolution_range_override)
    {
      lower_limit = res_bin_start;
      upper_limit = res_bin_end;
    }
  }

  if (bin == res_bins)
  {  
    outlog << "┠─────┼──────────────┼────────┼───────────────┼────────┼───────────────┼───────────────┨" << std::endl;

  }
    outlog << std::fixed;
    if (bin == res_bins)
    {outlog << "┃TOTAL│ ";} else
    {outlog << "┃ " << std::setprecision(0) << std::setw(3) << (bin+1) << " │ ";}
    outlog << std::setprecision(2) << std::setw(5) << lower_limit << " >" <<
              std::setprecision(2) << std::setw(5) << upper_limit << " │" <<
              std::setprecision(0) << std::setw(7) << data_i_counts_bins[bin] << " │" <<
              std::setprecision(0) << std::setw(14) << data_i_bins[bin] << " │" <<
              std::setprecision(0) << std::setw(7) << model_i_counts_bins[bin] << " │" <<
              std::setprecision(0) << std::setw(14) << model_i_bins[bin] <<  " │" <<
              std::setprecision(6) << std::setw(14) << (data_i_bins[bin] / model_i_bins[bin]) << " ┃" <<
    std::endl;
  }

  outlog << "●━━━━━┷━━━━━━━━━━━━━━┷━━━━━━━━┷━━━━━━━━━━━━━━━┷━━━━━━━━┷━━━━━━━━━━━━━━━┷━━━━━━━━━━━━━━━●" << std::endl;

  outlog << " Wilson Model < I(h) > =  < Ω(h)^2/k^2 * DWF(h)^2 * ε(h) * ∑f(h)^2 > where: " << std::endl;
  outlog << " Ω(h) is the Gaussian apodization function applied to the structure factor amplitudes " << std::endl; 
  outlog << " k is the linear scale factor that puts the structure factor amplitudes on an absolute scale " << std::endl;
  outlog << " DWF(h) is an overall anistropic Debye-Waller factor, accounting for the displacements of atoms from their mean positions in the crystal " << std::endl;
  outlog << " ε(h) is the standard statistical weight for an observation with index h " << std::endl;  
  outlog << " ∑f(h)^2 represents the sum of the squared atomic scattering factors of the atoms in the unit cell" << std::endl;
  outlog << std::endl;
  outlog << " Data Completeness" << std::endl;
  outlog << " The low-resolution completeness  ( 25.0 - 10.0 Å) is " << std::setprecision(2) << low_res_completeness*100.0 << "%" << std::endl;
  outlog << " The high-resolution completeness ( 10.0 - " <<  std::setprecision(1) << std::setw(4) << 1.0/sqrt(resord_model.ordinal( 1.0 )) << " Å) is " << std::setprecision(2) << high_res_completeness*100.0 << "%\n" << std::endl;
  
  return;
}

std::vector<std::string> util::convert_fasta_to_globule_list(const std::string input, bool sorted)
{
  /*
  Simple function which iterates through the
  one letter-codes of the provided sequence,
  and creates an appropriate atom list.
  In this case, for globic scatterers.
  */
  std::vector<std::string> new_list;
  int add = 0; // How many elements to add.

  for (int element = 0; element < input.size(); element ++)
  {
    // Nothing fancy here, just find out how many of these suckers to add, add two for bigger AA's.
    char c = input.at(element);

    std::cout << input.at(element);

    switch(c)
    {
      case 'G' : add++;
      case 'A' : add++;
      case 'C' : add++; add++;
      case 'S' : add++;
      case 'V' : add++;
      case 'T' : add++;
      case 'P' : add++;
      case 'I' : add++;
      case 'L' : add++;
      case 'M' : add++; add++;
      case 'N' : add++;
      case 'D' : add++;
      case 'Q' : add++; add++;
      case 'E' : add++; add++;
      case 'K' : add++;
      case 'H' : add++; add++;
      case 'F' : add++; add++;
      case 'R' : add++; add++;
      case 'Y' : add++; add++;
      case 'W' : add++; add++;
      case 'X' : add++;
      default: add+=0;
    }
  }
  std::cout << " Total number of globules: " << add << std::endl;
  new_list.resize(add,"Xbb");
  return new_list;
}

void util::version_fingerprint(std::ostream &outlog)
{

  //outlog << "\n●━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━●" << std::endl;
  //outlog << " IPA: Iterative Projection Algorithms for protein crystallography. " << std::endl;
  //outlog << " Version: " << xmake_string(release_version) << " Build: " << xmake_string(svn_version) << std::endl;
  //outlog << " Compiled on " << xmake_string(compile_date) << "." << std::endl;
  //outlog << "●━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━●\n" << std::endl;

  util::print_section_line(78,outlog);
  outlog << util::parse_line(78," IPA: Iterative Projection Algorithms for protein crystallography. ");
  outlog << util::parse_line(78," Version: " + std::string(xmake_string(release_version)) + " Build: " + xmake_string(svn_version));
  outlog << util::parse_line(78," Compiled on " + std::string(xmake_string(compile_date)) + ".");
  util::print_section_line(78,outlog);

}

