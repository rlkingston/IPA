// IPA: Iterative Projection Algorithms for protein crystallography
// Core function library 

#include <iomanip>

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

#include "ipa-lib.h"

#include <cmath> // std::pow 

#include <iomanip> // std::setw & std::setprecison & std::boolalpha

#include <Eigen/Dense>

void ipa_functions::simple_linear_scale(
                                        const clipper::HKL_data<clipper::data64::F_sigF>  &fourier_amplitudes_work_set,
                                        const clipper::HKL_data<clipper::data64::F_sigF>  &fourier_amplitudes_test_set,
                                        const clipper::HKL_data<clipper::data64::Phi_fom> &apodization_weights,
                                        const clipper::HKL_data<clipper::data64::F_phi>   &working_amplitude_and_phase,
                                        const bool &use_test_set,
                                        const double &threshold,
                                        const int    &res_bins,
                                        double       &scale_factor,
                                        double       &correlation_coefficient_work_set,
                                        double       &r_factor_work_set,
                                        double       &correlation_coefficient_test_set,
                                        double       &r_factor_test_set,
                                        std::ostream &outlog,
                                        bool verbose)
{
  // The function takes a working set of amplitudes (working_amplitude_and_phase)
  // and computes the linear scale factor which puts them on the same scale as a second set of amplitudes (fourier_amplitudes_work_set)
  // It is expected that both sets of amplitudes have been apodized in the same fashion, according to the apodization weights
  // The apodization weights passed into the function are simply used for reporting purposes, and not applied to the data. 
  
  // Observations where the amplitudes < (threshold * the estimated standard deviation) aren't used to calculate the linear scale factor
  // The scale factor, plus the overall r-factor and correlation coefficient between the two sets of amplitudes are returned
  // If requested (use_test_set is true) then statistics are also returned for the test set
  
  
  const clipper::Spacegroup spacegroup  = fourier_amplitudes_work_set.base_hkl_info().spacegroup();
  const clipper::Cell       cell        = fourier_amplitudes_work_set.base_hkl_info().cell();
  const clipper::Resolution resolution  = fourier_amplitudes_work_set.base_hkl_info().resolution();
  
  clipper::HKL_info hkl(spacegroup, cell, resolution);
  hkl.generate_hkl_list(); // initially there are no observations hkl in the list - better make some
  
  // Use standard linear regression through the origin.
  // TODO ... Take proper account of the likely errors in the data when computing the linear scale factor.
  // Moreno C. A least-squares-based method for determining the ratio between two measured quantities. Measurement Science and Technology. 1996 Jan 1;7(2):137–41.
  
  std::vector<double> sum_abs_diff_x_y(res_bins + 1, 0.0);
  std::vector<double> sum_y(res_bins + 1, 0.0);
  
  std::vector<double> sumxx(res_bins + 1, 0.0);
  std::vector<double> sumxy(res_bins + 1, 0.0);
  std::vector<double> sumyy(res_bins + 1, 0.0);
  
  std::vector<double> sumwa(res_bins + 1, 0.0);
  
  
  std::vector<int> observation_count_work_set(res_bins + 1, 0);
  std::vector<int> observation_count_test_set(res_bins + 1, 0);
  
  std::vector<double> cc_work_set(res_bins + 1, 0.0);
  std::vector<double> rf_work_set(res_bins + 1, 0.0);
  
  std::vector<double> cc_test_set(res_bins + 1, 0.0);
  std::vector<double> rf_test_set(res_bins + 1, 0.0);
  
  clipper::HKL_data_base::HKL_reference_index ih;
  
  clipper::Resolution_ordinal resord; // Resolution Ordinal - Once initiated this returns the approximate fractional ordinal within a dataset (in the range 0...1) for any specified value of 1 / s^2
  resord.init( hkl, 1.0 );
  
  
  // First pass through the data, compute  the overall linear scale factor
  
  scale_factor = 1.0;
  
  double sumxx_s=0.0;
  double sumxy_s=0.0;
  int n_scale = 0;
  
  for ( ih = working_amplitude_and_phase.first(); !ih.last(); ih.next() )
  {
    
    double w = 1.0 / ih.hkl_class().epsilon(); // Get the standard statistical weights
    
    // If the term is observed; Fmeas >= threshold*SigFmeas, and the apodization weight > 0.05 then use this term to compute a linear scale factor (slope)
    
    if ( ( !fourier_amplitudes_work_set[ih].missing() ) &&
        ( fourier_amplitudes_work_set[ih].f() >= (threshold*fourier_amplitudes_work_set[ih].sigf()) ) &&
        ( apodization_weights[ih].fom() >= 0.05) )
    {
      n_scale  += 1;
      sumxx_s +=   w *      working_amplitude_and_phase[ih].f() *  working_amplitude_and_phase[ih].f();
      sumxy_s +=   w *      working_amplitude_and_phase[ih].f() *  fourier_amplitudes_work_set[ih].f();
    }
  }
  
  // compute overall scale factor
  
  
  if (n_scale > 0) scale_factor = sumxy_s / sumxx_s;
  
  if (verbose)
  {
    outlog << "Number of Fourier amplitudes used for computing the linear scale factor: " << n_scale << std::endl;
    outlog << "Linear scale factor: " << scale_factor << std::endl;
  }
  
  // Second pass through the data, compute the R-factor and correlation coefficient as a function of resolution
  
  // Working Set
  
  for ( ih = working_amplitude_and_phase.first(); !ih.last(); ih.next() )
  {
    
    double w = 1.0 / ih.hkl_class().epsilon(); // Get the standard statistical weights associated with each observation
    double wa = apodization_weights[ih].fom(); // Get the weights associated with the apodization function
    
    if (! fourier_amplitudes_work_set[ih].missing())
    {
      
      int bin = std::min ( int(static_cast<double>(res_bins) * resord.ordinal(ih.invresolsq())) , res_bins-1 ); // Which bin does this observation belong to ?
      
      // statistics as a function of resolution bin
      
      observation_count_work_set[bin] += 1;
      sum_abs_diff_x_y[bin]  +=   w * std::abs(scale_factor * working_amplitude_and_phase[ih].f() - fourier_amplitudes_work_set[ih].f() );
      sum_y[bin]             +=   w * fourier_amplitudes_work_set[ih].f();
      sumxx[bin]             +=   w * working_amplitude_and_phase[ih].f()          * working_amplitude_and_phase[ih].f();
      sumxy[bin]             +=   w * working_amplitude_and_phase[ih].f()          * fourier_amplitudes_work_set[ih].f();
      sumyy[bin]             +=   w * fourier_amplitudes_work_set[ih].f() * fourier_amplitudes_work_set[ih].f();
      sumwa[bin]             +=   wa;
      
      // overall statistics
      
      observation_count_work_set[res_bins] += 1;
      sum_abs_diff_x_y[res_bins]  +=   w * std::abs(scale_factor * working_amplitude_and_phase[ih].f() - fourier_amplitudes_work_set[ih].f() );
      sum_y[res_bins]             +=   w * fourier_amplitudes_work_set[ih].f();
      sumxx[res_bins]             +=   w * working_amplitude_and_phase[ih].f()          * working_amplitude_and_phase[ih].f();
      sumxy[res_bins]             +=   w * working_amplitude_and_phase[ih].f()          * fourier_amplitudes_work_set[ih].f();
      sumyy[res_bins]             +=   w * fourier_amplitudes_work_set[ih].f() * fourier_amplitudes_work_set[ih].f();
      sumwa[res_bins]             +=   wa;
      
    }
    
  }
  
  // compute the correlation coefficient and r-factor
  for ( int bin = 0; bin <= (res_bins); bin++ )
  {
    if (observation_count_work_set[bin] > 0)
    {
      cc_work_set[bin] = sumxy[bin] / sqrt(sumxx[bin]*sumyy[bin]);
      rf_work_set[bin] = sum_abs_diff_x_y[bin] / sum_y[bin];
    }
  }
  
  // Test Set
  
  if (use_test_set)
  {
    sum_abs_diff_x_y.clear();
    sum_y.clear();
    sumwa.clear();
    sumxx.clear();
    sumxy.clear();
    sumyy.clear();
    
    sum_abs_diff_x_y.resize(res_bins + 1, 0.0);
    sum_y.resize(res_bins + 1, 0.0);
    sumwa.resize(res_bins + 1, 0.0);
    sumxx.resize(res_bins + 1, 0.0);
    sumxy.resize(res_bins + 1, 0.0);
    sumyy.resize(res_bins + 1, 0.0);
    
    for ( ih = working_amplitude_and_phase.first(); !ih.last(); ih.next() )
    {
      
      double w = 1.0 / ih.hkl_class().epsilon(); // Get the standard statistical weights associated with each observation
      double wa = apodization_weights[ih].fom(); // Get the weights associated with the apodization function
      
      if (! fourier_amplitudes_test_set[ih].missing())
      {
        int bin = std::min ( int(static_cast<double>(res_bins) * resord.ordinal(ih.invresolsq())) , res_bins-1 ); // Which bin does this observation belong to ?
        
        // statistics as a function of resolution bin
        
        observation_count_test_set[bin] += 1;
        sum_abs_diff_x_y[bin]  +=   w * std::abs(scale_factor * working_amplitude_and_phase[ih].f() - fourier_amplitudes_test_set[ih].f() );
        sum_y[bin]             +=   w * fourier_amplitudes_test_set[ih].f();
        sumxx[bin]             +=   w * working_amplitude_and_phase[ih].f()          * working_amplitude_and_phase[ih].f();
        sumxy[bin]             +=   w * working_amplitude_and_phase[ih].f()          * fourier_amplitudes_test_set[ih].f();
        sumyy[bin]             +=   w * fourier_amplitudes_test_set[ih].f() * fourier_amplitudes_test_set[ih].f();
        sumwa[bin]             +=   wa;
        
        // overall statistics
        
        observation_count_test_set[res_bins] += 1;
        sum_abs_diff_x_y[res_bins]  +=   w * std::abs(scale_factor * working_amplitude_and_phase[ih].f() - fourier_amplitudes_test_set[ih].f() );
        sum_y[res_bins]             +=   w * fourier_amplitudes_test_set[ih].f();
        sumxx[res_bins]             +=   w * working_amplitude_and_phase[ih].f()          * working_amplitude_and_phase[ih].f();
        sumxy[res_bins]             +=   w * working_amplitude_and_phase[ih].f()          * fourier_amplitudes_test_set[ih].f();
        sumyy[res_bins]             +=   w * fourier_amplitudes_test_set[ih].f() * fourier_amplitudes_test_set[ih].f();
        sumwa[res_bins]             +=   wa;
        
      }
      
    }
    
    // compute the correlation coefficient and r-factor
    for ( int bin = 0; bin <= (res_bins); bin++ )
    {
      if (observation_count_test_set[bin] > 0)
      {
        cc_test_set[bin] = sumxy[bin] / sqrt(sumxx[bin]*sumyy[bin]);
        rf_test_set[bin] = sum_abs_diff_x_y[bin] / sum_y[bin];
      }
    }
    
  }
  
  // output the agreement statistics
  
  resord.invert(); // invert the distribution so that now we specify the ordinal value and it returns 1/|s|**2
  
  double lower_limit;
  double upper_limit;
  
  //
  //outlog << " |s| lo   |s| hi  1/|s| hi         #    mean W   cc       R-factor " << std::endl;
  if (verbose)
  {
    outlog << "Amplitude Agreement" << std::endl;
    // Header record
    outlog << "   |s| lo    "
    << "   |s| hi    "
    << " 1/|s| hi    "
    << "   mean W    "
    << " nobs(work)  "
    << "   cc(work)  "
    << "   RF(work)  ";
    if (use_test_set)
    {
      outlog << " nobs(test)  "
      << "   cc(test)  "
      << "   RF(test)  ";
    }
  
  outlog << std::endl;

  
  // statistics
  
  for ( int bin = 0; bin <= res_bins; bin++ )
  {
    if (bin == 0)
    {
      lower_limit = 0.0;
      upper_limit = sqrt( resord.ordinal( static_cast<double>(bin+1)/static_cast<double>(res_bins) ) );
    }
    else if ((bin > 0) && (bin < res_bins))
    {
      lower_limit = sqrt( resord.ordinal( static_cast<double>(bin)/static_cast<double>(res_bins) ) );
      upper_limit = sqrt( resord.ordinal( static_cast<double>(bin+1)/static_cast<double>(res_bins) ) );
    }
    else
    {
      lower_limit = 0.0;
      upper_limit = sqrt(resord.ordinal( 1.0 ));
    }
    
    outlog << std::fixed
    << std::setw(12) << lower_limit << " "
    << std::setw(12) << upper_limit << " "
    << std::setw(12) << 1/upper_limit <<  " "
    << std::setw(12) << sumwa[bin]/(float)observation_count_work_set[bin] <<  " "
    << std::setw(12) << observation_count_work_set[bin] <<  " "
    << std::setw(12) << cc_work_set[bin] <<  " "
    << std::setw(12) << rf_work_set[bin] << " ";
    if (use_test_set)
    {
      outlog
      << std::setw(12) << observation_count_test_set[bin] <<  " "
      << std::setw(12) << cc_test_set[bin] <<  " "
      << std::setw(12) << rf_test_set[bin] << " ";
    }
    outlog <<  std::endl;
  }
  }  
  
  // return the correlation coefficient and r-factor for the work set
  
  correlation_coefficient_work_set = cc_work_set[res_bins];
  r_factor_work_set = rf_work_set[res_bins];
  
  // return the correlation coefficient and r-factor for the test set
  
  if (use_test_set)
  {
    correlation_coefficient_test_set = cc_test_set[res_bins];
    r_factor_test_set = rf_test_set[res_bins];
  }
  else
  {
    correlation_coefficient_test_set = 0.0;
    r_factor_test_set = 0.0;
  }
  
}

void ipa_functions::compute_fourier_coeffs_first_order_deriv (const clipper::HKL_data<clipper::data64::F_phi> &fp, clipper::HKL_data<clipper::data64::F_phi> &fp_dpdx, clipper::HKL_data<clipper::data64::F_phi> &fp_dpdy, clipper::HKL_data<clipper::data64::F_phi> &fp_dpdz)
{
  
  
  const clipper::Cell cell = fp.base_hkl_info().cell();
  
  //  double degtorad = atan(1.0)/45.0;
  
  /*
   clipper::Mat33<double> fm = hkl1.cell().matrix_frac();
   outlog << fm(0,0) << std::endl;
   outlog << fm(0,1) << std::endl;
   outlog << fm(0,2) << std::endl;
   outlog << fm(1,0) << std::endl;
   outlog << fm(1,1) << std::endl;
   outlog << fm(1,2) << std::endl;
   outlog << fm(2,0) << std::endl;
   outlog << fm(2,1) << std::endl;
   outlog << fm(2,2) << std::endl;
   */
  
  // double pi = acos(-1);
  
  // Here are the elements of the orthogonalization (change of basis) matrix.
  
  double m11 =  cell.a_star() * sin(cell.beta_star()) * sin(cell.gamma()); // = 1/a
  double m21 = -cell.a_star() * sin(cell.beta_star()) * cos(cell.gamma()); // = -cosγ/(asinγ)
  double m31 =  cell.a_star() * cos(cell.beta_star());
  double m12 =  0;
  double m22 =  cell.b_star() * sin(cell.alpha_star());
  double m32 =  cell.b_star() * cos(cell.alpha_star());
  double m13 = 0.0;
  double m23 = 0.0;
  double m33 =  cell.c_star();
  
  // Under a change of basis,  the first partial derivatives of a function expressed in crystallographic coordinates transforms in the same way as the basis vectors themselves
  // (they are covariant with the basis vectors). Hence we use the matrix appropriate for transformation of the basis vectors, to transform the vector of partial derivatives.
  
  
  /*
   outlog <<  m11 << std::endl;
   outlog <<  m21 << std::endl;
   outlog << m31 << std::endl;
   outlog << m22 << std::endl;
   outlog << m32 << std::endl;
   outlog << m33 << std::endl;
   */
  
  clipper::HKL_data_base::HKL_reference_index ih;
  for ( ih = fp_dpdx.first(); !ih.last(); ih.next() )
  {
    float h = float(ih.hkl().h());
    float k = float(ih.hkl().k());
    float l = float(ih.hkl().l());
    
    // Compute the purely imaginary multiplier required & put it into polar form
    
    fp_dpdx[ih].f() = -clipper::Util::twopi()*(h*m11);
    fp_dpdy[ih].f() = -clipper::Util::twopi()*(h*m21 + k*m22);
    fp_dpdz[ih].f() = -clipper::Util::twopi()*(h*m31 + k*m32 + l*m33);
    
    if (fp_dpdx[ih].f() > 0)
    {
      fp_dpdx[ih].phi() = clipper::Util::pi()/2;
    }
    else
    {
      fp_dpdx[ih].phi() = 3*clipper::Util::pi()/2;
      fp_dpdx[ih].f() = std::abs(fp_dpdx[ih].f());
    }
    
    if (fp_dpdy[ih].f() > 0)
    {
      fp_dpdy[ih].phi() = clipper::Util::pi()/2;
    }
    else
    {
      fp_dpdy[ih].phi() = 3*clipper::Util::pi()/2;
      fp_dpdy[ih].f() = std::abs(fp_dpdy[ih].f());
    }
    
    if (fp_dpdz[ih].f() > 0)
    {
      fp_dpdz[ih].phi() = clipper::Util::pi()/2;
    }
    else
    {
      fp_dpdz[ih].phi() = 3*clipper::Util::pi()/2;
      fp_dpdz[ih].f() = std::abs(fp_dpdz[ih].f());
    }
    
    // Now multiply by the structure factors to get the required Fourier coefficients. We are in polar form, so just multiply the amplitudes and sum the phases.
    
    //      outlog << h << " " << k << " " << l << " " << std::endl;
    //      outlog << fp_dpdz[ih].f() << " " << fp_dpdz[ih].phi() << std::endl;
    //      outlog << fp[ih].f() << " " <<   fp[ih].phi() << std::endl;
    
    fp_dpdx[ih].f() *= fp[ih].f();
    fp_dpdy[ih].f() *= fp[ih].f();
    fp_dpdz[ih].f() *= fp[ih].f();
    
    fp_dpdx[ih].phi() = fmod(fp_dpdx[ih].phi() +  fp[ih].phi(), clipper::Util::twopi());
    fp_dpdy[ih].phi() = fmod(fp_dpdy[ih].phi() +  fp[ih].phi(), clipper::Util::twopi());
    fp_dpdz[ih].phi() = fmod(fp_dpdz[ih].phi() +  fp[ih].phi(), clipper::Util::twopi());
    
    //      outlog << fp_dpdz[ih].f() << " " <<  fp_dpdz[ih].phi() << std::endl;
    
  }
  
  
}

void ipa_functions::compute_fourier_coeffs_second_order_deriv(const clipper::HKL_data<clipper::data64::F_phi> &fp, clipper::HKL_data<clipper::data64::F_phi> &fp_lapl)
{
  
  const clipper::Cell cell = fp.base_hkl_info().cell();
  
  // double pi = acos(-1);
  
  double m11 =  cell.a_star() * sin(cell.beta_star()) * sin(cell.gamma());
  double m21 = -cell.a_star() * sin(cell.beta_star()) * cos(cell.gamma());
  double m31 =  cell.a_star() * cos(cell.beta_star());
  double m22 =  cell.b_star() * sin(cell.alpha_star());
  double m32 =  cell.b_star() * cos(cell.alpha_star());
  double m33 =  cell.c_star();
  
  clipper::HKL_data_base::HKL_reference_index ih;
  for ( ih = fp_lapl.first(); !ih.last(); ih.next() )
  {
    float h = float(ih.hkl().h());
    float k = float(ih.hkl().k());
    float l = float(ih.hkl().l());
    
    // Compute the purely real multiplier required & put it into polar form
    
    fp_lapl[ih].f() = -4*clipper::Util::twopi()*( ( m11* m11 +  m21*m21 + m31*m31)*(h*h) + (m22*m22 + m32*m32)*(k*k) + (m33*m33)*(l*l) + 2*(m21*m22 + m31*m32)*(h*k) + 2*m31*m33*h*l + 2*m32*m33*k*l );
    
    if (fp_lapl[ih].f() > 0)
    {
      fp_lapl[ih].phi() = 0;
    }
    else
    {
      fp_lapl[ih].phi() = clipper::Util::pi()/2;
      fp_lapl[ih].f() = std::abs(fp_lapl[ih].f());
    }
    
    // Now multiply by the structure factors to get the required Fourier coefficients. We are in polar form, so just multiply the amplitudes and sum the phases.
    
    fp_lapl[ih].f() *= fp[ih].f();
    
    fp_lapl[ih].phi() = fmod(fp_lapl[ih].phi() +  fp[ih].phi(), clipper::Util::twopi());
  }
}

void ipa_functions::calculate_gradient(const clipper::HKL_data<clipper::data64::F_phi> &fp, clipper::Xmap<float>  &gradient_x, clipper::Xmap<float>  &gradient_y, clipper::Xmap<float>  &gradient_z, clipper::Xmap<float>  &gradient_magnitude, std::ostream &outlog, bool verbose)
{
  
  // NB - I think there are some issues associated with finite sampling here - need to do some reading and thinking about the calculation of derivatives using an FFT.
  
  const clipper::Spacegroup spacegroup  = fp.base_hkl_info().spacegroup();
  const clipper::Cell       cell        = fp.base_hkl_info().cell();
  const clipper::Resolution resolution  = fp.base_hkl_info().resolution();
  
  // Make data objects to hold the Fourier coefficients for the partial derivatives
  
  clipper::HKL_info hkl(spacegroup, cell, resolution);
  hkl.generate_hkl_list(); // initially there are no observations hkl in the list - better make some
  
  clipper::HKL_data<clipper::data64::F_phi> fp_dpdx (hkl); // create a data object of type  F + phi to hold Fourier coefficients for first partial derivative of the electron density wrt orthogonal axis X
  clipper::HKL_data<clipper::data64::F_phi> fp_dpdy (hkl); // create a data object of type  F + phi to hold Fourier coefficients for first partial derivative of the electron density wrt orthogonal axis Y
  clipper::HKL_data<clipper::data64::F_phi> fp_dpdz (hkl); // create a data object of type  F + phi to hold Fourier coefficients for first partial derivative of the electron density wrt orthogonal axis Z
  if (verbose)
  {
  outlog << "Constructing the Fourier coefficients to allow calculation of the gradient" << std::endl;
  }
  ipa_functions::compute_fourier_coeffs_first_order_deriv(fp, fp_dpdx, fp_dpdy, fp_dpdz); // calculate the Fourier coefficients required to evaluate the first derivative of the electron density with respect to an orthonormal basis
  
  if (verbose)
  {
  outlog << "Calculating the gradient components along orthormal axes X, Y and Z" << std::endl;
  }
  gradient_x.fft_from( fp_dpdx ); // Backtransform to generate first partial derivative of the electron density wrt orthogonal axis X
  gradient_y.fft_from( fp_dpdy ); // Backtransform to generate first partial derivative of the electron density wrt orthogonal axis Y
  gradient_z.fft_from( fp_dpdz ); // Backtransform to generate first partial derivative of the electron density wrt orthogonal axis Z
  
  if (verbose)
  {
  outlog << "Calculating the gradient magnitude " << std::endl;
  }
  
  clipper::Xmap_base::Map_reference_index ix;
  
  for ( ix =   gradient_magnitude.first(); !ix.last(); ix.next() )
  {
    gradient_magnitude[ix] = sqrt( std::pow(gradient_x[ix],2) +  std::pow(gradient_y[ix],2) + std::pow(gradient_z[ix],2) );
  }
  
}

void ipa_functions::mask_connections (const int &mask_id, const bool &apply_translational_symm, const bool &check_diagonals, const float &fractional_threshold, clipper::Xmap<int> &mask, int &n_connected_sets_eliminated, int &n_connected_sets_retained, std::ostream &outlog, bool verbose)

{
  // Takes a mask in space group P1 and creates an anchor set map (a mask with every locally connected set of points assigned a unique positive integer)
  // Then computes all connected sets from the anchor sets, and writes a mask back out, removing the smallest connected sets
  
  // This routine requires a binary valued integer mask as input.
  // Points in the mask must be denoted with an integer > 0. This is specified by mask_id
  // Points not in the mask must be denoted with an integer <= 0
  // Input mask values are overwritten with the anchor set ID
  
  // Spectacular failure will result if these input conditions are not met.
  
  // Loop over the unit cell of the P1 expanded mask and convert in into the anchor_set_map
  
  clipper::Xmap_base::Map_reference_coord is;
  clipper::Xmap_base::Map_reference_coord i_out;
  clipper::Xmap_base::Map_reference_coord i_mid;
  clipper::Xmap_base::Map_reference_coord i_inn;
  
  clipper::Coord_grid gs(0,0,0); // starting grid coordinates
  clipper::Coord_grid gf(mask.grid_sampling().nu()-1, mask.grid_sampling().nv()-1, mask.grid_sampling().nw()-1); // ending grid coordinates
  is = clipper::Xmap_base::Map_reference_coord( mask, gs );
  
  
  int current_anchor_set;
  int working_anchor_set;
  int number_of_edge_links;
  
  Eigen::VectorXi anchor_set_size; // Dynamic-size Vector of integers
  Eigen::MatrixXi connection_matrix; // Dynamic-size Matrix of integers
  connection_matrix = Eigen::MatrixXi::Constant(1,1,1); // initialize matrix of size 1x1, with entry 1
  anchor_set_size = Eigen::VectorXi::Constant(1,0); // initialize vector of size 1, with entry 0
  
  char io_map_filename[64];
  
  // Loop over the full unit cell, one row of grid points at a time.
  // Grid points are assigned to locally contiguous sets (anchor sets) which are written into the mask, making it an anchor set map
  // Spatial connections between the anchor sets are recorded in the binary connection matrix.
  if (verbose)
  {
    outlog << "\nConstructing Anchor sets and evaluating connections within the unit cell" << std::endl;
    if (check_diagonals)
    {
      outlog << "Edge and Diagonal connections will be checked" << std::endl;
    }
    else
    {
      outlog << "Only Edge connections will be checked" << std::endl;
    }
  }
  
  int points_in_mask = 0;
  int number_of_anchor_sets = 0;
  
  
  current_anchor_set = 0;
  for ( i_out = is; i_out.coord().u() <= gf.u(); i_out.next_u() ) // outer loop in u
  {
    current_anchor_set = 0;
    for ( i_mid = i_out; i_mid.coord().v() <= gf.v(); i_mid.next_v() ) // middle loop in v
    {
      current_anchor_set = 0;
      for ( i_inn = i_mid; i_inn.coord().w() <= gf.w(); i_inn.next_w() ) // inner loop in w
      {
        
        
        clipper::Xmap_base::Map_reference_coord p000;
        p000 = i_inn;
        
        //        outlog << p000.coord().u() << " " << p000.coord().v() << " " << p000.coord().w() <<  std::endl;
        //        outlog << "connection_matrix =" << std::endl << connection_matrix << std::endl;
        //        outlog << "number_of_anchor_sets = " << number_of_anchor_sets << std::endl;
        
        if (!(mask[p000] == mask_id )) // the current point is not within the mask ... ignore it.
        {
          current_anchor_set = 0;
        }
        else                                     // the current point is within the mask ... process it.
        {
          points_in_mask += 1;
          
          //          outlog << "points_in_mask" << " " << points_in_mask << std::endl;
          //          outlog << "current_anchor_set" << " " <<  current_anchor_set << std::endl;
          
          if (check_diagonals)
            // slow procedure, checks both edge and diagonal connections
            
          {
            
            // construct the coordinates of all the nodes we need to check
            
            clipper::Xmap_base::Map_reference_coord p100;
            p100 = p000;
            p100.prev_u();
            clipper::Xmap_base::Map_reference_coord p010;
            p010 = p000;
            p010.prev_v();
            clipper::Xmap_base::Map_reference_coord p001;
            p001 = p000;
            p001.prev_w();
            clipper::Xmap_base::Map_reference_coord p110;
            p110 = p000;
            p110.prev_u();
            p110.prev_v();
            clipper::Xmap_base::Map_reference_coord p101;
            p101 = p000;
            p101.prev_u();
            p101.prev_w();
            clipper::Xmap_base::Map_reference_coord p011;
            p011 = p000;
            p011.prev_v();
            p011.prev_w();
            clipper::Xmap_base::Map_reference_coord p111;
            p111 = p000;
            p111.prev_u();
            p111.prev_v();
            p111.prev_w();
            
            // If we are at the edge of the unit cell we don't check certain nodes ... work this out
            // NB We deal with the translational symmetry of the crystal after we work out the connectivity *within* the unit cell
            
            bool check_100 = true;
            bool check_010 = true;
            bool check_001 = true;
            bool check_110 = true;
            bool check_101 = true;
            bool check_011 = true;
            bool check_111 = true;
            
            if (p000.coord().u() == is.coord().u()) // we are at the edge of the cell in u
            {
              check_100 = false;
              check_110 = false;
              check_101 = false;
              check_111 = false;
            }
            
            if (p000.coord().v() ==  is.coord().v()) // we are at the edge of the cell in v
            {
              check_010 = false;
              check_110 = false;
              check_011 = false;
              check_111 = false;
            }
            
            if (p000.coord().w() == is.coord().w()) // we are at the edge of the cell in w
            {
              check_001 = false;
              check_101 = false;
              check_011 = false;
              check_111 = false;
            }
            
            /*
             outlog <<  std::endl;
             outlog << p000.coord().u() << " " << p000.coord().v() << " " << p000.coord().w() <<  std::endl;
             outlog <<  std::endl;
             outlog << p100.coord().u() << " " << p100.coord().v() << " " << p100.coord().w() <<  std::endl;
             outlog << check_100 << std::endl;
             outlog << p010.coord().u() << " " << p010.coord().v() << " " << p010.coord().w() <<  std::endl;
             outlog << check_010 << std::endl;
             outlog << p001.coord().u() << " " << p001.coord().v() << " " << p001.coord().w() <<  std::endl;
             outlog << check_001 << std::endl;
             outlog << p110.coord().u() << " " << p110.coord().v() << " " << p110.coord().w() <<  std::endl;
             outlog << check_110 << std::endl;
             outlog << p101.coord().u() << " " << p101.coord().v() << " " << p101.coord().w() <<  std::endl;
             outlog << check_101 << std::endl;
             outlog << p011.coord().u() << " " << p011.coord().v() << " " << p011.coord().w() <<  std::endl;
             outlog << check_011 << std::endl;
             outlog << p111.coord().u() << " " << p111.coord().v() << " " << p111.coord().w() <<  std::endl;
             outlog << check_111 << std::endl;
             
             */
            
            number_of_edge_links = 0;
            
            bool link_100 = false;
            bool link_010 = false;
            bool link_001 = false;
            //          bool link_110 = false;
            //          bool link_101 = false;
            //          bool link_011 = false;
            //          bool link_111 = false;
            
            // First check the 3 edge connections
            
            //          outlog << "Checking edge connections" <<  std::endl;
            
            if (check_100)
            {
              working_anchor_set = mask[p100];
              if (working_anchor_set > 0) // if true this is a link
              {
                link_100 = true;
                number_of_edge_links += 1;
                update_connection_matrix(current_anchor_set, working_anchor_set, connection_matrix);
              }
            }
            
            if (check_010)
            {
              working_anchor_set = mask[p010];
              if (working_anchor_set > 0) // if true this is a link
              {
                link_010 = true;
                number_of_edge_links += 1;
                update_connection_matrix(current_anchor_set, working_anchor_set, connection_matrix);
              }
            }
            
            if (check_001)
            {
              working_anchor_set = mask[p001];
              if (working_anchor_set > 0) // if true this is a link
              {
                link_001 = true;
                number_of_edge_links += 1;
                update_connection_matrix(current_anchor_set, working_anchor_set, connection_matrix);
              }
            }
            
            //          outlog << "number_of_edge_links" << " " << number_of_edge_links << std::endl;
            
            // if there are two or three edge connections, then no unique diagonal connections can exist ... we are done
            // If there are zero or one edge connections we have to examine the diagonals
            
            // one edge connection  - there is at most one unique diagonal connection, and we know where to look
            
            
            if (number_of_edge_links == 1)
            {
              if ( (link_100) && (check_011) ) // if the relvant diagonal node exists then check for link
              {
                working_anchor_set = mask[p011];
                if (working_anchor_set > 0) update_connection_matrix(current_anchor_set, working_anchor_set, connection_matrix);
              }
              if ( (link_010) && (check_101) ) // if the relvant diagonal node exists then check for link
              {
                working_anchor_set = mask[p101];
                if (working_anchor_set > 0) update_connection_matrix(current_anchor_set, working_anchor_set, connection_matrix);
              }
              if ( (link_001) && (check_110) ) // if the relvant diagonal node exists then check for link
              {
                working_anchor_set = mask[p110];
                if (working_anchor_set > 0) update_connection_matrix(current_anchor_set, working_anchor_set, connection_matrix);
              }
              
            }
            
            // no edge connections - there is at most one unique diagonal connection, but we do not know where to look
            
            
            if (number_of_edge_links == 0)
            {
              bool found = false;
              if(            (check_110))
              {
                working_anchor_set = mask[p110];
                if (working_anchor_set > 0) // if true this is a link
                {
                  found = true;
                  update_connection_matrix(current_anchor_set, working_anchor_set, connection_matrix);
                }
              }
              if(!(found) && (check_101))
              {
                working_anchor_set = mask[p101];
                if (working_anchor_set > 0) // if true this is a link
                {
                  found = true;
                  update_connection_matrix(current_anchor_set, working_anchor_set, connection_matrix);
                }
              }
              if(!(found) && (check_011))
              {
                working_anchor_set = mask[p011];
                if (working_anchor_set > 0) // if true this is a link
                {
                  found = true;
                  update_connection_matrix(current_anchor_set, working_anchor_set, connection_matrix);
                }
              }
              if(!(found) && (check_111))
              {
                working_anchor_set = mask[p111];
                if (working_anchor_set > 0) // if true this is a link
                {
                  found = true;
                  update_connection_matrix(current_anchor_set, working_anchor_set, connection_matrix);
                }
              }
              
            }
          }
          
          else
            // fast procedure, only edge connections
            
          {
            
            
            // If we are at the edge of the unit cell we don't check certain nodes ... work this out
            // NB We deal with the translational symmetry of the crystal after we work out the connectivity *within* the unit cell
            
            
            bool check_100 = true;
            if (p000.coord().u() == is.coord().u()) // we are at the edge of the cell in u
            {
              check_100 = false;
            }
            
            bool check_010 = true;
            if (p000.coord().v() ==  is.coord().v()) // we are at the edge of the cell in v
            {
              check_010 = false;
            }
            
            bool check_001 = true;
            if (p000.coord().w() == is.coord().w()) // we are at the edge of the cell in w
            {
              check_001 = false;
            }
            
            // Check the 3 edge connections
            
            if (check_100)
            {
              clipper::Xmap_base::Map_reference_coord p100;
              p100 = p000;
              p100.prev_u();
              working_anchor_set = mask[p100];
              if (working_anchor_set > 0) // if true this is a link
              {
                update_connection_matrix(current_anchor_set, working_anchor_set, connection_matrix);
              }
            }
            
            if (check_010)
            {
              clipper::Xmap_base::Map_reference_coord p010;
              p010 = p000;
              p010.prev_v();
              working_anchor_set = mask[p010];
              if (working_anchor_set > 0) // if true this is a link
              {
                update_connection_matrix(current_anchor_set, working_anchor_set, connection_matrix);
              }
            }
            
            if (check_001)
            {
              clipper::Xmap_base::Map_reference_coord p001;
              p001 = p000;
              p001.prev_w();
              working_anchor_set = mask[p001];
              if (working_anchor_set > 0) // if true this is a link
              {
                update_connection_matrix(current_anchor_set, working_anchor_set, connection_matrix);
              }
            }
            
          }
          // Now we have now effectively checked all the neighboring nodes for connections (edge and diagonal, or edge only, as requested).
          
          if (current_anchor_set == 0) // If this is true, it's either the very first point we have examined, or there are no connections. Either way, this must be a new anchor set
          {
            number_of_anchor_sets += 1;
            current_anchor_set = number_of_anchor_sets;
            if (number_of_anchor_sets > 1)
            {
              connection_matrix.conservativeResizeLike(Eigen::MatrixXi::Zero(number_of_anchor_sets,number_of_anchor_sets));
              connection_matrix(number_of_anchor_sets-1,number_of_anchor_sets-1) = 1;
              
              anchor_set_size.conservativeResizeLike(Eigen::VectorXi::Zero(number_of_anchor_sets));
            }
            
          }
          
          mask[p000] = current_anchor_set;
          anchor_set_size(current_anchor_set-1) += 1;
        }
        
      } // end of the inner loop
    } // end of in middle loop
  } // end of the outer loop
  
  if (verbose)
  {
    outlog << "Total number of Points in Mask = " << points_in_mask << std::endl;
    outlog << "Number of Anchor sets = " << number_of_anchor_sets << std::endl;
  }
  /*
   // Write out the anchor set map
   clipper::CCP4MAPfile map_out;
   snprintf(io_map_filename, sizeof(io_map_filename), "anchor_set_map.ccp4");
   map_out.open_write( io_map_filename );
   map_out.export_xmap(mask); // write anchor_set_map
   map_out.close_write();
   
   */
  
  
  // if requested loop over the xy, xz and yz faces of the unit cell and check for the additional connections between anchor sets resulting from crystallographic translational symmetry
  // We only check edge connections here and not diagonal connections  ... this should be close enough for government work
  
  if (apply_translational_symm)
  {
    if (verbose)
    {
      outlog << "Checking for the additional connections between anchor sets resulting from crystallographic translational symmetry" << std::endl;
    }
    // xy plane
    
    for ( i_out = is; i_out.coord().u() <= gf.u(); i_out.next_u() ) // outer loop in u
    {
      for ( i_inn = i_out; i_inn.coord().v() <= gf.v(); i_inn.next_v() ) // inner loop in v
      {
        
        clipper::Xmap_base::Map_reference_coord p000;
        p000 = i_inn;
        
        if (mask[p000] > 0) // only process if the point is within the mask and has been assigned to an anchor set
        {
          
          // construct the coordinates of the single node we need to check
          
          clipper::Xmap_base::Map_reference_coord p001;
          p001 = p000;
          p001.prev_w();
          
          if (mask[p001] > 0) // current point is connected to a set across the unit cell boundary
          {
            update_connection_matrix(mask[p000], mask[p001], connection_matrix);
          }
          
        }
        
      }
    }
    
    // xz plane
    
    for ( i_out = is; i_out.coord().u() <= gf.u(); i_out.next_u() ) // outer loop in u
    {
      for ( i_inn = i_out; i_inn.coord().w() <= gf.w(); i_inn.next_w() ) // inner loop in w
      {
        
        clipper::Xmap_base::Map_reference_coord p000;
        p000 = i_inn;
        
        if (mask[p000] > 0) // only process if the point is within the mask and has been assigned to an anchor set
        {
          
          // construct the coordinates of the single node we need to check
          
          clipper::Xmap_base::Map_reference_coord p010;
          p010 = p000;
          p010.prev_v();
          
          if (mask[p010] > 0) // current point is connected to a set across the unit cell boundary
          {
            update_connection_matrix(mask[p000], mask[p010], connection_matrix);
          }
          
        }
        
      }
    }
    
    // yz plane
    
    for ( i_out = is; i_out.coord().v() <= gf.v(); i_out.next_v() ) // outer loop in v
    {
      for ( i_inn = i_out; i_inn.coord().w() <= gf.w(); i_inn.next_w() ) // inner loop in w
      {
        
        clipper::Xmap_base::Map_reference_coord p000;
        p000 = i_inn;
        
        if (mask[p000] > 0) // only process if the point is within the mask and has been assigned to an anchor set
        {
          
          // construct the coordinates of the single node we need to check
          
          clipper::Xmap_base::Map_reference_coord p100;
          p100 = p000;
          p100.prev_u();
          
          if (mask[p100] > 0) // current point is connected to a set across the unit cell boundary
          {
            update_connection_matrix(mask[p000], mask[p100], connection_matrix);
          }
          
        }
        
      }
    }
    
  }
  
  /*
   outlog << "Anchor Set Sizes" << std::endl;
   outlog << "Before Row Reduction on the binary connection matrix" << std::endl;
   
   for ( int n_set = 1; n_set <= number_of_anchor_sets; n_set++ )
   {
   outlog << n_set << " " << anchor_set_size(n_set-1) << std::endl;
   }
   
   outlog << "connection_matrix =" << std::endl << connection_matrix << std::endl;
   */
  
  // Perform a row reduction procedure on the binary connection matrix
  
  int connected_row_number;
  
  for ( int row_number = number_of_anchor_sets-1; row_number >= 0; row_number-- )
  {
    // Scan for a connection
    for ( int col_number = row_number-1; col_number >= 0; col_number-- )
    {
      // outlog << row_number << " " << col_number << std::endl;
      if (connection_matrix(row_number,col_number) == 1) // if true, then a connection has been found
      {
        // combine the connected rows
        connected_row_number = col_number;
        for ( int col_number = number_of_anchor_sets-1; col_number >= 0; col_number-- )
        {
          if (connection_matrix(row_number,col_number) == 1) connection_matrix(connected_row_number,col_number) = 1;
        }
        anchor_set_size(connected_row_number) += anchor_set_size(row_number);
        anchor_set_size(row_number) = 0;
        break; // Leave the loop scanning for a connection
      }
    }
  }
  
  /*
   outlog << "After Row Reduction on the binary connection matrix" << std::endl;
   
   for ( int n_set = 1; n_set <= number_of_anchor_sets; n_set++ )
   {
   outlog << n_set << " " << anchor_set_size(n_set-1) << std::endl;
   }
   
   outlog << "connection_matrix =" << std::endl << connection_matrix << std::endl;
   */
  
  // Create a lookup table from the connection matrix. The lookup table allows the anchor sets to be readily combined into globally connected sets
  
  
  
  int number_of_connected_sets = 0;
  Eigen::VectorXi connected_set_size; // Dynamic-size Vector of integers
  connected_set_size = Eigen::VectorXi::Constant(1,0); // initialize vector of size 1, with entry 0
  
  Eigen::VectorXi lookup_table; // Dynamic-size Vector of integers
  lookup_table = Eigen::VectorXi::Constant(number_of_anchor_sets,0); // initialize vector of size number_of_anchor_sets, with entry 0
  
  
  for ( int n_row = 0; n_row < number_of_anchor_sets; n_row++ )
  {
    if ( anchor_set_size(n_row) != 0)
    {
      number_of_connected_sets += 1;
      connected_set_size(number_of_connected_sets-1)=anchor_set_size(n_row);
      connected_set_size.conservativeResizeLike(Eigen::VectorXi::Zero(number_of_connected_sets+1));
      for ( int n_col = 0; n_col < number_of_anchor_sets; n_col++ )
      {
        if ( connection_matrix(n_row,n_col) != 0)
        {
          lookup_table(n_col)   = number_of_connected_sets;
        }
      }
    }
  }
  if (verbose)
  {
    outlog << "Number of Connected Sets =" << number_of_connected_sets << std::endl;
  }
  /*
   outlog << "Lookup table - establishing the correspondence between the Anchor sets and Connected sets " << std::endl ;
   
   for ( int n_set = 1; n_set <= number_of_anchor_sets; n_set++ )
   {
   outlog << "Anchor Set: " << n_set << " " << "Connected set: " << lookup_table[n_set-1] << std::endl;
   }
   */
  
  // Figure out which of the connected sets need to be eliminated
  // Sort the anchor sets by size and figure the threshold for elimination
  
  
  int   min_retained_set_size;
  
  // Copy the data into an array for sorting and list the connected sets
  
  int *sorted_set_size; // Declare a pointer to an array of integer data
  sorted_set_size = new int[number_of_connected_sets]; // Use the new operator to dynamically allocate space for the array
  
  
  //  outlog << "Connected Set Sizes:" << std::endl;
  
  for ( int n_set = 1; n_set <= number_of_connected_sets; n_set++ )
  {
    //   outlog << n_set << " " <<  connected_set_size(n_set-1) << std::endl;
    sorted_set_size[n_set-1] = connected_set_size(n_set-1);
  }
  
  
  std::sort(sorted_set_size, sorted_set_size + number_of_connected_sets); // Sort. Note that the Standard template library sort() function requires the end to be indicated with the address of the element *beyond* the last element that is to be sorted.
  if (verbose)
  {
    outlog << "Connected Set sizes in rank order:" << std::endl;
    for ( int n_set = 1; n_set <= number_of_connected_sets; n_set++ )
    {
      outlog << n_set << " " << sorted_set_size[n_set-1] << std::endl;
    }
  }
  min_retained_set_size = lround( float(sorted_set_size[number_of_connected_sets-1])*fractional_threshold);
  
  if (verbose)
  {
    outlog << "\nThreshold for removal of connected sets, as a percentage of the largest set size: " << fractional_threshold*100 << std::endl ;
    outlog << "Minimum size for retention of a connected set: " << min_retained_set_size << std::endl;
  }
  delete [] sorted_set_size;
  
  // get some statistics on the size and number of the eliminated sets
  
  //  int n_connected_sets_eliminated = 0;
  
  int n_points_eliminated = 0;
  n_connected_sets_eliminated = 0;
  
  for ( int n_set = 1; n_set <= number_of_connected_sets; n_set++ )
  {
    if (connected_set_size(n_set-1) < min_retained_set_size)
    {
      n_connected_sets_eliminated += 1;
      n_points_eliminated +=  connected_set_size(n_set-1);
    }
  }
  
  n_connected_sets_retained = number_of_connected_sets - n_connected_sets_eliminated;
  if (verbose)
  {
    outlog << " Number of connected sets to be eliminated: " << n_connected_sets_eliminated << std::endl;
    outlog << " Number of connected sets to be retained: " << n_connected_sets_retained  << std::endl;
    outlog << " Number of points to be eliminated from mask: " << n_points_eliminated << std::endl;
    outlog << " Percentage reduction in mask volume: " << (float(n_points_eliminated)/ float( points_in_mask ))*100 << std::endl;
  }
  
  // finally perform a single pass through the map forming a union of the largest connected sets, and eliminating the smallest connected sets
  
  
  clipper::Xmap_base::Map_reference_index ix;
  
  for ( ix = mask.first(); !ix.last(); ix.next() )
  {
    int anchor_set = mask[ix];
    
    if (anchor_set > 0) // the point is part of the mask
    {
      int connected_set = lookup_table(anchor_set-1); // lookup the connected set associated with this anchor set
      if (connected_set_size(connected_set-1) > min_retained_set_size) {mask[ix] = 1;} // if the connected set is big enough, retain it in the mask
      else                                                             {mask[ix] = 0;} // otherwise remove it from the mask
    }
  }
  
  /*
   // Write out the connected set map
   snprintf(io_map_filename, sizeof(io_map_filename), "connected_set_map.ccp4");
   map_out.open_write( io_map_filename );
   map_out.export_xmap(mask); // write anchor_set_map
   map_out.close_write();
   */
}

void ipa_functions::update_connection_matrix (int &current_anchor_set, const int &working_anchor_set, Eigen::MatrixXi &connection_matrix)

{
  if (current_anchor_set == 0)
  {
    current_anchor_set = working_anchor_set;
    //    anchor_set_size(current_anchor_set-1) += 1;
    return;
  }
  if (current_anchor_set != working_anchor_set)
  {
    connection_matrix(current_anchor_set-1, working_anchor_set-1) = 1;
    connection_matrix(working_anchor_set-1, current_anchor_set-1) = 1;
  }
}

void ipa_functions::compute_alternate_origins (const int &space_group_number, const clipper::Grid_sampling &grid, std::vector<clipper::Coord_frac> &origins)

{
  // Given the space group, numbered according to IUCR convention, returns a vector holding the fractional coordinates of all alternate origin choices in that space group
  // Based on the document authored by Ian Tickle: Alternate origins for the 65 enantiomorphic space groups
  // http://www.ccp4.ac.uk/html/alternate_origins.html
  // https://smb.slac.stanford.edu/facilities/software/ccp4/html/alternate_origins.html
  
  // In cases were the origin may be fixed arbitrarily along a certain axis, the origin choices are sampled along that axis, as dictated by the user-supplied grid_sampling
  
  // The underpinning theory is discussed at length in Chapter 2 of "Direct Phasing in Crystallography" (Giacovazzo).
  // In particular
  // Table 2.2 in Giacovazzo gives the results for the non-centrosymmetric Primitive Space groups
  // Table 2.5 in Giacovazzo gives the results for the non-centrosymmetric Centered Space groups
  //
  
  clipper::Coord_frac c0;
  
  switch (space_group_number){
      
    case 1 : // space group P1
      
      for ( int iu = 0; iu <= grid.nu(); iu++ )
      {
        for ( int iv = 0; iv <= grid.nv(); iv++ )
        {
          for ( int iw = 0; iw <= grid.nw(); iw++ )
          {
            c0 = clipper::Coord_frac( 0.0 + float(iu)/float(grid.nu() + 1),
                                     0.0 + float(iv)/float(grid.nv() + 1),
                                     0.0 + float(iw)/float(grid.nw() + 1));
            origins.push_back(c0);
          }
        }
      }
      break;
      
    case 3 : // space group P 2 (unique axis b) 
    case 4 : // space group P 21 (unique axis b)
      
      for ( int iv = 0; iv <= grid.nv(); iv++ )
      {
        c0 = clipper::Coord_frac( 0.0, 0.0 + float(iv)/float(grid.nv() + 1), 0.0  );
        origins.push_back(c0);
        c0 = clipper::Coord_frac( 0.0, 0.0 + float(iv)/float(grid.nv() + 1), 0.5  );
        origins.push_back(c0);
        c0 = clipper::Coord_frac( 0.5, 0.0 + float(iv)/float(grid.nv() + 1), 0.0  );
        origins.push_back(c0);
        c0 = clipper::Coord_frac( 0.5, 0.0 + float(iv)/float(grid.nv() + 1), 0.5 );
        origins.push_back(c0);
      }
      break;
      
    case 5 : // space group C2 (unique axis b). A2 and I2 are not dealt with here
      
      for ( int iv = 0; iv <= grid.nv(); iv++ )
      {
        c0 = clipper::Coord_frac( 0.0, 0.0 + float(iv)/float(grid.nv() + 1), 0.0  );
        origins.push_back(c0);
        c0 = clipper::Coord_frac( 0.0, 0.0 + float(iv)/float(grid.nv() + 1), 0.5  );
        origins.push_back(c0);
      }
      break;
      
    case 16 : // space group P222
    case 17 : // space group P2122 / P2212 / P2221 
    case 18 : // space group P22121 / P21221 / P21212 
    case 19 : // space group P212121
      
      c0 = clipper::Coord_frac(0.0, 0.0, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.0, 0.0, 0.5);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.0, 0.5, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.0, 0.5, 0.5);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.5, 0.0, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.5, 0.0, 0.5);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.5, 0.5, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.5, 0.5, 0.5);
      origins.push_back(c0);
      break;
      
    case 20 : // space group C2221
    case 21 : // space group C222
      
      c0 = clipper::Coord_frac(0.0, 0.0, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.0, 0.0, 0.5);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.5, 0.0, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.5, 0.0, 0.5);
      origins.push_back(c0);
      break;
      
    case 22 : // space group F222
      
      c0 = clipper::Coord_frac(0.00, 0.00, 0.00);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.25, 0.25, 0.25);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.50, 0.50, 0.50);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.75, 0.75, 0.75);
      origins.push_back(c0);
      break;
      
      
    case 23 : // space group I222
    case 24 : // space group I212121
      
      c0 = clipper::Coord_frac(0.0, 0.0, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.0, 0.0, 0.5);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.0, 0.5, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.5, 0.0, 0.0);
      origins.push_back(c0);
      break;
      
    case 75 : // space group P4
    case 76 : // space group P41
    case 77 : // space group P42
    case 78 : // space group P43
      
      for ( int iw = 0; iw <= grid.nw(); iw++ )
      {
        c0 = clipper::Coord_frac( 0.0, 0.0, 0.0 + float(iw)/float(grid.nw() + 1)  );
        origins.push_back(c0);
        c0 = clipper::Coord_frac( 0.5, 0.5, 0.0 + float(iw)/float(grid.nw() + 1)  );
        origins.push_back(c0);
      }
      break;
      
      
    case 79 : // space group I4
    case 80 : // space group I41
      
      for ( int iw = 0; iw <= grid.nw(); iw++ )
      {
        c0 = clipper::Coord_frac( 0.0, 0.0, 0.0 + float(iw)/float(grid.nw() + 1)  );
        origins.push_back(c0);
      }
      break;
      
      
    case 89 : // space group P422
    case 90 : // space group P4212
    case 91 : // space group P4122
    case 92 : // space group P41212
    case 93 : // space group P4222
    case 94 : // space group P42212
    case 95 : // space group P4322
    case 96 : // space group P43212
      
      c0 = clipper::Coord_frac(0.0, 0.0, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.0, 0.0, 0.5);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.5, 0.5, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.5, 0.5, 0.5);
      origins.push_back(c0);
      break;
      
    case 97 : // space group I422
    case 98 : // space group I4122
      
      c0 = clipper::Coord_frac(0.0, 0.0, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.0, 0.0, 0.5);
      origins.push_back(c0);
      break;
      
    case 143 : // space group P3
    case 144 : // space group P31
    case 145 : // space group P32
      
      for ( int iw = 0; iw <= grid.nw(); iw++ )
      {
        c0 = clipper::Coord_frac( 0.000000, 0.000000, 0.0 + float(iw)/float(grid.nw() + 1)  );
        origins.push_back(c0);
        c0 = clipper::Coord_frac( 0.333333, 0.666667, 0.0 + float(iw)/float(grid.nw() + 1)  );
        origins.push_back(c0);
        c0 = clipper::Coord_frac( 0.666667, 0.333333, 0.0 + float(iw)/float(grid.nw() + 1)  );
        origins.push_back(c0);
      }
      break;
      
    case 146 : // space group R3:H "H3" (R3:Rhombohedral setting not dealt with here)
      
      for ( int iw = 0; iw <= grid.nw(); iw++ )
      {
        c0 = clipper::Coord_frac( 0.0, 0.0, 0.0 + float(iw)/float(grid.nw() + 1) );
        origins.push_back(c0);
      }
      break;
      
      /*
       // space group  R3: Rhombohedral setting
       for ( int iu = 0; iu <= grid.nu); iu++ )
       {
       c0 = clipper::Coord_frac( 0.0 + float(iu)/float(grid.nu() + 1), 0.0 + float(iu)/float(grid.nu() + 1), 0.0 + float(iu)/float(grid.nu() + 1) );
       origins.push_back(c0);
       }
       */
      
    case 149 : // space group  P312
    case 151 : // space group  P3112
    case 153 : // space group  P3212
      
      c0 = clipper::Coord_frac(0.000000, 0.000000, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.000000, 0.000000, 0.5);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.333333, 0.666667, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.333333, 0.666667, 0.5);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.666667, 0.333333, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.666667, 0.333333, 0.5);
      origins.push_back(c0);
      break;
      
    case 150 : // space group  P321
    case 152 : // space group  P3121
    case 154 : // space group  P3221
      c0 = clipper::Coord_frac(0.0, 0.0, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.0, 0.0, 0.5);
      origins.push_back(c0);
      break;
      
    case 155 : // space group  R32:Hexagonal setting "H32" (R32:Rhombohedral setting not dealt with here)
      
      c0 = clipper::Coord_frac(0.0, 0.0, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.0, 0.0, 0.5);
      origins.push_back(c0);
      break;
      
      /*
       // space group  R32: Rhombohedral setting
       c0 = clipper::Coord_frac(0.0, 0.0, 0.0);
       origins.push_back(c0);
       c0 = clipper::Coord_frac(0.5, 0.5, 0.5);
       origins.push_back(c0);
       */
      
    case 168 : // space group P6
    case 169 : // space group P61
    case 170 : // space group P65
    case 171 : // space group P62
    case 172 : // space group P64
    case 173 : // space group P63
      
      for ( int iw = 0; iw <= grid.nw(); iw++ )
      {
        c0 = clipper::Coord_frac( 0.0, 0.0, 0.0 + float(iw)/float(grid.nw() + 1) );
        origins.push_back(c0);
      }
      break;
      
    case 177 : // space group P622
    case 178 : // space group P6122
    case 179 : // space group P6522
    case 180 : // space group P6222
    case 181 : // space group P6422
    case 182 : // space group P6322
      
      c0 = clipper::Coord_frac(0.0, 0.0, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.0, 0.0, 0.5);
      origins.push_back(c0);
      break;
      
    case 195 : // space group P23
    case 198 : // space group P213
      
      c0 = clipper::Coord_frac(0.0, 0.0, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.5, 0.5, 0.5);
      origins.push_back(c0);
      break;
      
    case 196 : // space group F23
      
      c0 = clipper::Coord_frac(0.00, 0.00, 0.00);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.25, 0.25, 0.25);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.50, 0.50, 0.50);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.75, 0.75, 0.75);
      origins.push_back(c0);
      break;
      
    case 197 : // space group I23
    case 199 : // space group I213
      
      c0 = clipper::Coord_frac(0.0, 0.0, 0.0);
      origins.push_back(c0);
      break;
      
    case 207 : // space group P432
    case 208 : // space group P4232
    case 212 : // space group P4332
    case 213 : // space group P4132
      
      c0 = clipper::Coord_frac(0.0, 0.0, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.5, 0.5, 0.5);
      origins.push_back(c0);
      break;
      
    case 209 : // space group F432
    case 210 : // space group F4132
      
      c0 = clipper::Coord_frac(0.0, 0.0, 0.0);
      origins.push_back(c0);
      c0 = clipper::Coord_frac(0.5, 0.5, 0.5);
      origins.push_back(c0);
      break;
      
    case 211 : // space group I432
    case 214 : // space group I4132
      
      c0 = clipper::Coord_frac(0.0, 0.0, 0.0);
      origins.push_back(c0);
      break;
      
      /*
       // Search for the origin with "sub-voxel" accuracy - on a grid twice as fine as that input.
       
       for ( int iw = 0; iw <= 2*grid.nw()+1; iw++ )
       {
       c0 = clipper::Coord_frac( 0.0, 0.0, 0.0 + float(iw) / (2.0*float(grid.nw()) + 2.0) );
       origins.push_back(c0);
       }
       */
      
    default :
      std::cerr << "\n Uh-Oh - Looks like the origin choices for your space group aren't encoded yet \n"  << std::endl;
      // avert your eyes c++ aficionados
  }
}

void ipa_functions::phased_translation_function (const clipper::HKL_data<clipper::data64::F_phi> &fa, const clipper::HKL_data<clipper::data64::F_phi> &fb, const bool &allow_origin_shifts, clipper::Coord_frac &x, float &cc)
{
  // Computing the Phased Translation function (Determining the translational offset between two crystallographic images with the same symmetry)
  //
  // For discussion see:
  // Read RJ, Schierbeek AJ. A phased translation function. J Appl Cryst. 1988 Oct;21(5):490–5.
  // Bentley GA, Houdusse A. Some applications of the phased translation function in macromolecular structure determination. Acta Crystallogr A. 1992 May 1;48 ( Pt 3):312–22.  '
  // Colman, PM, Fehlhammer, H. & Bartels, K. 1976.  Patterson Search Methods in Protein Structure Determination: Beta_trypsin and Immunoglobulin fragments. Crystallographic Computing Techniques
  // See also the general discussion about "phase correlation" on Wikipedia https://en.wikipedia.org/wiki/Phase_correlation
  // The calculation is performed in Fourier space, so we can make use of the FFT
  //
  // Suppose we have image A with corresponding density rhoA and Structure Factors Fa(h)
  //             and image B with corresponding density rhoB and Structure Factors Fb(h)
  // The real space integral we want to compute is ∫rhoA(x)rhoB(x-t)dx, to determine the translational offset t that maximizes the overlap between the two images
  // This is the mean of the product of the two density functions, evaluated over the unit cell, under the shift t.
  // If we divide that by √(∫rhoA(x)^2dx∫rhoB(x)^2dx) we get the correlation coefficient, which is useful for comparative purposes
  // N.B = this is true so long as we omit the F(000) term from the Fourier summation so that the mean of the syntheses are all zero
  //
  // The integral can be conveniently computed in Fourier space as:
  // ∑Fa(h)Fb*(h)exp(-2πih.t) where Fb*(h) is the complex conjugate of Fb(h).
  // If the two functions have the same symmetry, then the summation can be performed over the reciprocal space asymmetric unit.
  // You just need to compute the complex coefficients Fa(h)Fb*(h)and do a FFT. This is what the Clipper originmatch function does
  //
  // The denominator, required for normalization can be computed by straight summation
  // ∫rhoA(x)^2dx = ∑Fa(h)^2 & ∫rhoB(x)^2dx = ∑Fb(h)^2
  // However if the summation is done over the reciprocal space asymmetric unit, you need to include the appropriate statistical weights epsilon.
  // You then need to multiply appropriately to obtain the sums that would be obtained over the full sphere of data in reciprocal space
  // Oh  ... and the cell volume pops up in there too as a normalizing factor
  
  // The maximum of the function corresponds to the optimum superposition of map A onto Map B, and its coordinates give the translation to be applied to map B
  // The change of basis is defined by:
  // x' = x + t where x is a vector of original coordinate system for map B and t is the translation vector.
  // x' = x - t where x is a vector of original coordinate system for map A and t is the translation vector.
  
  // fa and fb better have the same spacegroup and cell dimensons, or this will fail spectacularly
  
  const clipper::Spacegroup spacegroup_a  = fa.base_hkl_info().spacegroup();
  const clipper::Cell       cell_a        = fa.base_hkl_info().cell();
  const clipper::Resolution resolution_a  = fa.base_hkl_info().resolution();
  
  
  const clipper::Spacegroup spacegroup_b  = fb.base_hkl_info().spacegroup();
  const clipper::Cell       cell_b        = fb.base_hkl_info().cell();
  const clipper::Resolution resolution_b  = fb.base_hkl_info().resolution();
  
  //   clipper::HKL_info hkla(spacegroup_a, cell_a, resolution_a);
  //   hkla.generate_hkl_list(); // initially there are no reflections in the reflection list - better make some
  //   clipper::HKL_info hklb(spacegroup_b, cell_b, resolution_b);
  //   hklb.generate_hkl_list(); // initially there are no reflections in the reflection list - better make some
  
  if ( !(spacegroup_a.spacegroup_number() == spacegroup_b.spacegroup_number() ) or
      !(cell_a.equals(cell_b,1.0)) )
  {
    std::cerr << "Mismatching data sent to the phased translation function" << std::endl; return;
  }
  
  // Okay so the cell and resolution get carried directly from the input, while the space group of the PTF is constructed from the *point group symmetry*  of the input space group
  // This is following Kevin Cowtan - need to understand the logic here.
  
  const clipper::Cell       cell_ptf(cell_a);
  const clipper::Resolution resolution_ptf(resolution_a);
  
  //          const clipper::Spacegroup spacegroup_ptf(spacegroup_a);
  
  const clipper::Spgr_descr Spgr_descr( spacegroup_a.generator_ops().pgrp_ops() );
  const clipper::Spacegroup spacegroup_ptf(Spgr_descr);
  
  // Make a data object to hold the Fourier coefficients for the Phased Translation Function
  
  clipper::HKL_info hkl(spacegroup_ptf, cell_ptf, resolution_ptf);
  hkl.generate_hkl_list(); // initially there are no reflections in the reflection list - better make some
  
  clipper::HKL_data<clipper::data64::F_phi> PTF_coeffs(hkl);
  
  
  // outlog << " Spacegroup " << hkl.spacegroup().symbol_xhm() << std::endl;
  // outlog << cell.format() << std::endl;
  // outlog << "Cell Volume (Angstroms^3) "  << cell.volume() << std::endl;
  // outlog << hkl.resolution().limit() << std::endl;
  
  
  // Compute the Fourier coefficients
  
  //   outlog << "computing Fourier coefficients" << std::endl;
  
  clipper::HKL_data_base::HKL_reference_index ih;
  
  for ( ih = hkl.first(); !ih.last(); ih.next() )
  {
    
    // Have to set the F(0,0,0) terms to zero to achieve zero mean density, and simplify the correlation
    // Also a check for missing data, even though there shouldn't be any in this application
    
    if ( ( (ih.hkl().h() == 0) && (ih.hkl().k() == 0) && (ih.hkl().l() == 0) ) or fa[ih.hkl()].missing()  or fb[ih.hkl()].missing() )
    {
      PTF_coeffs[ih].f()   = 0.0;
      PTF_coeffs[ih].phi() = 0.0;
    }
    else
    {
      PTF_coeffs[ih].f()   = fa[ih.hkl()].f()   * fb[ih.hkl()].f();
      PTF_coeffs[ih].phi() = fa[ih.hkl()].phi() - fb[ih.hkl()].phi();
    }
    
    /*
     if (ih.index()  < 100)
     {
     outlog << ih.index() << " " << ih.hkl().h() << " " << ih.hkl().k() << " " <<  ih.hkl().l() << " "  << fa[ih].f()           <<  " "  << fa[ih].phi()         << std::endl;
     outlog << ih.index() << " " << ih.hkl().h() << " " << ih.hkl().k() << " " <<  ih.hkl().l() << " "  << fb[ih].f()           <<  " "  << fb[ih].phi()         << std::endl;
     outlog << ih.index() << " " << ih.hkl().h() << " " << ih.hkl().k() << " " <<  ih.hkl().l() << " "  << PTF_coeffs[ih].f()   <<  " "  << PTF_coeffs[ih].phi() << std::endl;
     }
     */
    
  }
  
  clipper::Grid_sampling grid(hkl.spacegroup(), hkl.cell(), hkl.resolution(), 2.0);
  clipper::Xmap<float> PTF(hkl.spacegroup(), hkl.cell(), grid);
  PTF.fft_from(PTF_coeffs);
  
  
  // Output the PTF, with a timestamp, for troubleshooting purposes
  
  /*
   clipper::CCP4MAPfile map_out;
   
   time_t t = time(0);   // get time now
   struct tm * now = localtime( & t );
   
   char io_map_filename[64];
   strftime (io_map_filename,sizeof(io_map_filename),"PTF_%H_%M_%S.ccp4",now);
   map_out.open_write( io_map_filename );
   map_out.export_xmap(PTF); // write the PTF
   map_out.close_write();
   
   */
  
  // Now normalize to convert the integral into the conventional correlation coefficient
  
  double sum_a_squared = 0.0;
  double sum_b_squared = 0.0;
  
  
  for ( ih = hkl.first(); !ih.last(); ih.next() )
  {
    if ( !((ih.hkl().h() == 0) && (ih.hkl().k() == 0) && (ih.hkl().l() == 0)) && !fa[ih.hkl()].missing())
    {
      sum_a_squared += std::pow(fa[ih.hkl()].f(), 2) / float(ih.hkl_class().epsilon());
    }
  }
  
  for ( ih = hkl.first(); !ih.last(); ih.next() )
  {
    if ( !((ih.hkl().h() == 0) && (ih.hkl().k() == 0) && (ih.hkl().l() == 0)) && !fb[ih.hkl()].missing())
    {
      sum_b_squared += std::pow(fb[ih.hkl()].f(), 2) / float(ih.hkl_class().epsilon());
    }
  }
  
  /*
   for ( ih = hkl.first(); !ih.last(); ih.next() )
   {
   if ( !((ih.hkl().h() == 0) && (ih.hkl().k() == 0) && (ih.hkl().l() == 0)) && !fa[ih.hkl()].missing() && !fb[ih.hkl()].missing() )
   {
   //     outlog << ih.index() << " " << ih.hkl().h() << " " << ih.hkl().k() << " " <<  ih.hkl().l() <<  std::endl;
   //     outlog << fa[ih].f()  <<  " " << float(ih.hkl_class().epsilon()) << " " << fa[ih].f()/float(ih.hkl_class().epsilon()) << std::endl;
   //     outlog << std::pow(fa[ih].f()/float(ih.hkl_class().epsilon()), 2) << std::endl;
   //     outlog << fb[ih].f()  <<  " " << float(ih.hkl_class().epsilon()) << " " << fb[ih].f()/float(ih.hkl_class().epsilon()) << std::endl;
   //     outlog << std::pow(fb[ih].f()/float(ih.hkl_class().epsilon()), 2) << std::endl;
   
   //      sum_a_squared += std::pow( fa[ih].f() / float(ih.hkl_class().epsilon() ), 2);
   //      sum_b_squared += std::pow( fb[ih].f() / float(ih.hkl_class().epsilon() ), 2);
   
   //         sum_a_squared += std::pow(fa[ih].f(), 2) / float(ih.hkl_class().epsilon());
   //         sum_b_squared += std::pow(fb[ih].f(), 2) / float(ih.hkl_class().epsilon());
   
   sum_a_squared += std::pow(fa[ih.hkl()].f(), 2) / float(ih.hkl_class().epsilon());
   sum_b_squared += std::pow(fb[ih.hkl()].f(), 2) / float(ih.hkl_class().epsilon());
   }
   }
   */
  
  
  //   outlog << "sum_a_squared " << sum_a_squared << std::endl;
  //   outlog << "sum_b_squared " << sum_a_squared << std::endl;
  
  //      sum_a_squared = 2.0*sum_a_squared*(float(n_asu_p1)/float(n_asu));
  //      sum_b_squared = 2.0*sum_b_squared*(float(n_asu_p1)/float(n_asu));
  
  //   sum_a_squared = 2.0*sum_a_squared*float( spacegroup_a.num_symops() );
  //   sum_b_squared = 2.0*sum_b_squared*float( spacegroup_b.num_symops() );
  
  //   outlog << hkl.cell().volume() << std::endl;
  //   outlog << hkl.spacegroup().num_symops() << std::endl;
  
  double weight = hkl.cell().volume() / (2.0*float( hkl.spacegroup().num_symops() )*std::pow(sum_a_squared*sum_b_squared, 0.5));
  //   outlog << " Weight " << weight << std::endl;
  
  // The allowed origin shifts are not arbitrary ... check the possibilities systematically
  // Generate the alternate origins for the target space group or set the origin shift to the zero vector
  
  //   outlog << "Getting origin shifts" << std::endl;
  
  std::vector<clipper::Coord_frac> origins;
  
  if (allow_origin_shifts)
  {
    ipa_functions::compute_alternate_origins(PTF.spacegroup().spacegroup_number(), grid, origins);
  }
  else
  {
    clipper::Coord_frac c0;
    c0 = clipper::Coord_frac(0.0, 0.0, 0.0);
    origins.push_back(c0);
  }
  
  //   outlog << "Figuring optimum shift " << std::endl;
  //   outlog <<  origins.size()  << std::endl;
  
  for ( int i = 0 ; i < origins.size() ; i++ )
  {
    clipper::Xmap_base::Map_reference_coord iy(PTF); // create a coordinate-like reference to the translation function map
    iy.set_coord(origins[i].coord_grid(grid)); // convert the fractional coordinates of the origin to grid coordinates and use that to set the map reference coordinate
    
    if (i == 0)
    {
      cc = PTF[iy];
      x = origins[i];
    }
    else
    {
      if ( PTF[iy] > cc )
      {
        cc = PTF[iy];
        x = origins[i];
      }
    }
    
    //  outlog << i << " " << origins[i].u() << " " << origins[i].v() <<  " " << origins[i].w() << " " << PTF[iy] << " " << cc << std::endl;
    
  }
  
  cc = cc*weight;
}

void ipa_functions::compute_intensity_statistics (const clipper::HKL_data<clipper::data64::F_sigF> &f, const int &res_bins, std::vector<double> &mean_I_acentric, std::vector<double> &mean_I_centric, std::vector<int> &observation_counts_acentric, std::vector<int> &observation_counts_centric, std::ostream &outlog, bool verbose)
{
  // Compute resolution-dependent mean transformed intensities Î = I(hkl)/ε(hkl), where ε are the statistical weighting factors for each observation hkl
  // TODO switch to calculation of an overall <Î>  ... while centric and acentric data have a different *distribution* for Î, they have the same *mean* Î (at least under Wilson statistics)
   
  clipper::HKL_info hkl =  f.base_hkl_info();
  clipper::Resolution_ordinal resord; // Resolution Ordinal - once initiated this returns the approximate fractional ordinal within a dataset (in the range 0...1) for any specified value of 1 / s^2
  resord.init(hkl, 1.0); // Initiating using HKL_info means the ordinal is constructed on the basis of all data. The binning here *must* be the same as the binning used in the function substitute_Fourier_amplitudes
  
  // According to Wilson statistics, the Intensities have a gamma distribution
  // Taking the arithmetic mean as an estimator o f the sample mean is not robust to outliers, but for now let's just get it done.
  // TODO ... Implement a robust estimator of the mean intensity 
  
  clipper::HKL_data_base::HKL_reference_index ih;
  
  for ( ih = f.first(); !ih.last(); ih.next() )
  {
    if ( !f[ih].missing())
    {
      int bin = std::min ( int(static_cast<double>(res_bins) * resord.ordinal(ih.invresolsq())) , res_bins-1 ); // Which bin does this observation belong to ?
      double w = 1.0 / ih.hkl_class().epsilon();  // Get the standard statistical weights for each observation
                                                  //double wa = pw[ih].fom(); // Get the weights associated with the apodization function
      
      //      mean_F[bin]      += w * wa * f[ih].f();
      //      mean_F[res_bins] += w * wa * f[ih].f();
      
      if ( ih.hkl_class().centric() )
      {
        observation_counts_centric[bin]      += 1;
        observation_counts_centric[res_bins] += 1;
        
        mean_I_centric[bin]       += w * pow(f[ih].f(), 2);
        mean_I_centric[res_bins]  += w * pow(f[ih].f(), 2);
        
        //mean_I_centric[bin]       += w * pow(wa * f[ih].f(), 2);
        //mean_I_centric[res_bins]  += w * pow(wa * f[ih].f(), 2);
      }
      else
      {
        observation_counts_acentric[bin]      += 1;
        observation_counts_acentric[res_bins] += 1;
        
        mean_I_acentric[bin]       += w * pow(f[ih].f(), 2);
        mean_I_acentric[res_bins]  += w * pow(f[ih].f(), 2);
        
        //mean_I_acentric[bin]       += w * pow(wa * f[ih].f(), 2);
        //mean_I_acentric[res_bins]  += w * pow(wa * f[ih].f(), 2);
      }
    }
  }
  
  for ( int bin = 0; bin <= res_bins; bin++ )
  {
    if (observation_counts_centric[bin] > 0)
    {
      mean_I_centric[bin] = mean_I_centric[bin] / static_cast<double>(observation_counts_centric[bin]);
    }
    if (observation_counts_acentric[bin] > 0)
    {
      mean_I_acentric[bin] = mean_I_acentric[bin] / static_cast<double>(observation_counts_acentric[bin]);
    }
  }
  
  // now output the intensity statistics, as a function of resolution
  
  resord.invert(); // invert the distribution so that now we specify the ordinal value and it returns 1/|s|**2
  
  double lower_limit;
  double upper_limit;
  
  
  
  for ( int bin = 0; bin <= res_bins; bin++ )
  {
    if (bin == 0)
    {
      lower_limit = 0.0;
      upper_limit = sqrt( resord.ordinal( static_cast<double>(bin+1)/static_cast<double>(res_bins) ) );
    }
    else if ((bin > 0) && (bin < res_bins))
    {
      lower_limit = sqrt( resord.ordinal( static_cast<double>(bin)/static_cast<double>(res_bins) ) );
      upper_limit = sqrt( resord.ordinal( static_cast<double>(bin+1)/static_cast<double>(res_bins) ) );
    }
    else
    {
      lower_limit = 0.0;
      upper_limit = sqrt(resord.ordinal( 1.0 ));
    }
    if (verbose)
    {
      outlog << std::fixed;
      if (bin < res_bins)
      { outlog << std::setw(3) << bin  << " ";}
      else
      { outlog << "    ";}
      outlog <<
      std::setw(8) << lower_limit << " " <<
      std::setw(8) << upper_limit << " " <<
      std::setw(8) << 1/upper_limit <<  " "  <<
      std::setw(8) << observation_counts_acentric[bin] <<  " "  <<
      std::scientific <<
      std::setw(10) << mean_I_acentric[bin] <<  " "  <<
      std::fixed <<
      std::setw(8) << observation_counts_centric[bin] <<  " "  <<
      std::scientific <<
      std::setw(10) << mean_I_centric[bin] << std::endl;
    }
  }
  
}

void ipa_functions::compute_phase_agreement (const clipper::HKL_data<clipper::data64::F_phi> &fp1, clipper::HKL_data<clipper::data64::F_phi> fp2, const clipper::HKL_data<clipper::data64::Phi_fom> &pw,
                                             const int &res_bins, const bool &allow_origin_shifts, const bool &invert, float &diff_all, float &diff_centric, float &diff_acentric, float &map_cc, clipper::Coord_frac &origin_shift, std::ostream &outlog, bool verbose)
{
  
  if (true)
  {
    // Here we're going to do a hash check to ensure that the inputs are valid across different things:
    clipper::HKL_data_base::HKL_reference_index ih;
    
    float hash_fp1 = 0.0;
    float hash_fp2 = 0.0;
    float hash_pw = 0.0;
    
    for (ih = fp1.first(); !ih.last() ; ih.next() )
    {
      if (!fp1[ih].missing())
      {
        hash_fp1 += fp1[ih].phi();
      }
      if (!fp2[ih.hkl()].missing())
      {
        hash_fp2 += fp2[ih].phi();
      }
      if (!pw[ih.hkl()].missing())
      {
        hash_pw += pw[ih].phi();
      }
    }
    
    outlog << "Hash: " << std::setprecision(4) << std::setw(10) << hash_fp1 << ":" << hash_fp2 << ":" << hash_pw << " : " << res_bins << " : " << allow_origin_shifts << std::endl;
  }
  
  clipper::HKL_data_base::HKL_reference_index ih;
  
  clipper::HKL_info hkl =  fp1.base_hkl_info();
  clipper::Resolution_ordinal resord; // Resolution Ordinal - once initiated this returns the approximate fractional ordinal within a dataset (in the range 0...1) for any specified value of 1 / s^2
  resord.init( hkl, 1.0 );
  
  // initialize vectors with the appropriate number of elements, all set to zero
  // bins 0,1, ...,n-1  store the individual bin statistics. bin n stores statistics for the entire data set.
  
  std::vector<double> observation_counts_all(res_bins + 1, 0.0);
  std::vector<double> phase_differences_all(res_bins + 1, 0.0);
  std::vector<double> weights_all(res_bins + 1, 0.0);
  
  std::vector<double> observation_counts_centric(res_bins + 1, 0.0);
  std::vector<double> phase_differences_centric(res_bins + 1, 0.0);
  std::vector<double> weights_centric(res_bins + 1, 0.0);
  
  std::vector<double> observation_counts_acentric(res_bins + 1, 0.0);
  std::vector<double> phase_differences_acentric(res_bins + 1, 0.0);
  std::vector<double> weights_acentric(res_bins + 1, 0.0);
  
  // First phase set is passed by reference. 2nd phase set is passed by value since we want to operate on it, and not mess with the original
  // Is this bad practice ? I leave this to the programming aficianados
  
  if (invert) // Invert the phases as requested before comparison. Clipper should take care of the tricky origin of inversion issues that exist in some space groups
  {
    for ( ih = fp2.first(); !ih.last(); ih.next() ) { fp2[ih].friedel(); }
  }
  
  if (allow_origin_shifts) // Use the FFT based image registration to determine the origin shift. This shift is applied to the second set of phases (fp2) to match the first (fp1)
  {
    float cc=0.0;
    ipa_functions::phased_translation_function(fp1, fp2, allow_origin_shifts, origin_shift, cc);
    if (verbose) {outlog << "Origin shift to maximize phase agreement: " << origin_shift.format() << std::endl;}
  }
  
  // Quantities required for calculation of the map correlation coefficient
  double cc_num = 0.0;
  double cc_denom1 = 0.0;
  double cc_denom2 = 0.0;
  
  // now loop over the asymmetric unit and perform the summations required to calculate the mean phase differences and the map correlation coefficient.
  for ( ih = fp1.first(); !ih.last(); ih.next() )
  {
    if ( !fp1[ih].missing() && !fp2[ih.hkl()].missing() )
    {
      const int bin = std::min ( int(static_cast<double>(res_bins) * resord.ordinal(ih.invresolsq())) , res_bins-1 ); // Which bin does this observation belong to ?
                                                                                                                      //const int bin = std::min ( int(double(res_bins) * resord.ordinal(ih.invresolsq())) , res_bins-1 ); // Which bin does this observation belong to ?
      const double w = 1.0 / ih.hkl_class().epsilon();  // Get the standard statistical weights
      const double wa = pw[ih].fom(); // Get the weights associated with the apodization function, used for computing the weighted statistics
      
      //  the absolute phase difference, accounting for periodicity - no phase shift applied
      //  double d = acos ( cos( fp2[ih].phi() - fp2[ih].phi() ) );
      // For computation of phase shift from the shift in origin, see Giacovazzo (2014) Phasing in Crystallography Eq 3.3
      
      clipper::Coord_reci_frac h( ih.hkl() ); // fractional coordinates
      const double hx = h * origin_shift;
      const double delta_phi = clipper::Util::twopi() * ( hx - floor(hx) );
      const double c_phase_diff = cos( fp1[ih].phi() - (fp2[ih.hkl()].phi() + delta_phi) );
      const double d = acos ( c_phase_diff );
      
      const double f1 = fp1[ih].f();
      const double f2 = fp2[ih.hkl()].f();
      
      // Augment the sums required for computation of the map correlation coefficient.
      // For details of the correlation coefficient calculation see 
      // Lunin VY, Woolfson MM. Mean phase error and the map-correlation coefficient. Acta Crystallogr D. 1993 Nov 1;49(Pt 6):530–3. 
      // also Bailey GD, Hyun J-KK, Mitra AK, Kingston RL. J Mol Biol. 2012 Mar 30;417(3):212–23. Supplementary Methods, Section 3.
      //
      // Note that the statistical weights must be included here as we are summing over the asymmetric unit only.
      
      if ( !((ih.hkl().h() == 0) && (ih.hkl().k() == 0) && (ih.hkl().l() == 0)) ) // F(000) must be excluded from these summations to ensure mean densities are zero
      {
        // outlog << ih.hkl().h() << " " << ih.hkl().k() << " " << ih.hkl().l() << " " << w << " " << f1 <<  " " << f2 <<  " " << c_phase_diff << std::endl;
        cc_num    += w*f1 * w*f2 * c_phase_diff;
        cc_denom1 += w*f1 * w*f1;
        cc_denom2 += w*f2 * w*f2;
      }    
      
      // Augment the sums required for computation of weighted mean absolute phase difference 
      
      observation_counts_all[bin] += 1.0;
      phase_differences_all[bin]  += w * wa * d;
      weights_all[bin]            += w * wa;
      
      observation_counts_all[res_bins] += 1.0;
      phase_differences_all[res_bins]  += w * wa * d;
      weights_all[res_bins]            += w * wa;
      
      if ( ih.hkl_class().centric() )
      {
        observation_counts_centric[bin] += 1.0;
        phase_differences_centric[bin]  += w * wa * d;
        weights_centric[bin]            += w * wa;
        
        observation_counts_centric[res_bins] += 1.0;
        phase_differences_centric[res_bins]  += w * wa * d;
        weights_centric[res_bins]            += w * wa;
      }
      else
      {
        observation_counts_acentric[bin] += 1.0;
        phase_differences_acentric[bin]  += w * wa * d;
        weights_acentric[bin]            += w * wa;
        
        observation_counts_acentric[res_bins] += 1.0;
        phase_differences_acentric[res_bins]  += w * wa * d;
        weights_acentric[res_bins]            += w * wa;
      }
    }
  }
  
  // Compute the map correlation coefficient
  
  map_cc = cc_num / std::sqrt(cc_denom1*cc_denom2);
  
  // Compute the weighted mean absolute phase differences
  for ( int bin = 0; bin <= res_bins; bin++ )
  {
    phase_differences_all[bin] /= std::max( weights_all[bin], 1.0 );
    phase_differences_centric[bin] /= std::max( weights_centric[bin], 1.0 );
    phase_differences_acentric[bin] /= std::max( weights_acentric[bin], 1.0 );
  }
  
  resord.invert(); // invert the distribution so that now we specify the ordinal value and it returns 1/|s|**2
  
  double lower_limit;
  double upper_limit;
  
  // report statistics as a function of resolution
  if (verbose)
  {
    outlog << "\nWeighted mean phase difference as a function of resolution" << std::endl;
    outlog <<   "----------------------------------------------------------" << std::endl;
    
    outlog << "                                     all data                centric data              acentric data    " << std::endl;
    outlog << " |s| lo   |s| hi  1/|s| hi       #   mean w    Δphi        #   mean w    Δphi        #   mean w    Δphi " << std::endl;
    
    for ( int bin = 0; bin < res_bins; bin++ )
    {
      if (bin == 0)
      {
        lower_limit = 0.0;
        upper_limit = sqrt( resord.ordinal( double(bin+1)/double(res_bins) ) );
      }
      else
      {
        lower_limit = sqrt( resord.ordinal( double(bin)/double(res_bins) ) );
        upper_limit = sqrt( resord.ordinal( double(bin+1)/double(res_bins) ) );
      }
      outlog << std::fixed <<
      std::setw(8) << lower_limit << " " <<
      std::setw(8) << upper_limit << " " <<
      std::setw(8) << 1/upper_limit <<  " "  <<
      std::setw(8) << int(observation_counts_all[bin]) <<  " "  <<
      std::setw(8) << weights_all[bin]/observation_counts_all[bin] <<  " "  <<
      std::setw(8) << clipper::Util::rad2d(phase_differences_all[bin]) <<
      std::setw(8) << int(observation_counts_centric[bin]) <<  " "  <<
      std::setw(8) << weights_centric[bin]/observation_counts_centric[bin] <<  " "  <<
      std::setw(8) << clipper::Util::rad2d(phase_differences_centric[bin]) <<
      std::setw(8) << int(observation_counts_acentric[bin]) <<  " "  <<
      std::setw(8) << weights_acentric[bin]/observation_counts_acentric[bin] <<  " "  <<
      std::setw(8) << clipper::Util::rad2d(phase_differences_acentric[bin]) << std::endl;
    }
    
    // report the overall statistics
    
    lower_limit = 0.0;
    upper_limit = sqrt(resord.ordinal( 1.0 ));
    
    outlog << std::fixed <<
    std::setw(8) << lower_limit << " " <<
    std::setw(8) << upper_limit << " " <<
    std::setw(8) << 1/upper_limit <<  " "  <<
    std::setw(8) << int(observation_counts_all[res_bins]) <<  " "  <<
    std::setw(8) << weights_all[res_bins]/observation_counts_all[res_bins] <<  " "  <<
    std::setw(8) << clipper::Util::rad2d(phase_differences_all[res_bins]) <<
    std::setw(8) << int(observation_counts_centric[res_bins]) <<  " "  <<
    std::setw(8) << weights_centric[res_bins]/observation_counts_centric[res_bins] <<  " "  <<
    std::setw(8) << clipper::Util::rad2d(phase_differences_centric[res_bins]) <<
    std::setw(8) << int(observation_counts_acentric[res_bins]) <<  " "  <<
    std::setw(8) << weights_acentric[res_bins]/observation_counts_acentric[res_bins] <<  " "  <<
    std::setw(8) << clipper::Util::rad2d(phase_differences_acentric[res_bins]) << std::endl;

    outlog << "\nMap Correlation Coefficient: " << map_cc << std::endl;
  }
  
  // Return the overall mean absolute phase differences in degrees 
  diff_all      = float(clipper::Util::rad2d(phase_differences_all[res_bins]));
  diff_centric  = float(clipper::Util::rad2d(phase_differences_centric[res_bins]));
  diff_acentric = float(clipper::Util::rad2d(phase_differences_acentric[res_bins]));

} 

void ipa_functions::smooth_filter_map_region(clipper::Xmap<float>& map_to_smooth, // Input map to be smoothed.
                              const clipper::Xmap<int>& mask, // An integer-valued mask, with the same dimensions as the map, which defines the region to be modified. 
                              const double& smoothfilter_radius, 
                              double& smoothed_map_mean,
                              const int affected_region_flag,
                              bool ignore_unaffected_region)
{
  /*
  This function is an alternative to solvent flattening, however can be applied to any region to an associated map, and works as a general smoothing function.
  Instead of flattening the solvent to a value, we apply a smoothing filter to the region determined by the int map, and provided solvent region flag. The smoothing is
  performed in fourier space, as it is more efficient generally.
  */

  clipper::Xmap<float>     smoothed_map( map_to_smooth.spacegroup(), map_to_smooth.cell(), map_to_smooth.grid_sampling() ); // Declare a temperary working map.
  MapFilterFn_triweight filter_function(smoothfilter_radius); // Construct a triweight filter function
  clipper::MapFilter_fft<float> map_filter( filter_function, 1.0, clipper::MapFilter_fft<float>::Relative ); // Construct the map filter

  clipper::Xmap<float> container_map; // Declare a temperary working map, may be needed later.


  if (ignore_unaffected_region)
  {
    // We minise the edge effects of the other region, by setting the values within that region to the mean of the affected region.
    int count;
    double cusum;
    
    container_map.init( map_to_smooth.spacegroup(), map_to_smooth.cell(), map_to_smooth.grid_sampling() ); // Init this now.

    // First copy the original non-affected parts, they will be needed later, and also calculate the mean of the solvent region.
    clipper::Xmap_base::Map_reference_index ix;
    for (ix = map_to_smooth.first(); !ix.last(); ix.next() )
    {
      if (mask[ix] == affected_region_flag)
      {
        // needed to calc mean.
        cusum += map_to_smooth[ix];
        count++;
      }
      container_map[ix] = map_to_smooth[ix]; // Just store the original map, unmodified.
    }

    double region_mean_temp = cusum / count; // calculate the mean of the affected region.

    // Set the unaffected region to the mean of the affected region.
    for (ix = map_to_smooth.first(); !ix.last(); ix.next() )
    {
      if (mask[ix] != affected_region_flag)
      {
          map_to_smooth[ix] = region_mean_temp;
      }
    }
  }

  // Smooth the map. Note, this WILL have edge effects of increasingness around the protein region.
  // We could theoretically try minimising this by replacing the unaffected region with the mean of the affected region.
  map_filter (smoothed_map, map_to_smooth); 

  int count = 0;
  double cusum = 0;

  // Replace the map with its smoothed counterpart for the affected region.
  clipper::Xmap_base::Map_reference_index ix;
  for (ix = map_to_smooth.first(); !ix.last(); ix.next() )
  {
    if (mask[ix] == affected_region_flag)
    {
      map_to_smooth[ix] = smoothed_map[ix];

      count++;
      cusum += map_to_smooth[ix];
    }
    if (ignore_unaffected_region)
    {
      // We need to copy back our unaffected region in this case:
      map_to_smooth[ix] = container_map[ix];
    }
  }

  smoothed_map_mean = cusum / count; // Output the global mean of the resulting maps smoothed region, this is useful if the map requires further scaling to an appropriate mean.
}

void ipa_functions::constant_bias_scale_map_region(clipper::Xmap<float>& map, // Input map 
                                    const clipper::Xmap<int>& mask, // An integer-valued mask, with the same dimensions as the map, which defines the region to be modified. 
                                    double& new_map_mean,
                                    const double& constant,
                                    const double& scale, 
                                    const int mask_flag)
{
  /*
  This function applies a "Constant Bias Scale" to the selected region of a map .. i.e. adds a constant to the map values, and then scales the result 
  With an input scale of 1, it becomes identical to shift_map_region() 
  */

  int count = 0;
  double cusum = 0;

  // Replace the map elements with their constant bias scale counterpart in the affected region.
  clipper::Xmap_base::Map_reference_index ix;
  for (ix = map.first(); !ix.last(); ix.next() )
  {
    if (mask[ix] == mask_flag)
    {
      map[ix] = (map[ix]+constant) * scale;

      count++;
      cusum += map[ix];
    }
  }

  new_map_mean = cusum / count; // Output the new mean of the region that has been operated on 
}

void ipa_functions::power_multiply_map_region(clipper::Xmap<float>& map, // Input map 
                                    const clipper::Xmap<int>& mask, // An integer-valued mask, with the same dimensions as the map, which defines the region to be modified. 
                                    double& new_map_mean,
                                    const double& power, // Power of
                                    const double& bias, // Multiplier
                                    const int mask_flag)
{
  // This function "Power multiplies" the selected region of a map .. i.e. raises the map values to a specified power, then scales the result. 
  

  int count = 0;
  double cusum = 0;

  // Replace the map with its power multiplied counterpart for the affected region.
  clipper::Xmap_base::Map_reference_index ix;
  for (ix = map.first(); !ix.last(); ix.next() )
  {
    if (mask[ix] == mask_flag)
    {
      map[ix] = std::pow(map[ix], power) * bias;

      count++;
      cusum += map[ix];
    }
  }

  new_map_mean = cusum / count; // Output the new mean of the region that has been operated on 

}



void ipa_functions::calculate_mask (const bool &compute_mask_from_local_mean, const bool &compute_mask_from_local_variance, const bool &compute_mask_from_local_mean_and_variance,
                                    const double &filter_radius, const float &mask_weighting_factor_mean,
                                    const clipper::Xmap<float> &map, const float &fraction_solvent, clipper::Xmap<int> &mask, std::ostream &outlog, bool verbose)
{
  if (verbose)
  {
    outlog << "\nCommencing mask calculation with filter radius: " << filter_radius << std::endl;
    outlog << "and solvent fraction: " << fraction_solvent << std::endl;
  }
  clipper::Xmap_base::Map_reference_index ix;
  
  // compute the local mean and local variance of the density , as required
  
  clipper::Xmap<float>     local_mean( map.spacegroup(), map.cell(), map.grid_sampling() );  // create map object to store the local mean
  clipper::Xmap<float> local_variance( map.spacegroup(), map.cell(), map.grid_sampling() );  // create map object to store the local variance
  
  //      clipper::MapFilterFn_linear      filter_function( filter_radius );            // Construct a triangular filter (weight) function
  //      clipper::MapFilterFn_step        filter_function( filter_radius );            // Construct a uniform filter function
  //      clipper::MapFilterFn_quadratic   filter_function( filter_radius );            // Construct a quadratic filter function

  MapFilterFn_triweight   filter_function( filter_radius );            // Construct a triweight filter function
  
  clipper::MapFilter_fft<float> map_filter( filter_function, 1.0, clipper::MapFilter_fft<float>::Relative ); // Construct the map filter
  
  //        outlog << "\nCalculating the local mean of the density with filter radius: " << filter_radius << std::endl;
  
  map_filter (local_mean, map);
  
  if ( (compute_mask_from_local_variance) || (compute_mask_from_local_mean_and_variance) ) // we need in addition the local variance
  {
    //            outlog << "\nCalculating the local variance of the density with filter radius: " << filter_radius << std::endl;
    
    clipper::Xmap<float>    map_squared( map.spacegroup(), map.cell(), map.grid_sampling() );  // create map object to store the estimate for the squared density
    
    for ( ix = map_squared.first(); !ix.last(); ix.next() )
    {
      map_squared[ix] = pow( map[ix], 2.0 );
    }
    
    // Calculate the local mean squared density ...
    
    map_filter(local_variance, map_squared);
    
    // ... and compute the local variance
    
    for ( ix = local_mean.first(); !ix.last(); ix.next() )
    {
      local_variance[ix] = local_variance[ix] - pow( local_mean[ix], 2.0 );
    }
    
  }
  
  // Now calculate mask by thresholding
  
  if (compute_mask_from_local_mean)
  {
    if (verbose)
    {
      outlog << "\nComputing mask from the local mean" << std::endl;
    }
    // In order to determine the density threshold which defines the protein & solvent regions we need to find the appropriate quantile of the map array
    // We will do this exactly via a simple (but costly) sort operation
    
    // TODO: ... Handle multiplicity correctly during this thresholdong step
    // Points in the map array may have non-unitary multiplicity which is not accounted for when creating the sort array, where every point is inserted once irrespective of multiplicity. 
    // However in our defence, the resulting errors appear likely to be inconsequential.
    
    std::vector<float> sort_array;
    
    for ( ix = local_mean.first(); !ix.last(); ix.next() )
    {
      sort_array.push_back(local_mean[ix]); // transfer the data from the map object into the sort array
    }
    
    std::sort(sort_array.begin(), sort_array.end());
    
    int n_asu = sort_array.size();
    
    double threshold = sort_array[ lround(double(n_asu)*fraction_solvent) - 1 ]; // calculate the threshold density value. Since the sort is in ascending order, all points below this threshold are designated as solvent
    if (verbose)
    {
      outlog << "Threshold density value for solvent mask generation: " << threshold << std::endl;
    }
    // create mask by thresholding the map
    
    for ( ix = local_mean.first(); !ix.last(); ix.next() )
    {
      mask[ix] = ( (local_mean[ix] < threshold) ? 0 : 1 ); // use the C++ conditional operator to assign solvent as 0.0 and protein as 1.0
    }
    
  }
  
  
  if (compute_mask_from_local_variance)
  {
    if (verbose)
    {
      outlog << "\nComputing mask from the local variance" << std::endl;
    }
    // In order to determine the density threshold which defines the protein & solvent regions we need to find the appropriate quantile of the map array
    // We will do this exactly via a simple (but costly) sort operation
    // Note that the multiplicity is not handled correctly here - points in the map array may have non-unitary multiplicity and this is not accounted for when creating the sort array, where every point is inserted once. However the resulting errors appear likely to be slight.
    
    std::vector<float> sort_array;
    
    for ( ix = local_variance.first(); !ix.last(); ix.next() )
    {
      sort_array.push_back(local_variance[ix]); // transfer the data from the map object into the sort array
    }
    
    std::sort(sort_array.begin(), sort_array.end());
    
    int n_asu = sort_array.size();
    
    double threshold = sort_array[ lround(double(n_asu)*fraction_solvent) - 1 ]; // calculate the threshold variance value. Since the sort is in ascending order, all points below this threshold are designated as solvent
    
    if (verbose)
    {
      outlog << "Threshold variance for solvent mask generation: " << threshold << std::endl;
    }
    // create mask by thresholding the map
    
    for ( ix = local_variance.first(); !ix.last(); ix.next() )
    {
      mask[ix] = ( (local_variance[ix] < threshold) ? 0 : 1 ); // use the C++ conditional operator to assign solvent as 0.0 and protein as 1.0
    }
    
    // Finally while we are here, compute a histogram from the local variance map and output it, as it's diagnostic of a solution
    // Currently Disabled
    /*
     outlog << "\nHistogram of the Local Variance map\n" << std::endl;
     
     const int n_local_variance_bins = 220;
     histogram local_variance_histogram;
     local_variance_histogram.number_of_intervals = n_local_variance_bins;
     local_variance_histogram.interval = new float[local_variance_histogram.number_of_intervals];
     
     nx = 0;
     for ( ix = local_variance.first(); !ix.last(); ix.next() )
     {
     map_f[nx] = local_variance[ix];
     nx ++;
     }
     
     calculate_moments = false;
     calculate_statistics(&all_regions_flag, &calculate_moments, &n_asu, &map_f[0], &msk_f[0], &weight_f[0], &ed_mean, &ed_variance, &ed_skewness, &ed_kurtosis, &ed_minimum, &ed_maximum);
     local_variance_histogram.minimum = ed_minimum;
     local_variance_histogram.maximum = ed_maximum;
     
     calc_1d_histogram(&all_regions_flag, &n_asu, &map_f[0], &msk_f[0], &weight_f[0], &local_variance_histogram.number_of_intervals, &local_variance_histogram.minimum, &local_variance_histogram.maximum, &local_variance_histogram.interval[0]); // compute histogram of local variance map
     
     delete [] local_variance_histogram.interval;
     */
    
  }
  
  if (compute_mask_from_local_mean_and_variance)
  {
    if (verbose)
    {
      outlog << "\nComputing mask from the local mean and variance with weighting factor: " << mask_weighting_factor_mean << std::endl;
    }
    // To compute a mask based on both the local mean *and* the local variance we need to perform some feature scaling (data normalization).
    // First rescale the local mean and local variance maps, so they have zero mean and unit variance.
    // Add the two maps, weighted appropriately, and then threshold them based on the resultant map to produce a binary mask.
    
    // two pass computation of the mean and variance
    
    // Although we are calculating statistics over the asu I don't think we need to worry about the varying multiplicity of the map points here, as the counting errors will be self-cancelling
    
    
    int n_asu = 0;
    float mean_lm;
    float mean_lv;
    
    for ( ix = map.first(); !ix.last(); ix.next() )
    {
      n_asu += 1;
      mean_lm += (local_mean[ix] );
      mean_lv += (local_variance[ix] );
    }
    mean_lm /= float(n_asu);
    mean_lv /= float(n_asu);
    
    //          int n_asu = 0;
    float var_lm;
    float var_lv;
    
    for ( ix = map.first(); !ix.last(); ix.next() )
    {
      //            n_asu += 1;
      var_lm += pow( (local_mean[ix]     - mean_lm), 2.0);
      var_lv += pow( (local_variance[ix] - mean_lv), 2.0);
    }
    
    var_lm /= float(n_asu-1);
    var_lv /= float(n_asu-1);
    
    float sd_lm = sqrt(var_lm);
    float sd_lv = sqrt(var_lv);
    
    // renormalize
    
    for ( ix = map.first(); !ix.last(); ix.next() )
    {
      local_mean[ix] = (local_mean[ix]         - mean_lm) / sd_lm;
      local_variance[ix] = (local_variance[ix] - mean_lv) / sd_lv;
    }
    
    // construct the "combined" map and put it into local_mean
    
    for ( ix = map.first(); !ix.last(); ix.next() )
    {
      local_mean[ix] = mask_weighting_factor_mean*local_mean[ix] + (1.0-mask_weighting_factor_mean)*local_variance[ix];
    }
    
    std::vector<float> sort_array;
    
    for ( ix = local_mean.first(); !ix.last(); ix.next() )
    {
      sort_array.push_back(local_mean[ix]); // transfer the data from the map object into the sort array
    }
    
    std::sort(sort_array.begin(), sort_array.end());
    
    n_asu = sort_array.size();
    
    double threshold = sort_array[ lround(double(n_asu)*fraction_solvent) - 1 ]; // calculate the threshold variance value. Since the sort is in ascending order, all points below this threshold are designated as solvent
    if (verbose)
    {
      outlog << "Threshold for solvent mask generation: " << threshold << std::endl;
    }
    // create mask by thresholding the map
    
    for ( ix = local_mean.first(); !ix.last(); ix.next() )
    {
      mask[ix] = ( (local_mean[ix] < threshold) ? 0 : 1 ); // use the C++ conditional operator to assign solvent as 0.0 and protein as 1.0
    }
    
    /*
     float mask_weighting_factor_variance = 1.0 - mask_weighting_factor_mean;
     int n_asu = 0;
     for ( ix = map.first(); !ix.last(); ix.next() )
     {
     n_asu += 1;
     }
     
     std::pair<float,int> *rank_mean;
     rank_mean = new std::pair<float,int>[n_asu];
     std::pair<float,int> *rank_variance;
     rank_variance = new std::pair<float,int>[n_asu];
     
     float *combined_rank;
     combined_rank = new float[n_asu];
     
     
     for ( ix = map.first(); !ix.last(); ix.next() )
     {
     rank_mean[ix.index()] = std::make_pair(local_mean[ix], ix.index());
     rank_mean[ix.index()] = std::make_pair(local_variance[ix], ix.index());
     }
     
     //   sort the values. note that by default std::sort prioritizes the first element of a pair, which is what we want.
     
     std::sort(rank_mean,     rank_mean     + n_asu);
     std::sort(rank_variance, rank_variance + n_asu);
     
     for ( int i = 0; i<n_asu; i++ )
     {
     combined_rank[rank_mean.second]     +=     mask_weighting_factor_mean*float(i);
     combined_rank[rank_variance.second] += mask_weighting_factor_variance*float(i);
     }
     
     
     delete [] rank_mean;
     delete [] rank_variance;
     delete [] combined_rank;
     */
    
  }
  
}

void ipa_functions::enforce_mask_connectivity (const bool &erase_islands, const bool &fill_voids, const float &erase_islands_threshold, const float &fill_voids_threshold,
                                               clipper::Xmap<int> &mask,
                                               int &n_connected_sets_eliminated_protein, int &n_connected_sets_retained_protein, int &n_connected_sets_eliminated_solvent, int &n_connected_sets_retained_solvent, std::ostream &outlog,
                                               bool verbose)
{
  
  // First expand the mask into space group P1
  
  const clipper::Spacegroup space_group_p1( clipper::Spacegroup::P1 );
  if (verbose)
  {
    outlog << "\nExpanding mask into space group P1" << std::endl;
  }
  clipper::Xmap<int> expanded_mask(space_group_p1, mask.cell(), mask.grid_sampling() ); // define a map with the same cell and grid as the target mask, but space group symmetry P1
  
  clipper::Xmap_base::Map_reference_index ix;
  
  for ( ix = expanded_mask.first(); !ix.last(); ix.next() )
  {
    //             clipper::Xmap_base::Map_reference_coord iy(mask);
    //             iy.set_coord( ix.coord() );
    //             expanded_mask[ ix ] = mask[ iy ];
    expanded_mask[ ix ] = mask.get_data( ix.coord() );
  }
  
  int n_connected_sets_eliminated;
  int n_connected_sets_retained;
  
  // erase islands
  
  if (erase_islands)
  {
    if (verbose)
    {
      outlog << "\nConnectivity reinforcement - will erase the smallest connected sets from the protein mask" << std::endl;
    }
    // Now compute the connected sets of points, and edit the mask
    int mask_id = 1;
    bool apply_translational_symm = true;
    bool check_diagonals = false;
    
    ipa_functions::mask_connections(mask_id, apply_translational_symm, check_diagonals, erase_islands_threshold, expanded_mask, n_connected_sets_eliminated_protein, n_connected_sets_retained_protein, outlog, verbose);
  }
  
  // fill voids
  
  if (fill_voids)
  {
    
    // complement the mask
    
	// TODO replace this with a call to the library function
    for ( ix = expanded_mask.first(); !ix.last(); ix.next() )
    {
      int current_point = expanded_mask[ ix ];
      expanded_mask[ ix ] = ( (current_point == 0) ? 1 : 0 );
    }
    
    if (verbose)
    {
      outlog << "\nConnectivity reinforcement - will erase the smallest connected sets from the solvent mask" << std::endl;
    }
    // Now compute the connected sets of points, and edit the mask
    int mask_id = 1;
    bool apply_translational_symm = true;
    bool check_diagonals = false;
    
    ipa_functions::mask_connections(mask_id, apply_translational_symm, check_diagonals, fill_voids_threshold, expanded_mask, n_connected_sets_eliminated_solvent, n_connected_sets_retained_solvent, outlog, verbose);
    
    // complement the mask
	// TODO replace this with a call to the library function
    
    for ( ix = expanded_mask.first(); !ix.last(); ix.next() )
    {
      int current_point = expanded_mask[ ix ];
      expanded_mask[ ix ] = ( (current_point == 0) ? 1 : 0 );
    }
    
  }
  
  // Now put the contents of the edited mask back into the orginal map object
  if (verbose)
  {
    outlog << "\nContracting mask to original space group\n" << std::endl;
  }
  for ( ix = mask.first(); !ix.last(); ix.next() )
  {
    //             clipper::Xmap_base::Map_reference_coord iy(expanded_mask);
    //             iy.set_coord( ix.coord() );
    //             mask[ix] = expanded_mask[iy];
    mask[ix] = expanded_mask.get_data(ix.coord());
  }
  
}

void ipa_functions::compute_mask_agreement (const clipper::Xmap<int> &mask_a, const clipper::Xmap<int> &mask_b, clipper::Coord_frac &origin_shift, const bool &compute_real_space_agreement, float &mask_cc, float &simple_matching_coefficient, float &cohen_kappa, float &mcc, std::ostream &outlog, bool verbose)
{
  // The function returns the correlation coefficient between two integer-valued binary masks, computed in Fourier space using a "phased translation function".
  // We assume we are comparing two crystallographic masks generated in the same space group, with the same cell dimensions, the same grid, and the same integer representation
  // However the two masks may be origin shifted, as permitted by the space group symmetry. The function returns the origin shift.
  // If compute_real_space_agreement is true, the function also returns the simple matching coefficient, Cohen's kappa coefficient, and the Matthews correlation coefficient between the masks, after application of the calculated origin shift.
  // If compute_real_space_agreement is false, the function returns only the correlation coefficient, computed in Fourier space, and the associated origin shift
  
  // Note that the Matthews correlation coefficient and Cohen's kappa become equivalent when the masks have exactly the same solvent volume fraction.
  // When the volume fractions are not the same the quantities are distinct with Cohen's kappa always <= correlation coefficient
  
  clipper::Xmap_base::Map_reference_index ix;
  
  // use a fft based method (phased translation function) to get the required origin shift
  
  clipper::Resolution mask_res(3.0);
  
  clipper::HKL_info hkl_a(mask_a.spacegroup(), mask_a.cell(), mask_res);
  hkl_a.generate_hkl_list(); // initially there are no reflections in the reflection list - better make some
  clipper::HKL_data<clipper::data64::F_phi> f_a (hkl_a); // create a data object of type  F + phi to hold Fourier coefficients for the first mask
  clipper::Xmap<float> mask_a_real(mask_a.spacegroup(), mask_a.cell(), mask_a.grid_sampling()); // floating point map object associated with the first mask
  for ( ix = mask_a.first(); !ix.last(); ix.next() )
  {
    mask_a_real[ix] = (float)mask_a[ix];
  }
  mask_a_real.fft_to (f_a); // FFT of the real-valued mask to get the Crystallographic structure factors
  
  
  clipper::HKL_info hkl_b(mask_b.spacegroup(), mask_b.cell(), mask_res);
  hkl_b.generate_hkl_list(); // initially there are no reflections in the reflection list - better make some
  clipper::HKL_data<clipper::data64::F_phi> f_b ( hkl_b ); // create a data object of type  F + phi to hold Fourier coefficients for the second mask
  clipper::Xmap<float> mask_b_real(mask_b.spacegroup(), mask_b.cell(), mask_b.grid_sampling()); // floating point map object associated with the known mask
  for ( ix = mask_b.first(); !ix.last(); ix.next() )
  {
    mask_b_real[ix] = (float)mask_b[ix];
  }
  mask_b_real.fft_to (f_b); // FFT of the real-valued mask to get the Crystallographic structure factors
  
  bool allow_origin_shifts = true;
  
  //clipper::Coord_frac origin_shift( 0.0, 0.0, 0.0 );
  //origin_shift = clipper::Coord_frac(0.0, 0.0, 0.0);
  
  mask_cc=0.0;
  
  ipa_functions::phased_translation_function(f_a, f_b, allow_origin_shifts, origin_shift, mask_cc);
  if (verbose)
  {
    outlog << "Origin shift to maximize agreement between the masks: " << origin_shift.format() << "\n" << "Correlation coefficient from the phased translation function: " << mask_cc << std::endl;
  }
  // set return values to zero, for all real space agreement statistics
  
  simple_matching_coefficient=0.0;
  cohen_kappa = 0.0;
  mcc = 0.0;
  
  // If requested, apply the origin shift and figure out the actual real space agreement statistics
  
  
  if (compute_real_space_agreement)
  {
    
    // The old way
    /*
     int nt =0;
     int nm =0;
     int ns_a =0;
     int ns_b =0;
     
     for ( ix =  mask_a.first(); !ix.last(); ix.next() )
     {
     
     nt ++;
     if (mask_a[ix] == 0) {ns_a += 1;} // this voxel is solvent
     if (mask_b[ix] == 0) {ns_b += 1;} // this voxel is solvent
     
     clipper::Coord_frac cf = ix.coord().coord_frac(mask_a.grid_sampling()); // get the fractional coordinates of the current grid point
     cf -= origin_shift; // Apply the required origin shift. We want to transform the coordinates of mask B so they are expressed with respect to the same origin as mask A
     int mask_b_value = lroundf( mask_b_real.interp<clipper::Interp_linear>( cf ) );
     if ( mask_a[ix] == mask_b_value ) {nm ++;}
     }
     
     simple_matching_coefficient = float(nm)/float(nt);
     float fraction_solvent_a =    float(ns_a)/float(nt);
     float fraction_solvent_b =    float(ns_b)/float(nt);
     
     float expected_value_smc = fraction_solvent_a*fraction_solvent_b +(1.0-fraction_solvent_a)*(1.0-fraction_solvent_b);
     
     outlog << "Expected value of the Simple Matching coefficient (SMC) for two *random* masks: " << expected_value_smc << std::endl;
     
     
     outlog << "SMC between the masks: " << simple_matching_coefficient << std::endl;
     cohen_kappa = (simple_matching_coefficient - expected_value_smc) / (1.0 - expected_value_smc);
     outlog << "Cohen Coefficient (chance-corrected SMC) between the masks: " << cohen_kappa << std::endl;
     */
    
    // the new way
    /*
     
     int n11 = 0;
     int n00 = 0;
     int n10 = 0;
     int n01 = 0;
     
     int n0_a =0;
     int n0_b =0;
     
     for ( ix =  mask_a.first(); !ix.last(); ix.next() )
     {
     int mask_a_value = mask_a[ix];
     clipper::Coord_frac cf = ix.coord().coord_frac(mask_a.grid_sampling()); // get the fractional coordinates of the current grid point
     cf -= origin_shift; // Apply the required origin shift. We want to transform the coordinates of mask B so they are expressed with respect to the same origin as mask A
     int mask_b_value = lroundf( mask_b_real.interp<clipper::Interp_linear>( cf ) );
     
     if (mask_a_value == 0) {n0_a += 1;} // the point in mask a is solvent
     if (mask_b_value == 0) {n0_b += 1;} // correspondent point in mask b is solvent
     
     if ((mask_a_value == 1) && (mask_b_value == 1)) {n11 += 1;}
     if ((mask_a_value == 0) && (mask_b_value == 0)) {n00 += 1;}
     if ((mask_a_value == 1) && (mask_b_value == 0)) {n10 += 1;}
     if ((mask_a_value == 0) && (mask_b_value == 1)) {n01 += 1;}
     
     }
     
     // Compute bivariate realtive frequencies of agreement/disagreement from the raw counts
     
     int ntotal = n11 + n00 + n10 + n01;
     
     float f11 = (float)n11 / (float)ntotal ;
     float f00 = (float)n00 / (float)ntotal ;
     float f10 = (float)n10 / (float)ntotal ;
     float f01 = (float)n01 / (float)ntotal ;
     
     // Compute fraction solvent in each envelope from the raw counts
     
     float frac_solvent_a =    float(n0_a) / (float)ntotal;
     float frac_solvent_b =    float(n0_b) / (float)ntotal;
     */
    
    // the even newer way, with correct statistical weights for each point
    // I have not checked this exhaustively. In space groups without special positions it returns the unweighted result, as it should. Therefore probably fine.
    
    float n11 = 0;
    float n00 = 0;
    float n10 = 0;
    float n01 = 0;
    
    float n0_a =0;
    float n0_b =0;
    
    for ( ix =  mask_a.first(); !ix.last(); ix.next() )
    {
      // probably we should pass the weights array to this function rather than looking it up here every time, which is costly
      float statistical_weight = 1.0 / mask_a.multiplicity( ix.coord() );
      
      int mask_a_value = mask_a[ix];
      clipper::Coord_frac cf = ix.coord().coord_frac(mask_a.grid_sampling()); // get the fractional coordinates of the current grid point
      cf -= origin_shift; // Apply the required origin shift. We want to transform the coordinates of mask B so they are expressed with respect to the same origin as mask A
      int mask_b_value = roundf( mask_b_real.interp<clipper::Interp_linear>( cf ) );
      
      if (mask_a_value == 0) {n0_a += statistical_weight;} // the point in mask a is solvent
      if (mask_b_value == 0) {n0_b += statistical_weight;} // correspondent point in mask b is solvent
      
      if ((mask_a_value == 1) && (mask_b_value == 1)) {n11 += statistical_weight;}
      if ((mask_a_value == 0) && (mask_b_value == 0)) {n00 += statistical_weight;}
      if ((mask_a_value == 1) && (mask_b_value == 0)) {n10 += statistical_weight;}
      if ((mask_a_value == 0) && (mask_b_value == 1)) {n01 += statistical_weight;}
    }
    
    // Compute bivariate relative frequencies of agreement/disagreement from the raw counts
    
    float ntotal = n11 + n00 + n10 + n01;
    
    float f11 = n11 / ntotal ;
    float f00 = n00 / ntotal ;
    float f10 = n10 / ntotal ;
    float f01 = n01 / ntotal ;
    
    // Compute fraction solvent in each envelope from the raw counts
    
    float frac_solvent_a = n0_a / ntotal;
    float frac_solvent_b = n0_b / ntotal;
    
    // compute expectation value for the simple matching coefficient
    
    float e_smc = frac_solvent_a*frac_solvent_b +(1.0-frac_solvent_a)*(1.0-frac_solvent_b);
    if (verbose) {
    outlog << "Expected value of the Simple Matching coefficient (SMC) for two independent masks: " << e_smc << std::endl;
    
    // Compute real space agreement stats
    }
    simple_matching_coefficient =  f11 +  f00;
    cohen_kappa = 2*(f00*f11 - f01*f10) / ( (f00 + f01)*(f01 + f11) + (f00 + f10)*(f10 + f11) );
    mcc = (f00*f11 - f01*f10) / sqrt(( f00 + f01)*( f00 + f10)*( f10 + f11)*( f01 + f11));
    
    if (verbose)
    {
    outlog << "SMC between the masks: " << simple_matching_coefficient << std::endl;
    outlog << "Cohen kappa coefficient (chance-corrected SMC) between the masks: " << cohen_kappa << std::endl;
    outlog << "Matthews correlation coefficient between the masks: " << mcc << std::endl;
    }
    
  }
}

clipper::Xmap<int> ipa_functions::complement_mask(const clipper::Xmap<int> &input_mask)
{
  
  // Returns the complement of a binary-valued mask.
  
  clipper::Xmap<int>    complement_mask(input_mask.spacegroup(), input_mask.cell(), input_mask.grid_sampling());
  clipper::Xmap_base::Map_reference_index ix;
  
  
  for ( ix = input_mask.first(); !ix.last(); ix.next() )
  {
    if (input_mask[ix] == 0)
    {
      complement_mask[ix] = 1;
    }
    else
    {
      complement_mask[ix] = 0;
    }
  }
  
  return complement_mask;
  
}

clipper::Xmap<int> ipa_functions::shift_mask(const clipper::Xmap<int> &input_mask, const clipper::Coord_frac &origin_shift)
{
  
  // Applies an origin shift to a binary-valued mask.
  // No checks on the validity of the origin shift are perfomed - that's up to the user. GIGO
  
  //clipper::Xmap<float>  working_mask(input_mask.spacegroup(), input_mask.cell(), input_mask.grid_sampling());
  clipper::Xmap<int>    shifted_mask(input_mask.spacegroup(), input_mask.cell(), input_mask.grid_sampling());
  clipper::Xmap_base::Map_reference_index ix;
  clipper::Xmap_base::Map_reference_coord is(input_mask);
  
  /*
   for ( ix = input_mask.first(); !ix.last(); ix.next() )
   {
   working_mask[ix] = (float)input_mask[ix]; // making a floating point version of the integer-valued mask
   }
   */
  
  for ( ix = input_mask.first(); !ix.last(); ix.next() )
  {
    clipper::Coord_frac cf = ix.coord().coord_frac(input_mask.grid_sampling()); // get the fractional coordinates of the current grid point
    clipper::Coord_frac tf = cf + origin_shift ;
    clipper::Coord_grid tg = tf.coord_grid(input_mask.grid_sampling()); // now convert the shifted coordinates back to grid coordinates
    
    shifted_mask[ix] =  input_mask[is.set_coord(tg)];
    
    //cf -= origin_shift; // Apply the required origin shift.
    //shifted_mask[ix] = roundf(working_mask.interp<clipper::Interp_linear>( cf ) );
    
  }
  
  return shifted_mask;
  
}

void ipa_functions::shift_phases(const clipper::HKL_data<clipper::data64::F_phi> &input_fp, clipper::HKL_data<clipper::data64::F_phi> &output_fp, const clipper::Coord_frac &origin_shift)
{
  // Applies the specified origin shift input_fp, passes back the result in output_fp
  // No checks on the validity of the origin shift are perfomed - that's up to the user. GIGO
  
  /*
   const clipper::Spacegroup spacegroup   = input_fp.base_hkl_info().spacegroup();
   const clipper::Cell       cell         = input_fp.base_hkl_info().cell();
   const clipper::Resolution resolution   = input_fp.base_hkl_info().resolution();
   
   clipper::HKL_info hkl(spacegroup, cell, resolution);
   hkl.generate_hkl_list(); // initially there are no observations in the list - better make some
   
   clipper::HKL_data<clipper::data64::F_phi> fp_shifted(hkl);
   */
  
  clipper::HKL_data_base::HKL_reference_index ih;
  
  for ( ih = input_fp.first(); !ih.last(); ih.next() )
  {
    // Convert the supplied origin shift into a phase difference
    
    double delta_phi;
    clipper::Coord_reci_frac h( ih.hkl() ); // fractional coordinates of current observation ih
    const double hx = h * origin_shift;
    delta_phi = clipper::Util::twopi() * ( hx - floor(hx) );
    
    //fp_shifted[ih].f()   = input_fp[ih].f(); // Leave the amplitude unchanged
    //fp_shifted[ih].phi() = input_fp[ih].phi() + delta_phi;  // Apply the phase shift
    
    output_fp[ih].f()   = input_fp[ih].f(); // Leave the amplitude unchanged
    output_fp[ih].phi() = input_fp[ih].phi() + delta_phi;  // Apply the phase shift
  }
  
  //return fp_shifted;
}

float ipa_functions::correlation_based_distance_measure(const float &cc, const int &distance_measure)
{
  float distance;
  
  switch (distance_measure)
  {
    case 1 :
      distance = 1.0 - abs(cc);
      break;
    case 2 :
      distance = sqrt( 1.0 - std::pow(cc,2) );
      break;
  }
  
  return distance;
  
}

clipper::Xmap<int> ipa_functions::invert_mask(const clipper::Xmap<int> &input_mask)
{
  
  // Returns a binary-valued crystallographic mask that has been inverted through the appropriate point in the unit cell (generally, but not always, the origin)
  // This procedure should only be performed in one of the 43 achiral space groups relevant to protein crystallography
  // When applied to one of the 22 chiral space groups it will produce nonsense, as in these cases, inversion of the image changes the space group (i.e you can't invert in a chiral space group)
  
  // For space groups I41, I4122 and I4132 the required inversion is not through the origin. Yet more trouble created by centering
  // For discussion see McCoy AJ, Read RJ. Experimental phasing: best practice and pitfalls. Acta Crystallogr D. 2010 Apr;66(Pt 4):458–69.
  
  clipper::Xmap<int> invert_mask(input_mask.spacegroup(), input_mask.cell(), input_mask.grid_sampling());
  clipper::Xmap_base::Map_reference_index ix;
  clipper::Xmap_base::Map_reference_coord is(input_mask);
  
  // TODO ... insert a safety check- if space group is chiral return the unmodified mask
  
  // Figure out what point we are inverting through, which is spacegroup dependent
  
  const int spacegroup_number = input_mask.spacegroup().spacegroup_number();
  
  clipper::Coord_frac inversion_point;
  
  if      (spacegroup_number == 80) // space group I41
  {
    // if we invert through the origin (x -> -x) then we have to shift the origin of the resultant function to (1/2, 0, 0)
    // So original point x,y,z -> -x,-y,-z (inversion through origin) -> -x-1/2, -y, -z (origin shift)
    // Alternately we could perform the inversion directly through the point p, for which the formula is (x -> 2p - x)
    // The required point p is (-1/4, 0, 0)
    inversion_point = clipper::Coord_frac(-0.250,  0.000,  0.000);
    
  }
  else if (spacegroup_number == 98) // space group I4122
  {
    // if we invert through the origin then we have to shift the origin of the resultant function to (1/2, 0, 1/4)
    // So original point x,y,z -> -x,-y,-z (inversion through origin) -> -x-1/2, -y, -z-1/4 (origin shift)
    // Alternately we could perform the inversion directly through the point p, for which the formula is (x -> 2p - x)
    // The required point p is (-1/4, 0, -1/8)
    inversion_point = clipper::Coord_frac(-0.250,  0.000, -0.125);
  }
  else if (spacegroup_number == 214) // space group I4132
  {
    // if we invert through the origin then we have to shift the origin of the resultant function to  (1/4, 1/4, 1/4).
    // So original point x,y,z -> -x,-y,-z (inversion through origin) -> -x-1/4, -y-1/4, -z-1/4 (origin shift)
    // Alternately we could perform the inversion directly through the point p, for which the formula is (x -> 2p - x)
    // The required point p is (-1/8, -1/8, -1/8)
    inversion_point = clipper::Coord_frac(-0.125, -0.125, -0.125);
  }
  else
  {
    // inversion through the origin - nice and simple
    inversion_point = clipper::Coord_frac( 0.000,  0.000,  0.000);
    
  }
  
  for ( ix = invert_mask.first(); !ix.last(); ix.next() )
  {
    clipper::Coord_frac cf = ix.coord().coord_frac(input_mask.grid_sampling()); // get the current fractional coordinates
    clipper::Coord_frac tf = 2*inversion_point - cf; // and here are our fractional coordinates following inversion through the relevant point
    clipper::Coord_grid tg = tf.coord_grid(input_mask.grid_sampling()); // now convert our fractional coordinates back to grid coordinates
    
    invert_mask[ix] =  input_mask[is.set_coord(tg)];
    
    // And here's the simpler version - pure inversion through the origin
    //clipper::Coord_grid cg = ix.coord(); // get the  the current grid point
    //invert_mask[ix] =  input_mask[is.set_coord(-cg)];
  }
  
  return invert_mask;
  
  
}

void ipa_functions::perform_ncs_density_averaging(clipper::Xmap<float> &map, const clipper::Xmap<int> &ncs_avg_mask, const std::vector< std::vector< clipper::RTop_frac > > &ncs_operations)
{
  clipper::Xmap_base::Map_reference_index ix;
  clipper::Xmap<float>  averaged_map( map.spacegroup(), map.cell(), map.grid_sampling() );  // create map object to temporarily store the averaged electron density
  
  for ( ix = map.first(); !ix.last(); ix.next() )
  {
    if (ncs_avg_mask[ix] == 0) // we are outside the mask, do nothing
    {
      averaged_map[ix] = map[ix];
    }
    else                    // we are inside the mask, average the density using the appropriate symmetry operators
    {
      averaged_map[ix] = map[ix];
      clipper::Coord_frac cf = ix.coord().coord_frac(map.grid_sampling()); // get the fractional coordinates of the current grid point
      
      int assembly_id = ncs_avg_mask[ix]-1; // recall that the mask subregion for assembly i is denoted with integer i+1.
      int pg_order = ncs_operations[0].size(); // Vector of Vectors is rectangular so this is what we need.
      
      //      outlog << pg_order << " " << assembly_id << std::endl;
      
      for (int ncs_op = 1; ncs_op<=pg_order-1 ; ncs_op++ )
      {
        clipper::Coord_frac cft = cf.transform(ncs_operations[assembly_id][ncs_op]); // apply the relevant symmetry operation to get the transformed fractional coordinates
        clipper::Coord_map  cmt = cft.coord_map(map.grid_sampling()); // convert the transformed fractional coordinate into a map coordinate
        float interpolated_density;
        clipper::Interp_linear::interp(map, cmt, interpolated_density);
        averaged_map[ix] += interpolated_density;
      }
      averaged_map[ix] /= float(pg_order);
    }
  }
  for ( ix = map.first(); !ix.last(); ix.next() )
  {
    map[ix] = averaged_map[ix];
  }
  
}

void ipa_functions::calculate_summary_statistics(const clipper::Xmap<float> &map, const clipper::Xmap<int> &mask, const clipper::Xmap<float> &weights, const int &region_flag, const bool &calculate_higher_moments, 
                                                 double &mean, double &variance, double &skewness, double &kurtosis, double &minimum, double &maximum)
{
  // Evaluates summary statistics for some region of a crystallographic map.
  // The region is defined by an accompanying integer-valued mask
  // IPA convention:
  // Solvent region is defined by mask values = 0
  // Protein region(s) are defined by mask values  > 0
  //
  // If the region_flag >= 0  -> all points in the map with correspondent mask value = region_flag will be used for calculating statistics
  // If the region_flag  = -1 -> all points in the map will be used for calculating statistics
  //
  // In the calculation, map points are weighted according to their multiplicity
  // These weights compensate for the fact that we retain only the asymmetric unit of the electron density map in real space and we want any statistics we calculate to relate to a complete unit cell.
  // Points lying on crystallographic "special positions" (non-translational symnmetry elements of the crystal) need to be downweighted so that they end up being "counted" the correct number of times.
  //
  // The function always returns the mean, minimum and maximum
  // If requested  (calculate_higher_moments - true) it returns the variance, skewness and kurtosis.
  //
  //  Employs a two-pass algorithm which is numerically stable but slow (see e.g. Numerical Recipes in Fortran 77)
  //  TODO ... Implement a stable single pass algorithm to compute the variance
  //
  
  clipper::Xmap_base::Map_reference_index ix;
  
  
  // First pass through the map - compute the mean, minimum, maximum
  
  double sumx = 0.0;
  double sumw = 0.0;
  double n;
  
  ix = map.first();
  minimum = map[ix];
  maximum = map[ix];
  
  for ( ix = map.first(); !ix.last(); ix.next() )
  {
    //float weight = 1.0 / map.multiplicity(ix.coord());
    if ( (mask[ix] == region_flag) || (region_flag == -1) )
    {
      if (map[ix] < minimum)
      {minimum = map[ix];}
      if (map[ix] > maximum)
      {maximum = map[ix];}
      sumx +=  weights[ix]*map[ix];
      sumw +=  weights[ix];
    }
  }
  
  n = sumw;
  mean = sumx/n;
  
  // If requested, make a second pass through the map, to compute the second, third and fourth sample central moments about the mean
  
  if (calculate_higher_moments)
  {
    //double m1 = mean;
    double m2 = 0.0;
    double m3 = 0.0;
    double m4 = 0.0;
    
    for ( ix = map.first(); !ix.last(); ix.next() )
    {
      //float weight = 1.0 / map.multiplicity(ix.coord());
      if ( (mask[ix] == region_flag) || (region_flag == -1) )
      {
        double s = weights[ix]*map[ix] - mean;
        double p = s*s;
        m2 = m2 + p;
        p = p*s;
        m3 = m3 + p;
        p = p*s;
        m4 = m4 + p;
      }
    }
    
    m2 = m2/n;
    m3 = m3/n;
    m4 = m4/n;
    
    // Compute the k-statistics - the unique symmetric unbiased estimators of the cumulants
    // see e.g. http://mathworld.wolfram.com/k-Statistic.html
    
    // double k1 =  m1;
    double k2 = (n*m2) / (n-1);
    double k3 = ((n*n)*m3) / ((n-1)*(n-2));
    double k4 =  (n*n)*((n+1)*m4-3.0*(n-1)*(m2*m2)) / ((n-1)*(n-2)*(n-3));
    
    //  Finally - output the generally accepted estimators of the population variance, skewness and kurtosis
    
    // mean = k1;
    variance = k2;
    skewness = k3/(std::pow(k2,3/2));
    kurtosis = k4/(std::pow(k2,2  ));
  }
  
}

void ipa_functions::calculate_histogram(const clipper::Xmap<float> &map, const clipper::Xmap<int> &mask, const clipper::Xmap<float> &weights, const int &region_flag, histogram &hg, std::ostream &outlog, bool verbose)
{
  
  // Calculates a frequency histogram for some region of a crystallographic map.
  
  // compute the histogram scale factor (the inverse of the bin width)
  
  double n_hist_bins = (double)hg.number_of_intervals;
  double histogram_scale_factor = (n_hist_bins-1.0)/(hg.maximum - hg.minimum);
  
  // For safety, ensure that the input histogram is zero valued when we start
  
  for ( int i= 0; i < hg.number_of_intervals; i++ )
  {
    hg.interval[i] = 0.0;
  }
  
  // Make a single pass through the map and compute the histogram
  
  
  clipper::Xmap_base::Map_reference_index ix;
  
  float sumw=0.0;
  
  for ( ix = map.first(); !ix.last(); ix.next() )
  {
    //float weight = 1.0 / map.multiplicity(ix.coord());
    if ( (mask[ix] == region_flag) || (region_flag == -1) )
    {
      sumw += weights[ix];
      
      int index = std::trunc(histogram_scale_factor * (map[ix] - hg.minimum)); // trunc rounds the argument towards zero and returns the nearest integral value that is not larger in magnitude than the argument.
      
      // safe navigation if the maximum and minimum values for the histogram sit outside the values found in the map.
      if (index < 0) {index = 0;}
      if (index > hg.number_of_intervals-1) {index = hg.number_of_intervals-1;}
      
      hg.interval[index] += weights[ix]; // this accounts correctly for the multiplicity of each point
      
    }
  }
  
  // Normalize and report the histogram
  if (verbose)
  {
    outlog << std::endl;
    outlog << "Normalized Histogram" << std::endl;
    outlog << std::fixed;
  }
  
  for ( int i= 0; i < hg.number_of_intervals; i++ )
  {
    hg.interval[i] = hg.interval[i]/sumw;
    
    if (verbose)
    {
      outlog << std::setprecision(0) << i << " " << std::setprecision(6) << hg.minimum + (float)i / histogram_scale_factor << " " << hg.interval[i] << "\n";
    }
  }

  if (verbose)
  {
    outlog << std::endl;
  }
}

void ipa_functions::flatten_map_region(clipper::Xmap<float> &map, const clipper::Xmap<int> &mask, const int &region_flag, const float &value)
{
  // Set the requested region of a map to a constant value
  clipper::Xmap_base::Map_reference_index ix;
  for ( ix = map.first(); !ix.last(); ix.next() )
  {
    if ( (mask[ix] == region_flag) || (region_flag == -1) )
    {
      map[ix] = value;
    }
  }
}

void ipa_functions::shift_map_region(clipper::Xmap<float> &map, const clipper::Xmap<int> &mask, const int &region_flag, const float &value)
{
  // Add a constant value to the requested region of a map
  
  clipper::Xmap_base::Map_reference_index ix;
  for ( ix = map.first(); !ix.last(); ix.next() )
  {
    if ( (mask[ix] == region_flag) || (region_flag == -1) )
    {
      map[ix] += value;
    }
  }
  
}

void ipa_functions::histogram_match_map_region(clipper::Xmap<float> &map, const clipper::Xmap<int> &mask, const int &region_flag, const histogram &observed_histogram, const histogram &reference_histogram, const double &a, const double &b, const bool &integrate_from_top)
{
  
  //  Compute the modifying function (lookup table).
  //  For each bin in the observed histogram, the modifying function tells us which bin it maps to in the reference histogram.
  
  double ccdf_ref = 0.0; // This is the complementary cumulative distribution function of the reference histogram
  double ccdf_obs = 0.0; // This is the complementary cumulative distribution function of the observed histogram
  
  std::vector<int> modifying_function(observed_histogram.number_of_intervals, 0); // This is the modifying function, initialize with all values zero
  
  if (integrate_from_top)
  {
    // Integrate each histogram from the top
    
    modifying_function[observed_histogram.number_of_intervals-1] = observed_histogram.number_of_intervals-1;
    int integration_bound = observed_histogram.number_of_intervals-1;
    
    for ( int ir = reference_histogram.number_of_intervals-1; ir >= 0; ir--)
    {
      
      ccdf_ref += (double)reference_histogram.interval[ir];
      //outlog << "ccdf_ref" << ccdf_ref << std::endl;
      
      while ( (ccdf_ref > ccdf_obs) && (integration_bound >= 0) )
      {
        ccdf_obs = ccdf_obs + (double)observed_histogram.interval[integration_bound];
        //outlog << "ccdf_obs" << ccdf_obs << std::endl;
        modifying_function[integration_bound] = ir;
        integration_bound = integration_bound - 1;
      }
      
    }
  }
  else
  {
    // Integrate each histogram from the bottom
    
    int integration_bound = 0;
    
    for ( int ir = 0 ; ir < reference_histogram.number_of_intervals; ir++)
    {
      
      ccdf_ref += (double)reference_histogram.interval[ir];
      //outlog << "ccdf_ref" << ccdf_ref << std::endl;
      
      while ( (ccdf_ref > ccdf_obs) && (integration_bound < observed_histogram.number_of_intervals) )
      {
        ccdf_obs = ccdf_obs + (double)observed_histogram.interval[integration_bound];
        //outlog << "ccdf_obs" << ccdf_obs << std::endl;
        modifying_function[integration_bound] = ir;
        integration_bound = integration_bound + 1;
      }
    }
    
  }
  
  //  for ( int ir = reference_histogram.number_of_intervals-1; ir >= 0; ir--)
  //    {
  //        outlog << ir << " " << modifying_function[ir] << std::endl;
  //    }
  
  // Now perform the remapping
  
  // First apply linear rescaling to the reference histogram
  
  double reference_minimum = a*reference_histogram.minimum + b;
  double reference_maximum = a*reference_histogram.maximum + b;
  
  // Now compute the histogram scale factors
  
  double histogram_scale_factor_obs = ((double)observed_histogram.number_of_intervals-1.0)/(observed_histogram.maximum - observed_histogram.minimum);
  double histogram_scale_factor_ref = ((double)reference_histogram.number_of_intervals-1.0)/(reference_maximum  - reference_minimum);
  
  //  outlog << histogram_scale_factor_obs  <<  std::endl;
  //  outlog << histogram_scale_factor_ref  <<  std::endl;
  
  // Report the rescaled reference histogram
  /*
   outlog << std::endl;
   outlog << "Reference Histogram after linear rescaling has been applied " << std::endl;
   outlog << std::endl;
   
   for ( int i= 0; i < reference_histogram.number_of_intervals; i++ )
   {
   outlog << i << " " <<  reference_minimum + ((float)i / histogram_scale_factor_ref) << " " << reference_histogram.interval[i] << std::endl;
   }
   */
  
  
  clipper::Xmap_base::Map_reference_index ix;
  
  for ( ix = map.first(); !ix.last(); ix.next() )
  {
    if ( (mask[ix] == region_flag) || (region_flag == -1) )
    {
      int index_obs = std::trunc( histogram_scale_factor_obs * (map[ix] - observed_histogram.minimum) ); // This is the bin index in the observed (original) histogram
      
      //    if ((index_obs < 0) || (index_obs > observed_histogram.number_of_intervals-1))
      //    {
      //     outlog << "Danger Will Robinson " <<  index_obs <<  std::endl;
      //     outlog << histogram_scale_factor_obs <<  std::endl;
      //     outlog << map[ix] << std::endl;
      //     outlog << observed_histogram.minimum << std::endl;
      //    }
      
      map[ix] =  reference_minimum + ( (double)modifying_function[index_obs] / histogram_scale_factor_ref); // This remaps the density according to the modifying function
                                                                                                            //    outlog << index_obs << " " << map[ix] << std::endl;
    }
  }
  
}

double ipa_functions::earth_movers_distance(const histogram &hg1, const histogram &hg2)
{
  
  // Computes the earth mover's (Wasserstein) distance between two histograms with the same bin intervals
  // See https://math.stackexchange.com/questions/714476/how-do-you-compute-numerically-the-earth-movers-distance-emd
  // No checks are performed that the input histograms are sensibly binnned, that's up to the user.
  
  std::vector<double> cdf1(hg1.number_of_intervals, 0.0); // cumulative density function for histogram 1
  std::vector<double> cdf2(hg1.number_of_intervals, 0.0); // cumulative density function for histogram 2
  
  
  /// first compute the cdfs
  
  cdf1[0] = hg1.interval[0];
  cdf2[0] = hg2.interval[0];
  
  for ( int i= 1; i < hg1.number_of_intervals; i++)
  {
    cdf1[i] = cdf1[i-1] + hg1.interval[i];
    cdf2[i] = cdf2[i-1] + hg2.interval[i];
  }
  
  // compute the histogram bin width (must be the same for both input histograms or we are in deep trouble)
  
  // double histogram_bin_width1 = (hg1.maximum - hg1.minimum)/( (double)(hg1.number_of_intervals-1) );
  
  /// now compute the earth movers distance
  
  double earth_movers_distance = 0.0;
  
  for ( int i= 0; i < hg1.number_of_intervals; i++)
  {
    cdf1[i] = cdf1[i]/cdf1[hg1.number_of_intervals-1]; // normalize cdf1
    cdf2[i] = cdf2[i]/cdf2[hg1.number_of_intervals-1]; // normalize cdf2
    earth_movers_distance += std::abs(cdf1[i] - cdf2[i]);
  }
  
  // earth_movers_distance *= histogram_bin_width1;
  
  return earth_movers_distance;
}

void ipa_functions::reconstruct_image_from_gradient(const clipper::HKL_info &hkl_info, const clipper::Xmap<float>  &gradient_x, const clipper::Xmap<float>  gradient_y, const clipper::Xmap<float> &gradient_z, clipper::Xmap<float> output_map)
{
  // TBC    
}

// TODO  break this function into parts, which perform the individual projections onto the constraints (solvent flattening, histogram matching etc) 
void ipa_functions::density_modification (const bool &solvent_flattening,
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
                                          bool verbose)
{
  
  
  
  
  bool calculate_higher_moments;
  
  // These are currently also defined as member variables of the ipa-worker object, TODO: Remove hard coding.
  const int solvent_region_flag = 0;
  const int protein_region_flag = 1;
  
  map_region_summary_statistics solvent_density_stats;
  map_region_summary_statistics protein_density_stats;
  
  map_region_summary_statistics protein_gradient_stats;
  
  clipper::Xmap_base::Map_reference_index ix;
  
  // First make an unmodified copy of the input map, to which we'll apply the density modification operations
  for ( ix =  output_map.first(); !ix.last(); ix.next() )
  {
    output_map[ix] = input_map[ix];
  }
  
  //  Solvent flattening
  
  if (solvent_flattening)
  {
    
    outlog << "Solvent Flattening " << std::endl;
    
    // First compute some statistics for the solvent region of the input map
    
    calculate_higher_moments = false;
    ipa_functions::calculate_summary_statistics(input_map, mask, weights, solvent_region_flag, calculate_higher_moments, solvent_density_stats.mean, solvent_density_stats.variance, solvent_density_stats.skewness, solvent_density_stats.kurtosis, solvent_density_stats.minimum, solvent_density_stats.maximum);
    
    outlog << "Mean density value in the solvent region: " << solvent_density_stats.mean << std::endl;
    
    // Impose some reasonable physical constraints on the mean solvent density
    
    if (solvent_density_stats.mean > 1.1*expected_mean_density_solvent)
    {
      solvent_density_stats.mean =  1.1*expected_mean_density_solvent;
      outlog << "Mean in the solvent region higher than expected - capping the density at 1.1 * the expected value" << std::endl;
    }
    
    if (solvent_density_stats.mean < 0.9*expected_mean_density_solvent)
    {
      solvent_density_stats.mean =  0.9*expected_mean_density_solvent;
      outlog << "Mean in the solvent region lower than expected - capping the density at 0.9 * the expected value "  <<  std::endl;
    }
    
    ipa_functions::flatten_map_region(output_map, mask, solvent_region_flag, solvent_density_stats.mean);
    
    
    outlog << "Assigned density value in the solvent region: " << solvent_density_stats.mean << std::endl;
  }
  
  if      (histogram_matching)
  {
    
    outlog << "Density Histogram Matching " << std::endl;
    
    // We need to calculate the electron density histogram  for the input map
    
    // First compute some statistics for the protein region of the input  map
    
    calculate_higher_moments = true;
    ipa_functions::calculate_summary_statistics(input_map, mask, weights, protein_region_flag, calculate_higher_moments, protein_density_stats.mean, protein_density_stats.variance, protein_density_stats.skewness, protein_density_stats.kurtosis, protein_density_stats.minimum, protein_density_stats.maximum);
    
    outlog <<  "Mean in the Protein Region:     " << protein_density_stats.mean << std::endl;
    outlog <<  "Variance in the Protein Region: " << protein_density_stats.variance << std::endl;
    outlog <<  "Minimum in the Protein Region:  " << protein_density_stats.minimum << std::endl;
    outlog <<  "Maximum in the Protein Region:  " << protein_density_stats.maximum << std::endl;
    
    density_histogram.minimum  = protein_density_stats.minimum;
    density_histogram.maximum  = protein_density_stats.maximum;
    
    // Now compute the histogram
    
    if (verbose) { outlog << "Electron density histogram for protein region of the map, prior to modification" << std::endl; }
    ipa_functions::calculate_histogram(input_map, mask, weights, protein_region_flag, density_histogram, outlog, verbose);
    
    // compute the scaling parameters for the reference ed histogram
    
    double a;
    double b;
    
    // Suppose we have reference density with mean E(X) and Variance Var(X)
    // We need a linear transformation of the reference density ( Y=aX+b ), such that we end up with the required mean E(Y) and variance Var(Y).
    // Under this linear transformation E(Y) = aE(X)+ b and Var(Y) = a^2 * Var(X)
    // Hence the parameters of the linear transfornation are  a = Sqrt( Var(Y)/Var(X) ), and b = E(Y)- aE(X)
    
    // In  this instance
    // E(X) = mean of the reference density in the protein region
    // Var(X) = variance of the reference density in the protein region
    // E(Y) = expected mean density in the protein region
    // Var(Y) = variance of the current map in the protein region
    
    
    a = sqrt( protein_density_stats.variance / protein_density_stats_reference.variance); //  a = Sqrt( Var(Y)/Var(X) )
    b = expected_mean_density_protein - (a * protein_density_stats_reference.mean); // b = E(Y)- aE(X)
    
    outlog << std::endl;
    outlog << "Linear rescaling factors for the reference electron density, so it has the expected mean protein density, and the variance of the input map " << std::endl;
    outlog <<" slope (a):      " << a << std::endl;
    outlog <<" intercept (b):  " << b << std::endl;
    
    outlog << "Following rescaling" << std::endl;
    outlog << "Reference Histogram Minimum Density: " << a*density_histogram_reference.minimum + b << std::endl;
    outlog << "Reference Histogram Maximum Density: " << a*density_histogram_reference.maximum + b << std::endl;
    
    outlog << std::endl;
    outlog << "Remapping the electron density histogram" << std::endl;
    
    bool integrate_from_top = true;
    ipa_functions::histogram_match_map_region(output_map, mask, protein_region_flag, density_histogram, density_histogram_reference, a, b, integrate_from_top); // perform the histogram matching on the protein region
    
  }
  else if (gradient_matching)
  {
    
    outlog << "Gradient Matching " << std::endl;
    
    
    
    clipper::Xmap<float>  gradient_x   ( input_map.spacegroup(), input_map.cell(), input_map.grid_sampling() );  // create map object to store gradient component along X
    clipper::Xmap<float>  gradient_y   ( input_map.spacegroup(), input_map.cell(), input_map.grid_sampling() );  // create map object to store gradient component along Y
    clipper::Xmap<float>  gradient_z   ( input_map.spacegroup(), input_map.cell(), input_map.grid_sampling() );  // create map object to store gradient component along Z
    clipper::Xmap<float>  gradient_magnitude( input_map.spacegroup(), input_map.cell(), input_map.grid_sampling() );  // create map object to store the gradient magnitude
    clipper::Xmap<float>  gradient_magnitude_remapped( input_map.spacegroup(), input_map.cell(), input_map.grid_sampling() );  // create map object to store the remapped gradient magnitude
    
    
    // First we need to calculate the gradient associated with the input map
    
    clipper::HKL_data<clipper::data64::F_phi> work_fp (hkl_target);
    
    input_map.fft_to( work_fp ); // FFT of the input map to get the Fourier Coefficients
    ipa_functions::calculate_gradient(work_fp, gradient_x, gradient_y, gradient_z, gradient_magnitude, outlog);
    
    // Now we need to calculate the gradient magnitude histogram
    
    // First compute some statistics for the gradient magnitude map
    
    bool calculate_higher_moments = true;
    ipa_functions::calculate_summary_statistics(gradient_magnitude, mask, weights, protein_region_flag, calculate_higher_moments,
                                                protein_gradient_stats.mean, protein_gradient_stats.variance, protein_gradient_stats.skewness, protein_gradient_stats.kurtosis, protein_gradient_stats.minimum, protein_gradient_stats.maximum);
    
    gradient_histogram.minimum = protein_gradient_stats.minimum;
    gradient_histogram.maximum = protein_gradient_stats.maximum;
    
    // outlog << gradient_histogram.minimum << std::endl;
    // outlog << gradient_histogram.maximum << std::endl;
    
    // Now compute the histogram
    
    if (verbose) { outlog << "Gradient magnitude histogram for protein region of the map, prior to modification" << std::endl; }
    ipa_functions::calculate_histogram(gradient_magnitude, mask, weights, protein_region_flag, gradient_histogram, outlog);
    
    // compute the scaling parameter for the reference gradient magnitude histogram
    
    double a;
    double b;
    
    // Suppose the reference gradient magnitude map has Variance Var(X) in the protein region.
    // We need a linear transformation of the reference gradient magnitudes ( Y=aX), such that we end up with the variance Var(Y) of our current gradient magnitude map
    // Under this linear transformation E(Y) = aE(X) and Var(Y) = a^2 * Var(X)
    // Hence the parameter of the linear transfornation a = Sqrt( Var(Y)/Var(X) )
    
    // In this instance
    // Var(X) = variance of the reference gradient magnitude map in the protein region
    // Var(Y) = variance of the current gradient magnitude map in the protein region
    
    a = sqrt( protein_gradient_stats.variance / protein_gradient_stats_reference.variance); //  a = Sqrt( Var(Y)/Var(X) )
    b = 0;
    
    outlog << std::endl;
    outlog << "Linear rescaling factor for the reference gradient magnitude, so it has the variance of the input gradient magnitude map " << std::endl;
    outlog <<" slope (a):      " << a << std::endl;
    
    
    outlog << "Following rescaling" << std::endl;
    outlog << "Reference gradient magnitude histogram minimum: " << a*gradient_histogram_reference.minimum << std::endl;
    outlog << "Reference gradient magnitude histogram maximum: " << a*gradient_histogram_reference.maximum << std::endl;
    
    outlog << std::endl;
    outlog << "Remapping the gradient magnitude histogram" << std::endl;
    
    bool integrate_from_top = false;
    
    // make a copy of the gradient magnitude map that we can modify
    
    for ( ix = gradient_magnitude_remapped.first(); !ix.last(); ix.next() )
    {
      gradient_magnitude_remapped[ix] =  gradient_magnitude[ix];
    }
    
    ipa_functions::histogram_match_map_region(gradient_magnitude_remapped, mask, protein_region_flag, gradient_histogram, gradient_histogram_reference, a, b, integrate_from_top); // remap the gradient magnitude in the protein region
    
    // rescale the gradient components according to the gradient magnitude remapping
    
    outlog << "Rescaling the individual gradient magnitude components Gx, Gy and Gz so the gradient magnitude histogram matches the reference" << std::endl;
    
    for ( ix = gradient_magnitude_remapped.first(); !ix.last(); ix.next() )
    {
      float scale_factor = gradient_magnitude_remapped[ix]/gradient_magnitude[ix];
      if (std::isinf(scale_factor)) {scale_factor = 0.0;}  // Catch division by zero. If input gradient magnitude was zero, it will be remapped to zero, so scale factor should be zero.
      gradient_x[ix] *= scale_factor;
      gradient_y[ix] *= scale_factor;
      gradient_z[ix] *= scale_factor;
    }
    
    outlog << "Reconstructing the image having the required Gradient magnitudes" << std::endl;
    ipa_functions::reconstruct_image_from_gradient(work_fp.base_hkl_info(),gradient_x, gradient_y, gradient_z, output_map);
    
    // Examine the output
    // This can be omitted once we are sure everything is working as it should
    
    output_map.fft_to( work_fp ); // FFT of the output map to get the Fourier Coefficients
    
    ipa_functions::calculate_gradient(work_fp, gradient_x, gradient_y, gradient_z, gradient_magnitude, outlog);
    
    ipa_functions::calculate_summary_statistics(gradient_magnitude, mask, weights, protein_region_flag, calculate_higher_moments,
                                                protein_gradient_stats.mean, protein_gradient_stats.variance, protein_gradient_stats.skewness, protein_gradient_stats.kurtosis, protein_gradient_stats.minimum, protein_gradient_stats.maximum);
    
    gradient_histogram.minimum = protein_gradient_stats.minimum;
    gradient_histogram.maximum = protein_gradient_stats.maximum;
    
    
    outlog << "Gradient magnitude histogram for protein region of the map, after modification" << std::endl;
    ipa_functions::calculate_histogram(gradient_magnitude, mask, weights, protein_region_flag, gradient_histogram, outlog);
    
    
  }
  else // If we are not histogram or gradient matching, we simply shift the mean value in the protein region to its expected value
  {
    outlog << "\n" << std::endl;
    outlog << "Shifting density in the protein region to achieve the expected mean value" << std::endl;
    outlog << "\n" << std::endl;
    
    double shift = expected_mean_density_protein - protein_density_stats.mean;
    
    ipa_functions::shift_map_region(output_map, mask, solvent_region_flag, shift);
    
  }
  
}    

double ipa_functions::rmsd_between_maps(const clipper::Xmap<float> &map1, const clipper::Xmap<float> &map2, const clipper::Xmap<float> &weights)
{
  double sumw = 0.0;
  double sumsq = 0.0;
  
  clipper::Xmap_base::Map_reference_index ix;
  
  for ( ix = map1.first(); !ix.last(); ix.next() )
  {
    sumsq += std::pow((map1[ix] - map2[ix]),2)*weights[ix];
    sumw  += weights[ix];
  }
  
  double rmsd =  sqrt(sumsq/sumw);
  
  return rmsd;
  
}

double ipa_functions::flattened_sinusoid(const double b, const double phase, const double integral)
{
  /*
   The basic function here is a "flattened" sinusoidal wave √(1+b^2)/(1+b^2*cos^2(2πxf)) * cos(2πxf)
   see https://math.stackexchange.com/questions/100655/cosine-esque-function-with-flat-peaks-and-valleys
   with a phase offset added in
   
   Here f is the frequency (1/iterations) and x is the advancement (iterations)
   
   The unique parameter of the function is the dimensionless (b) which controls the "squareness" of the sinusoid
   
   The function returns a value between 1 amd -1
   
   An additional complexity is that the function may be "chirped" - with a linear change in either the period or the frequency as the function advances.
   
   To achieve this we need to integrate the function describing the change in frequency, and use that as an argument to the cosinusoidal function.
   
   The necessary integral must be computed and supplied by the calling function
   
   
   */
  
  
  double b_squared = b * b;
  
  double output = std::cos(phase + clipper::Util::twopi()*integral);
  
  output = sqrt ( (1.0 + b_squared) / ( 1.0 + ( b_squared * output * output )) ) * output;     // Apply the modifier which imposes "squareness" on the function
  
  return output;
  
}

double ipa_functions::pulse_wave(const double pulse_duration, const double phase, const double integral)

{
  /*
   The basic function here is a rectangular pulse wave
   
   H[ cos(2*π*x*f) - cos(π*t) ]
   
   (see e.g. https://en.wikipedia.org/wiki/List_of_periodic_functions)
   
   with a phase offset added in
   
   Here f is the frequency (1/iterations) and x is the advancement (iterations)
   
   H is the Heaviside step function
   
   The unique parameter of the function is the pulse duration, t, specified as a fraction of the wave period (iterations)
   
   The first bracketed term is a function of x and produces a cosine wave of the desired frequency.
   
   The second bracketed term produces a constant baseline offset, making the function positive, only for the desired pulse duration t.
   
   Application  of the Heaviside function coverts this into a rectangular pulse train
   
   (N.B. This function returns a value of either 1 or -1, wheras exact definition in terms of the Heaviside would return either 0 or 1)
   
   An additional complexity is that the function may be "chirped" - with a linear change in either the period or the frequency as the function advances.
   
   To achieve this we need to integrate the function describing the change in frequency, and use that as an argument to the cosinusoidal function.
   
   The necessary integral must be computed and supplied by the calling function
   
   */
  
  double output = std::cos(phase + clipper::Util::twopi()*integral) - std::cos( clipper::Util::pi()*pulse_duration);
  
  
  if(output >  0.0) {output =  1.0;}
  else              {output = -1.0;}
  
  return output;
  
}

bool ipa_functions::compute_consensus_mask(consensus_mask_settings settings,
                                           std::ostream& outlog,
                                           std::string out_directory,
                                           std::vector<std::string>& concensus_list)
{
  /*
   This function use to be its own independant c++ program, and has been quickly adapted for use within the library. 
   It is overlong and needs tidying up
   On the plus side, it's simple, it works, and it doesn't need to be integrated into the threading process like phase clustering
   */
  if (settings.verbose == true)
  {
    outlog << "\nSettings Input:\n";
    outlog << "min_pts: " << settings.min_pts << std::endl;
    outlog << "epsilon: " << settings.epsilon << std::endl; 
    outlog << "epsilon_is_absolute: " << settings.epsilon_is_absolute << std::endl; 
    outlog << "distance_measure: " << settings.distance_measure << std::endl;
    outlog << "erase_islands: " << settings.erase_islands << std::endl;
    outlog << "fill_voids: " << settings.fill_voids << std::endl;
    outlog << "erase_islands_threshold: " << settings.erase_islands_threshold << std::endl;
    outlog << "fill_voids_threshold: " << settings.fill_voids_threshold << std::endl;
    outlog << "input_known_mask: " << settings.input_known_mask << std::endl;
    outlog << "input_filename_known_mask: " << settings.input_filename_known_mask << std::endl;
    outlog << "input_Fourier_amplitudes: " << settings.input_Fourier_amplitudes << std::endl;
    outlog << "input_filename_target_data: " << settings.input_filename_target_data << std::endl;
    outlog << "target_col_fo: " << settings.target_col_fo << std::endl;
    outlog << "apodization: " << settings.apodization << std::endl;
    outlog << "apodization_limit: " << settings.apodization_limit << std::endl;
    outlog << "shift_to_known_origin: " << settings.shift_to_known_origin << std::endl;
    outlog << "allow_inversion: " << settings.allow_inversion << std::endl;
    outlog << "fraction_solvent_target: " << settings.fraction_solvent_target << std::endl;
    outlog << "verbose: " << settings.verbose << std::endl;
    
    outlog << "\nThe mask files that will be analyzed: " << std::endl;
    for (int i = 0; i < settings.file_list.size(); i ++ )
    {
      outlog << settings.file_list[i] << std::endl;
    }
  }
  
  clipper::HKL_data_base::HKL_reference_index ih;
  clipper::Xmap_base::Map_reference_index ix;
  
  clipper::CCP4MAPfile map_out;
  
  // Some Rudimentary checks on input
  
  if ( (settings.erase_islands_threshold < 0.0) or ( settings.erase_islands_threshold > 1.0) )
  {
    outlog << "Threshold value for island erasure should be > 0 and < 1\n";
    return false;
  }
  
  if ( (settings.fill_voids_threshold < 0.0) or ( settings.fill_voids_threshold > 1.0) )
  {
    outlog << "Threshold value for void filling should be > 0 and < 1\n";
    return false;
  }
  
  if ( (settings.fraction_solvent_target < 0.0)  )
  {
    outlog << "Fraction solvent has not been specified - this is compulsory input \n";
    return false;
  }

  
  if ( settings.epsilon_is_absolute && ((settings.epsilon >= 1.0) || (settings.epsilon <= 0.0)) )
  {
    outlog << "Absolute value for ε should be > 0 and < 1\n";
    return false;
  }
 
  if ( !settings.epsilon_is_absolute && ((settings.epsilon >= 100.0) || (settings.epsilon <= 0.0)) )
  {
    outlog << "Value for ε expressed as a percentile of the k-distance distribution should be  > 0 and < 100\n";
    return false;
  }
  
  if (settings.min_pts < 2)
  {
    outlog << "Minimum number of points cannot be < 2\n";
    return false;
  }
  
  outlog << "Number of mask files to be analyzed: " << settings.file_list.size() << "\n\n";
    
  std::vector< clipper::Xmap<int> > input_masks; // create a vector of integer valued map objects to hold the input masks
  clipper::Xmap<int>   work_mask;    // create a integer valued map object to hold the working mask
  
  
  // Read the mask files one by one
  // Need to implement some checks for mask existence and consistency here
  
  for(size_t t = 0; t < settings.file_list.size(); t++) // what is a size_t?!
  {
    outlog << "Reading mask: " << t << " " << settings.file_list[t] << std::endl;
    clipper::CCP4MAPfile maskin;
    maskin.open_read( settings.file_list[t] );
    maskin.import_xmap( work_mask ); // read in the current mask
    maskin.close_read();
    input_masks.push_back(work_mask); // augment the vector
  }
  
  bool space_group_is_achiral = !(util::space_group_is_chiral(input_masks[0].spacegroup().spacegroup_number()));
  
  // Get the Fourier amplitudes for the target structure.
  // TODO Get the Fourier amplitudes direct from settings rather than reading them here in standalone fashion
  // This is just a throwback to the origins of this function as a independent program
  
  clipper::HKL_info hkl_target; // Uninitialized clipper::HKL_info object to hold Spacegroup, Cell etc for the target structure
  
  outlog << "\nReading Fourier Amplitudes for the target structure" <<  std::endl;
  
  clipper::CCP4MTZfile mtzin; // create a CCP4MTZfile object
  mtzin.open_read( settings.input_filename_target_data ); // open the mtz file for reading
  mtzin.import_hkl_info( hkl_target ); // read sg, cell, reso
  clipper::HKL_data<clipper::data64::F_sigF>  meas_fo ( hkl_target ); // create a data object of type f + sigf for the measured Structure Factor amplitudes
  clipper::HKL_data<clipper::data64::F_sigF>  apod_fo ( hkl_target ); // create a data object of type f + sigf for the measured Structure Factor amplitudes, following apodization
  
  if (settings.target_col_fo != "NONE") mtzin.import_hkl_data( meas_fo, settings.target_col_fo ); // Here are the structure factor amplitudes  ...
  mtzin.close_read(); // close the mtz file after reading.
  
  // Report some basic information
  
  outlog << "\n TARGET STRUCTURE -- CRYSTALLOGRAPHIC INFO:\n\n";
  
  outlog << "Spacegroup " << hkl_target.spacegroup().symbol_xhm() << std::endl;
  outlog << "Cell Dimensions " << hkl_target.cell().format() << std::endl;
  outlog << "Cell Volume (Angstroms^3) "  << hkl_target.cell().volume() << std::endl;
  
  outlog << "\n TARGET STRUCTURE -- DATA INFO:\n\n";
  
  float target_min_resn = 1.0/sqrt(meas_fo.invresolsq_range().min());
  float target_max_resn = 1.0/sqrt(meas_fo.invresolsq_range().max());
  
  outlog << " Resolution range " <<  target_min_resn << " - " << target_max_resn << std::endl;
  outlog << " Number of data - possible " << hkl_target.num_reflections() << " - observed " << meas_fo.num_obs() << std::endl;
  
  // Create the weights to be used during correlation of Fourier amplitudes
  // These are either Gaussian or unitary
  
  clipper::HKL_data<clipper::data64::Phi_fom> apodization_weights(hkl_target);
  
  float sigma = 1.0/(2.4477*settings.apodization_limit);
  float sigma_squared = std::pow(sigma , 2);
  
  for ( ih = hkl_target.first(); !ih.last(); ih.next() )
  {
    if (settings.apodization)
    {
      float s_squared = ih.invresolsq();
      float weight = exp(-s_squared/(2*sigma_squared));
      
      apodization_weights[ih].phi() = 0.0;
      apodization_weights[ih].fom() = weight;
      apod_fo[ih].f()    =  weight*meas_fo[ih].f();
      apod_fo[ih].sigf() =  weight*meas_fo[ih].sigf();
    }
    else
    {
      apodization_weights[ih].phi() = 0.0;
      apodization_weights[ih].fom() = 1.0;
      apod_fo[ih].f()    =  meas_fo[ih].f();
      apod_fo[ih].sigf() =  meas_fo[ih].sigf();
    }
  }
  
  // If an known mask is input, get it ready for comparison with the cluster prototypes
  
  clipper::Xmap<int> known_mask(input_masks[0].spacegroup(), input_masks[0].cell(), input_masks[0].grid_sampling() );  // create a integer valued map object to hold the mask
  
  if (settings.input_known_mask)
  {
    outlog << "\nA known mask is input \n"  << std::endl;
    
    clipper::Xmap<float> input_target_mask_known; 
    
    // Read the mask
    outlog << "\nReading the known mask\n\n";
    clipper::CCP4MAPfile mapin1;
    mapin1.open_read(settings.input_filename_known_mask);
    mapin1.import_xmap(input_target_mask_known);   // import known target_mask
    mapin1.close_read();
    
    // Check that the input mask has the correct symmetry
    // Should also really check cell dimensions and extent ...
    
    outlog << "Grid Dimensions of the mask: nu " << input_target_mask_known.grid_sampling().nu() << " nv " << input_target_mask_known.grid_sampling().nv() << " nw " << input_target_mask_known.grid_sampling().nw() << std::endl;
    outlog << "Symmetry of the mask:" << input_target_mask_known.spacegroup().symbol_hm() << std::endl;
    
    if (input_target_mask_known.spacegroup().spacegroup_number() != input_masks[0].spacegroup().spacegroup_number())
    {
      outlog << "\nSymmetry mismatch between the input masks" << std::endl; return false;
    }
    else
    {
      outlog << "\nLooks okay. Proceeding ..." << std::endl;
    }
    
    // Interpolate the mask onto the same grid as the original input masks and the prototypes
    
    outlog << "Performing linear interpolation of the known mask onto the required grid" << std::endl;
    
    // nb input_target_mask_known is real valued, known_mask is integer valued
    for ( ix = known_mask.first(); !ix.last(); ix.next() )
    {
      clipper::Coord_frac cf = ix.coord().coord_frac(known_mask.grid_sampling()); // get the fractional coordinates of the current grid point
      known_mask[ix] =  lroundf( input_target_mask_known.interp<clipper::Interp_linear>( cf ) );
    }
  }
  
  std::vector< std::vector<float> > distance(input_masks.size(), std::vector<float>(input_masks.size(), 0.0)); // a matrix (vector-of-vectors) to hold the pairwise distance measure, initialized to zero
  std::vector< std::vector<int> > inversion_flags(input_masks.size(), std::vector<int>(input_masks.size(), 0)); // a matrix (vector-of-vectors) to hold flags indicating if masks need to be inverted to achieve maximal agreement, initialized to zero
  clipper::Coord_frac origin_shift( 0.0, 0.0, 0.0 );
  clipper::Coord_frac origin_shift_inv( 0.0, 0.0, 0.0 );
  std::vector< std::vector< clipper::Coord_frac > > shift_coords(input_masks.size(), std::vector< clipper::Coord_frac > (input_masks.size(), origin_shift) ); // a matrix (vector-of-vectors) of clipper data objects to hold the origin shifts required to maximize map/phase agreement
  
  outlog << "\nCalculating pairwise correlations of the input masks" << std::endl;
  outlog << "Distance will be calculated using distance measure " << settings.distance_measure << std::endl;
  if      (space_group_is_achiral &&  settings.allow_inversion) {outlog << "Distances will be assessed with allowance for mask inversion\n" <<  std::endl;}
  else if (space_group_is_achiral && !settings.allow_inversion) {outlog << "Distances will be assessed without allowance for mask inversion\n" <<  std::endl;}
  
  // Note that we don't bother to evaluate the diagonal of the distance matrix (the self-distances)  ...
  
  int n_masks = input_masks.size();
  for(int i = 0; i < n_masks; i++)
  {
    /*
    clipper::Xmap<int> integer_mask_i(input_masks[0].spacegroup(), input_masks[0].cell(), input_masks[0].grid_sampling() );
    for ( ix = integer_mask_i.first(); !ix.last(); ix.next() )
    {
      integer_mask_i[ix] = roundf(input_masks[i][ix]);
    }
    */
    for(int j = i+1; j < n_masks; j++)
    {
      /*
      clipper::Xmap<int> integer_mask_j(input_masks[0].spacegroup(), input_masks[0].cell(), input_masks[0].grid_sampling() );
      for ( ix = integer_mask_j.first(); !ix.last(); ix.next() )
      {
        integer_mask_j[ix] = roundf(input_masks[j][ix]);
      }
      */
      
      bool compute_real_space_agreement = true;

      float mask_cc=0.0;
      float smc = 0.0;
      float cc_smc = 0.0;
      float mcc = 0.0;
      
      ipa_functions::compute_mask_agreement (input_masks[i], input_masks[j], origin_shift, compute_real_space_agreement, mask_cc, smc, cc_smc, mcc, outlog);
      
      outlog << i << " " << j << " " << "Origin shift to maximize agreement between the masks: " << origin_shift.format() << " Matthews Correlation coefficient between masks: " << mcc <<  std::endl;

      float mask_cc_inv;
      float smc_inv;
      float cc_smc_inv;
      float mcc_inv = -1.0; // Set to -1 here, so that if inversion is not permitted,  mcc must be > mcc_inv 


      if(space_group_is_achiral && settings.allow_inversion)
      {
        outlog << "Repeating with inversion of mask: " << j <<  std::endl;
      
        clipper::Xmap<int> inverted_mask_j(input_masks[0].spacegroup(), input_masks[0].cell(), input_masks[0].grid_sampling() );
        inverted_mask_j = ipa_functions::invert_mask(input_masks[j]);
        ipa_functions::compute_mask_agreement (input_masks[i], inverted_mask_j, origin_shift_inv, compute_real_space_agreement, mask_cc_inv, smc_inv, cc_smc_inv, mcc_inv, outlog);

        outlog << i << " " << j << " " << "Origin shift to maximize agreement between the masks: " << origin_shift_inv.format() << " Matthews Correlation coefficient between masks: " << mcc_inv <<  std::endl;

      }
      
      // Compute the distances and store distances, origin shifts and inversion flags
      
      if (mcc > mcc_inv) // Correlation is maximized without inversion (or inversion is not possible)
      {
        distance[i][j] = ipa_functions::correlation_based_distance_measure(mcc, settings.distance_measure);
        shift_coords[i][j]     = origin_shift; 
        inversion_flags[i][j] = 0; // Indicates no inversion is required to maximize agreement. This is the default, but just for safety we set here again
      }
      else // Correlation is maximized with inversion 
      {
        distance[i][j] = ipa_functions::correlation_based_distance_measure(mcc_inv, settings.distance_measure);
        shift_coords[i][j]     = origin_shift_inv; 
        inversion_flags[i][j] = 1; // Indicates inversion is required to maximize agreement
      }
      // Fill in the values across the matrix diagonal.
      
      if (i != j)
      {
        distance[j][i] = distance[i][j];
        shift_coords[j][i] = -1.0*shift_coords[i][j];
        inversion_flags[j][i] = inversion_flags[i][j];
      }

      outlog << "distance: " << distance[i][j];
      
      if (space_group_is_achiral && settings.allow_inversion)
      {
        if (inversion_flags[i][j] == 0) {outlog << " was obtained without mask inversion\n" << std::endl;}
        else                            {outlog << " was obtained with mask inversion\n" << std::endl;}
      }
      else
      {
        outlog << "\n" << std::endl;
      }
            
    } // end loop on mask j
  } // end loop on mask i
  
  
  // Compute "k-distances" and "k-average-distances"
  // This section is experimental and can be tidied up once we figure out what works.

  outlog << "\nCalculating k-distances         = the distance between each mask and its k'th closest neighbour" << std::endl;
  outlog << "      and     k-average-distances = the mean distance between each mask and its k closest neighbours" << std::endl;
  outlog << "\nwhere k = MinPoints" << std::endl;

  
  std::vector<float> k_distances(input_masks.size(), 0.0);
  std::vector<float> k_average_distances(input_masks.size(), 0.0);
  
  //const int k = minPoints;
  
  for(int i = 0; i < n_masks; i++) // outer loop over all input masks
  {
    outlog << "\nMask " << i << std::endl;
    
    int counter =0;

    std::vector<float> current_distances((input_masks.size()-1), 0.0);    
    for(int j = 0; j < n_masks; j++) // inner loop over all input masks 
    {
      if (i != j) // exclude self-agreement 
      {
        current_distances[counter] = distance[i][j];
        counter += 1;
      }
    }

    // Now we have the relevant distances for mask i, sort them in ascending order
    
    std::sort(current_distances.begin(), current_distances.begin() + current_distances.size(), std::less<float>());
    
    /*
    outlog << "\nDistances to all other masks, in rank increasing order" << std::endl;
    for(size_t j = 0; j < current_distances.size(); j++)
    {
      outlog << current_distances[j] << std::endl;
    }
    */
    
    // Compute the k-distance and the k-average-distance
    
    k_distances[i] = current_distances[settings.min_pts-1]; 
  
    for(size_t j = 0; j < settings.min_pts; j++)
    {
      k_average_distances[i] += current_distances[j];
    }
    k_average_distances[i] /= float(settings.min_pts);
    
    
    
    outlog << "\nk-distance (k = " << settings.min_pts << " ) is " << k_distances[i] << std::endl;
    outlog << "k-average-distance (k = " << settings.min_pts << " ) is " << k_average_distances[i] << std::endl;
    
  }
  
  // Now we have all the k-distances and k-average-distances, sort them, compute a simple curvature measure, and output the results
  
  std::vector<float> theta_k_distances(input_masks.size(), 0.0);
  std::vector<float> theta_k_average_distances(input_masks.size(), 0.0);
  

  std::sort(k_distances.begin(), k_distances.begin() + k_distances.size(), std::less<float>());

  outlog << "\nk-distances for all masks, in rank increasing order" << std::endl;
  outlog << "and correspondent curvature measure, Θ(°), of the discrete plot\n" << std::endl;
   
  // Simple minded procedure for trying to locate the elbows in the plot, which would indicate what values of epsilon might be appropriate
  
  for(size_t j = 1; j < (k_distances.size() - 1); j++)
  {
    // Define θj to be the signed angle between the vector a = pj − pj−1 and the vector b = pj+1 − pj.
    // This seems like a reasonable measure of curvature ... Note that when the two vectors are colinear then θj must be zero
    
    // compute the components of two vectors a and b describing the relevant line segments
    double xa = 1.0;
    double ya = k_distances[j]-k_distances[j-1];;
    double xb = 1.0; 
    double yb = k_distances[j+1]-k_distances[j];
    
    // Compute the (directed) angle between vector a and vector b which is a measure of the "curvature" of the plot.
    // This is the angle of vector b with the x axis - the angle of vector a with the x axis 
    // Points with very large positive theta will correspond to the "elbows" in the plot.

    theta_k_distances[j] =  atan2(yb, xb) - atan2(ya, xa);
    theta_k_distances[j] *= 180.0/clipper::Util::pi(); // Degrees are easier ?
    
  }

  for(size_t j = 0; j < k_distances.size(); j++)
  {
    if ( (j == 0) || (j == (k_distances.size()-1)) )
    {
      outlog << j << " " << k_distances[j] << std::endl;      
    }
    else
    {
      outlog << j << " " << k_distances[j] <<  " " << theta_k_distances[j] << std::endl;
    }
      
  }

     
  std::sort(k_average_distances.begin(), k_average_distances.begin() + k_average_distances.size(), std::less<float>());

  outlog << "\nk-average-distances for all masks, in rank increasing order" << std::endl;
  outlog << "and correspondent curvature measure, Θ(°), of the discrete plot\n" << std::endl;
     
  // Simple minded procedure for trying to locate the elbows in the plot, which would indicate what values of epsilon might be appropriate
  
  for(size_t j = 1; j < (k_average_distances.size() - 1); j++)
  {
    // Define θj to be the signed angle between the vector a = pj − pj−1 and the vector b = pj+1 − pj.
    // This seems like a reasonable measure of curvature ... Note that when the two vectors are colinear then θj must be zero
    
    // compute the components of two vectors a and b describing the relevant line segments
    double xa = 1.0;
    double ya = k_average_distances[j]-k_average_distances[j-1];;
    double xb = 1.0; 
    double yb = k_average_distances[j+1]-k_average_distances[j];
    
    // Compute the (directed) angle between vector a and vector b which is a measure of the "curvature" of the plot.
    // This is the angle of vector b with the x axis - the angle of vector a with the x axis 
    // Points with very large positive theta will correspond to the "elbows" in the plot.
    
    theta_k_average_distances[j] = atan2(yb, xb) - atan2(ya, xa);
    theta_k_average_distances[j] *= 180.0/clipper::Util::pi(); // Degrees are easier ?
    
  }
  
  for(size_t j = 0; j < k_average_distances.size(); j++)
  {
    if ( (j == 0) || (j == (k_average_distances.size()-1)) )
    {
      outlog << j << " " << k_average_distances[j] << std::endl;      
    }
    else
    {
      outlog << j << " " << k_average_distances[j] <<  " " << theta_k_average_distances[j] << std::endl;
    }
      
  }
  
  /*
  
  // Put distance measure values in rank order, to get a feeling for the best and worst agreement, ignoring self-agreement (the diagonal terms)
  // Number of triangular matrix entries in a 𝑁×𝑁 matrix - without the diagonal it is N*(N-1)/2 - with the diagonal it is N*(N+1)/2
  
  std::vector<float> dm_values(( (input_masks.size()  * (input_masks.size() - 1) ) / 2 ), 0.0);
  
  int ntotal = -1;
  for(size_t i = 0; i < input_masks.size(); i++)
  {
    for(size_t j = i+1; j < input_masks.size(); j++)
    {
      ntotal += 1;
      dm_values[ntotal] = distance[i][j];
    }
  }
  
  //outlog << dm_values.size() << std::endl;
  
  std::sort(dm_values.begin(), dm_values.begin() + dm_values.size(), std::less<float>());
  
  outlog << "\nPairwise distances in rank increasing order, excluding self-agreement " << std::endl;
  
  for(size_t i = 0; i < dm_values.size(); i++)
  {
    outlog << i << " " << dm_values[i]  << std::endl;
  }
  
  size_t size = dm_values.size();
  
  outlog << "\nMinimum distance: " << dm_values[0] << std::endl;
  outlog << "Maximum distance: " << dm_values[size-1] << std::endl;
  
  float median;
  
  if ( (size % 2) == 0) // even number of values
  {median = (dm_values[(size/2)] + dm_values[(size/2)-1])/2.0;}
  
  else                  // odd number of values, note that integer division in c++ returns the floor
  {median =  dm_values[(size/2)];}
  
  outlog << "Median distance: " << median << std::endl;
  
  */
  
  double threshold_absolute;
  double threshold_as_percentile;  
   
  if (settings.epsilon_is_absolute)
    {
     // threshold distance for clustering is specified in absolute terms. Work out the approximate corresponding percentile in the k-distance distribution, for subsequent reporting

      threshold_absolute = settings.epsilon;
      
      /*
      int i=-1;
      do
      {
        i += 1;
      }
      while ( (dm_values[i] < settings.epsilon) && (i <= size) );
    
      threshold_as_percentile = 100*(float(i)/size);
      */
      
      int i=-1;
      do
      {
        i += 1;
      }
      while ( (k_distances[i] < settings.epsilon) && (i < k_distances.size() ) );
    
      threshold_as_percentile = 100*(float(i)/k_distances.size());

      
    }  
  else
    {
      // threshold distance for clustering is specified as a percentile of the k-distance distribution. Work out the corresponding absolute distance for the subsequent clustering operation
    
      threshold_as_percentile = settings.epsilon;

      int t_index = static_cast<int>( (settings.epsilon/100)*k_distances.size() ) - 1;
      
      // Some bounds checks, for safety
      if (t_index < 0)                        {t_index = 0;} 
      if (t_index > (k_distances.size() - 1)) {t_index = (k_distances.size() - 1);}
            
      threshold_absolute = k_distances[t_index];

    }
   
  struct cluster_data // Only ever used here so we're going to keep in scope.
  {
    float percentile;
    float  epsilon;
    int cluster_id;
    int n_members;
    int n_connected_sets_eliminated_protein;
    int n_connected_sets_retained_protein;
    int n_connected_sets_eliminated_solvent;
    int n_connected_sets_retained_solvent;
    float fraction_solvent_pre_edit;
    float fraction_solvent_post_edit;
    float mean_mcc_of_cluster_members_with_known;
    float mean_mcc_of_inverted_cluster_members_with_known;
    float smc_with_known;
    float cc_smc_with_known;
    float mcc_with_known;
    float smc_inverted_with_known;
    float cc_smc_inverted_with_known;
    float mcc_inverted_with_known;
    float cc_with_measured_amplitudes;
    float rf_with_measured_amplitudes;
  };
  
  std::vector<cluster_data> all_prototypes;
  
  outlog << "\nPerforming density based clustering (DBSCAN) on the input masks." << std::endl;

    double epsilon = threshold_absolute;
    double percentile  = threshold_as_percentile;    
    
    outlog << "Threshold distance (ε) expressed as an absolute distance: " << epsilon << std::endl;
    outlog << "Threshold distance (ε) expressed as a percentile of the k-distance distribution: " << percentile  << std::endl;
    outlog << "Parameter minPoints: " << settings.min_pts  << std::endl;
    outlog << std::endl;    

    int n_clusters;
    std::vector<int> cluster_id;
    
    util::db_scan(distance, epsilon, settings.min_pts, n_clusters, cluster_id, outlog);
    
    outlog << "\nNumber of clusters: " << n_clusters << std::endl;
    outlog << "\nMask  Cluster ID " << std::endl;
    for(int i = 0; i < n_masks; i++) {outlog << i << "      " << cluster_id[i] << std::endl;}
    outlog << "\nA Cluster ID of -2 indicates the mask is classified as noise" << std::endl;
    outlog << std::endl;
    
    if (n_clusters > 0)
    {
      outlog << "\nNow Generating the consensus masks (cluster prototypes)" << std::endl;
      outlog << std::endl;
    }
    
    clipper::Xmap<int> cluster_prototype(input_masks[0].spacegroup(), input_masks[0].cell(), input_masks[0].grid_sampling());
    clipper::Xmap<int> cluster_prototype_edited(input_masks[0].spacegroup(), input_masks[0].cell(), input_masks[0].grid_sampling());
    clipper::Xmap<int> cluster_prototype_edited_shifted(input_masks[0].spacegroup(), input_masks[0].cell(), input_masks[0].grid_sampling());
    clipper::Xmap<int> cluster_prototype_edited_inverted(input_masks[0].spacegroup(), input_masks[0].cell(), input_masks[0].grid_sampling());
    clipper::Xmap<int> cluster_prototype_edited_inverted_shifted(input_masks[0].spacegroup(), input_masks[0].cell(), input_masks[0].grid_sampling());
    
    for(int ic = 0; ic < n_clusters; ic++) // loop over all the clusters
    {
      
      clipper::Xmap<float> mask_summation; 
      mask_summation.init( input_masks[0].spacegroup(), input_masks[0].cell(), input_masks[0].grid_sampling() );
      
      cluster_data current_prototype;
      
      outlog << "Generating consensus mask for cluster: " << ic << std::endl;
      
      outlog << "The Members of this cluster: ";
      
      // Count the number of members in each cluster and list them
      current_prototype.n_members = 0;
      for(size_t im = 0; im < input_masks.size(); im++) // loop over all the input masks
      {
        if (cluster_id[im] == ic)
        {
          outlog << im << ", ";
          current_prototype.n_members += 1;
        }
      }
      outlog << std::endl;
      
      // DBSCAN can produce a "cluster" consisting of a single core point. Exclude these "clusters" from analysis
      
      if (current_prototype.n_members == 1)
      {
        outlog << "This cluster consists of only a single core point. Excluding from subsequent analysis" << std::endl;
        for(size_t im = 0; im < input_masks.size(); im++) // loop over all the input masks
        {
          if (cluster_id[im] == ic) {cluster_id[im] = -2;}
        }
        continue; // this will take us to the end of the loop over all clusters
      }
      
      current_prototype.percentile = percentile;
      current_prototype.epsilon = epsilon;
      current_prototype.cluster_id = ic;

      int n_members_identified =0; // a counter for the number of masks in this cluster
      
      for(size_t im = 0; im < input_masks.size(); im++) // loop over all input masks to sum the members of the cluster
      {
        int im_ref; // identifies the reference mask
        
        if (cluster_id[im] == ic) // if true, the mask is is a member of the current cluster
        {
          n_members_identified +=1; // augment the membership counter
          if (n_members_identified == 1)  {im_ref = im;} // the first mask in the cluster is the reference, which is used to define both the origin and the hand

          if (settings.verbose) {outlog << "Processing mask: " << im << " for which the inversion flag is: " << inversion_flags[im_ref][im] << " and the required origin shift is: " << shift_coords[im_ref][im].format() << std::endl;}
          
          clipper::Xmap<int> inverted_mask(input_masks[0].spacegroup(), input_masks[0].cell(), input_masks[0].grid_sampling()); // a place to store an inverted copy of the mask
          if (inversion_flags[im_ref][im] == 1) //  The mask needs to be inverted to give it the same hand as the reference 
          {
            inverted_mask = ipa_functions::invert_mask(input_masks[im]); // Now generate the inverted mask
          }

          
          // Loop over all points in the mask and augment the sum over all members of the cluster
          
          for ( clipper::Xmap_base::Map_reference_index ix = mask_summation.first(); !ix.last(); ix.next() )
          {
            clipper::Coord_frac cf = ix.coord().coord_frac(cluster_prototype.grid_sampling()); // get the fractional coordinates of the current grid point
            clipper::Coord_frac cf_shifted = cf - shift_coords[im_ref][im]; // apply the needed origin shift
            
             if (inversion_flags[im_ref][im] == 1) // we need to use the inversion of the current mask im
             {
               mask_summation[ix] += inverted_mask.interp<clipper::Interp_linear>( cf_shifted );           

             }
             else // we need to use the original mask im
             {
               mask_summation[ix] += input_masks[im].interp<clipper::Interp_linear>( cf_shifted );           
             }
          }
        }
        
      }
      
      // Now we have summed over all masks in the cluster, allowing for origin shifts and inversion 
      // Operate on the sum to generate the final mask (plus a copy for editing) 
      
      for ( clipper::Xmap_base::Map_reference_index ix = cluster_prototype.first(); !ix.last(); ix.next() )
      {
        if (float(mask_summation[ix])/float(n_members_identified) >= 0.5 ) // the final output is the mode of the input
        {
          // protein region
          cluster_prototype[ix] = 1;
          cluster_prototype_edited[ix] = 1; // a copy that will be edited
        }
        else
        {
          // solvent region
          cluster_prototype[ix] = 0;
          cluster_prototype_edited[ix] = 0; // a copy that will be edited
        }
        
      }      
      
      // Output the cluster prototype
      
      std::string io_map_filename;
      
      // Construct the output filename
      // First convert the floating point epsilon into a string.
      // Use std::stringstream Class and str() Method to do this
      std::stringstream sstream;
      sstream << epsilon;
      std::string eps_as_str = sstream.str();

      eps_as_str = eps_as_str.substr( eps_as_str.find(".") + 1);  // Retain only what follows the decimal point.

      // Build the filename (added appropriate out_directory to this)
      //outlog << "DEBUG: " << out_directory << std::endl;
      io_map_filename = out_directory + "min_points_" + std::to_string(settings.min_pts) + "_epsilon_0p" + eps_as_str.substr(0,5) + "_prototype_" + std::to_string(ic) + ".ccp4";
      
      map_out.open_write( io_map_filename );
      map_out.export_xmap(cluster_prototype); // write ouput_mask
      map_out.close_write();
      
      /*
       There are two options here, write the file name to a file_list.txt - similar to how we do the outputs of envelope_determination, however I think opting for an internal solution with no file_list.txt file makes more sense:
       a) You're only going to be running this once, or if you do run it multiple times, its likely that its because the masks input/output has changed, and you will be inputting the desired file name directly in the experiment.params file.
       b) You only need to pass the current experiments output if we're trying to run all jobs at once, and that will be because there are no other mask concensus's to consider anyway.

       TODO ... Related point - when the procedure generates more than one mask, we can devise some criteria to put the masks in order of likely correctness, before being fed into phase determination
       For example, correct masks are consistent with the specified solvent fraction, incorrect masks often aren't (the averaging process means that the consensus mask doesn't have to have the same solvent fraction as the imputs)
 
       */
      
      
      // calculate the fraction solvent before envelope editing
	  
      float pts_total = 0;
      float pts_solvent = 0;
      
      for ( ix = cluster_prototype.first(); !ix.last(); ix.next() )
      {
        float multiplicity = cluster_prototype.multiplicity( ix.coord());
        pts_total += (1.0 / multiplicity);
        if (cluster_prototype[ix] == 0)
        {
          pts_solvent += (1.0 / multiplicity);
        }
      }
      current_prototype.fraction_solvent_pre_edit = pts_solvent/pts_total;
      
      // Now perform Connectivity reinforcement on the cluster prototype - erasing islands and filling voids
      
      outlog << "\nPerforming Connectivity reinforcement on the cluster prototype" << std::endl;
      
      int n_connected_sets_eliminated_protein;
      int n_connected_sets_retained_protein;
      int n_connected_sets_eliminated_solvent;
      int n_connected_sets_retained_solvent;
      
      ipa_functions::enforce_mask_connectivity (settings.erase_islands, settings.fill_voids, settings.erase_islands_threshold, settings.fill_voids_threshold, cluster_prototype_edited,
                                                n_connected_sets_eliminated_protein, n_connected_sets_retained_protein, n_connected_sets_eliminated_solvent, n_connected_sets_retained_solvent, outlog, settings.verbose);
      
      current_prototype.n_connected_sets_eliminated_protein = n_connected_sets_eliminated_protein;
      current_prototype.n_connected_sets_retained_protein = n_connected_sets_retained_protein;
      current_prototype.n_connected_sets_eliminated_solvent = n_connected_sets_eliminated_solvent;
      current_prototype.n_connected_sets_retained_solvent = n_connected_sets_retained_solvent;
      
      // recalculate the fraction solvent after envelope editing
      
      pts_total = 0;
      pts_solvent = 0;
      
      for ( ix = cluster_prototype_edited.first(); !ix.last(); ix.next() )
      {
        float multiplicity = cluster_prototype_edited.multiplicity( ix.coord());
        pts_total += (1.0 / multiplicity);
        if (cluster_prototype_edited[ix] == 0)
        {
          pts_solvent += (1.0 / multiplicity);
        }
      }
      current_prototype.fraction_solvent_post_edit = pts_solvent/pts_total;
      
      
      
      // Output the cluster prototype after editing
      
      // Build the filename
      io_map_filename = out_directory + "min_points_" + std::to_string(settings.min_pts) + "_epsilon_0p" + eps_as_str.substr(0,5) + "_prototype_" + std::to_string(ic) + "_edited.ccp4";
      
      map_out.open_write( io_map_filename );
      map_out.export_xmap( cluster_prototype_edited); // write ouput_mask
      map_out.close_write();
      
	  // Add this consensus envelope (cluster prototype) to the list of files to be used for phase determination
      concensus_list.push_back(io_map_filename); 
      
      // Output the inverse of the cluster prototype after editing
      
      if (space_group_is_achiral) // space group is invariant under change of hand
      {
        cluster_prototype_edited_inverted = ipa_functions::invert_mask(cluster_prototype_edited);
        
        // Build the filename
        io_map_filename = out_directory + "min_points_" + std::to_string(settings.min_pts) + "_epsilon_0p" + eps_as_str.substr(0,5) + "_prototype_" + std::to_string(ic) + "_edited_inverted.ccp4";
        
        map_out.open_write( io_map_filename );
        map_out.export_xmap(cluster_prototype_edited_inverted); // write out the inverted mask
        map_out.close_write();
                
      }
      
      // Compute the agreement between each member of the cluster and the known mask, if given
      
      
      if (settings.input_known_mask)
      {
        
        outlog << "\nComputing the agreement between the known mask and each member of cluster: " <<  ic <<  "\n" << std::endl;
        
        clipper::Xmap<int>  working_mask(input_masks[0].spacegroup(), input_masks[0].cell(), input_masks[0].grid_sampling());   // create an object to hold a integer valued copy of the mask
        
        bool compute_real_space_agreement = true;
        
        float mask_cc=0.0;
        float smc = 0.0;
        float cc_smc = 0.0;
        float mcc = 0.0;

        float mask_cc_inv=0.0;
        float smc_inv = 0.0;
        float cc_smc_inv = 0.0;
        float mcc_inv = 0.0;

       
        int n_members_identified = 0;
        float mean_mcc = 0.0;
        float mean_dist = 0.0;
        
        current_prototype.mean_mcc_of_cluster_members_with_known = 0;
        current_prototype.mean_mcc_of_inverted_cluster_members_with_known = 0;
        
        for(size_t im = 0; im < input_masks.size(); im++) // loop over all the masks in the dataset
        {
          if (cluster_id[im] == ic) // if true, the mask is a member of the current cluster
          {
            n_members_identified += 1;
            
            /*
            for ( ix = cluster_prototype.first(); !ix.last(); ix.next() )
            {
              working_mask[ix] = roundf(input_masks[im][ix]);
            }
            */
            
            if (settings.verbose == true) {outlog << "\nMask: " << im << std::endl;}
            clipper::Coord_frac origin_shift;
            ipa_functions::compute_mask_agreement (known_mask, input_masks[im], origin_shift, compute_real_space_agreement, mask_cc, smc, cc_smc, mcc, outlog, settings.verbose);
                        
            if (space_group_is_achiral) // space group is invariant under change of hand
            {
              if (settings.verbose == true) {outlog << "\nInverse of Mask: " << im << std::endl;}
              clipper::Xmap<int> inverted_mask(input_masks[0].spacegroup(), input_masks[0].cell(), input_masks[0].grid_sampling());
              
              inverted_mask = ipa_functions::invert_mask(input_masks[im]);
              ipa_functions::compute_mask_agreement (known_mask, inverted_mask, origin_shift_inv, compute_real_space_agreement, mask_cc_inv, smc_inv, cc_smc_inv, mcc_inv, outlog, settings.verbose);
            }

            mean_mcc += std::max(mcc, mcc_inv);
            mean_dist += ipa_functions::correlation_based_distance_measure(std::max(mcc, mcc_inv), settings.distance_measure);

          }
        }
        mean_mcc /= (float)n_members_identified;
        mean_dist /= (float)n_members_identified;
        
        current_prototype.mean_mcc_of_cluster_members_with_known = mean_mcc;
        
        outlog << std::endl;
        if (space_group_is_achiral) {outlog << "Allowing for inversion: " << std::endl;}
        outlog << "The mean correlation of the clustered masks with the known mask is: " << mean_mcc << std::endl;
        outlog << "The mean distance of the clustered masks from the known mask is: " << mean_dist << std::endl;
              
        // Compute the agreement of the cluster prototype with the known mask, if given
      
        current_prototype.smc_with_known = 0;
        current_prototype.cc_smc_with_known = 0;
        current_prototype.mcc_with_known = 0;
        current_prototype.smc_inverted_with_known = 0;
        current_prototype.cc_smc_inverted_with_known = 0;
        current_prototype.mcc_inverted_with_known = 0;
              
        
        if (settings.verbose == true) {outlog << "\nAgreement of the known mask with the cluster " << ic << " prototype " << std::endl;}

        clipper::Coord_frac origin_shift;

        if (settings.verbose == true) {outlog << "\nAgreement with the edited Prototype:\n";}
        ipa_functions::compute_mask_agreement (known_mask, cluster_prototype_edited, origin_shift, compute_real_space_agreement, mask_cc, smc, cc_smc, mcc, outlog, settings.verbose);
        
        current_prototype.smc_with_known = smc;
        current_prototype.cc_smc_with_known = cc_smc;
        current_prototype.mcc_with_known = mcc;
        
        if (settings.shift_to_known_origin) // write out the edited cluster prototype after applying the origin shift
        {
          cluster_prototype_edited_shifted = ipa_functions::shift_mask(cluster_prototype_edited, -origin_shift);
          // Build the filename
          io_map_filename = out_directory + "min_points_" + std::to_string(settings.min_pts) + "_epsilon_0p" + eps_as_str.substr(0,5) + "_prototype_" + std::to_string(ic) + "_edited_shifted.ccp4";
          map_out.open_write( io_map_filename );
          map_out.export_xmap(cluster_prototype_edited_shifted);
          map_out.close_write();
          
        }
        
        // should only do the following for the 43 achiral space groups. Inversion in the 22 chiral space groups can't be done as it changes the space group ...
        
        if (space_group_is_achiral) // space group is invariant under change of hand
        {
          outlog << "\nAgreement with the *inverse* of the edited Prototype: \n";
          ipa_functions::compute_mask_agreement (known_mask, cluster_prototype_edited_inverted, origin_shift, compute_real_space_agreement, mask_cc, smc, cc_smc, mcc, outlog, settings.verbose);
          
          current_prototype.smc_inverted_with_known = smc;
          current_prototype.cc_smc_inverted_with_known = cc_smc;
          current_prototype.mcc_inverted_with_known = mcc;
          
          if (settings.shift_to_known_origin) // write out the edited and inverted cluster prototype after applying the origin shift
          {
            
            // clipper::Xmap<int>  working_mask(input_masks[0].spacegroup(), input_masks[0].cell(), input_masks[0].grid_sampling());   // create an object to hold a integer valued copy of the mask
            
            cluster_prototype_edited_inverted_shifted = ipa_functions::shift_mask(cluster_prototype_edited_inverted, -origin_shift);
            // Build the filename
            io_map_filename = out_directory + "min_points_" + std::to_string(settings.min_pts) + "_epsilon_0p" + eps_as_str.substr(0,5) + "_prototype_" + std::to_string(ic) + "_edited_inverted_shifted.ccp4";
            map_out.open_write( io_map_filename );
            map_out.export_xmap( cluster_prototype_edited_inverted_shifted );
            map_out.close_write();
            
          }
          
        }
        else
        {
          current_prototype.smc_inverted_with_known = 0.0;
          current_prototype.cc_smc_inverted_with_known = 0.0;
          current_prototype.mcc_inverted_with_known = 0.0;
        }
        
        outlog << std::endl;
      }
      
	  
      // Compare Fourier amplitudes calculated from the mask with the measured Fourier amplitudes
      
      if (settings.verbose == true) {outlog << "\n Comparing Fourier amplitudes calculated from the envelope with the measured Fourier amplitudes: \n";}
      
      clipper::Xmap<int> real_valued_mask(input_masks[0].spacegroup(), input_masks[0].cell(), input_masks[0].grid_sampling());
      clipper::HKL_data<clipper::data64::F_phi> fmask (hkl_target); // create a data object of type  F + phi to hold Fourier coefficients of the mask
      
      for ( ix = cluster_prototype_edited.first(); !ix.last(); ix.next() )
      {
        real_valued_mask[ix] = (float)cluster_prototype_edited[ix];
      }
      
      real_valued_mask.fft_to (fmask); // FFT of the real-valued mask to get the Crystallographic structure factors of the mask
      
      float esd_threshold = 2.0;
      int res_bins = 25;
      bool use_test_set = false;
      double  correlation_coefficient;
      double  r_factor;
      double  correlation_coefficient_test_set;
      double  r_factor_test_set;
      double scale_factor;
      
      ipa_functions::simple_linear_scale(apod_fo,  apod_fo,  apodization_weights,
                                         fmask, use_test_set, esd_threshold, res_bins, scale_factor, correlation_coefficient, r_factor, correlation_coefficient_test_set, r_factor_test_set, outlog);
      
      current_prototype.cc_with_measured_amplitudes = correlation_coefficient;
      current_prototype.rf_with_measured_amplitudes = r_factor;
      	  
      // Store all the information about this prototype
      
      all_prototypes.push_back(current_prototype);
      
      
    } // end of loop over all clusters
    
    
  
  outlog << "\nSummmary of results" << std::endl;
  outlog << "     Percentile         DBSCAN ε      Cluster ID        #members sets eliminated   sets retained  sets eliminated   sets retained   ΔSolventfrac    ΔSolventfrac           SMC            cc-SMC            MCC         Inv SMC       Inv cc-SMC       Inv  MCC      cc with      r-factor with  " << std::endl;
  outlog << "                                                                    protein           protein        solvent           solvent        pre edit        post edit         with known      with known      with known    with known      with known      with known   measured |F|    measured |F|  " << std::endl;
  
  for ( int i = 0; i < all_prototypes.size(); i++ )
  {
    outlog <<
    std::setw(15) << all_prototypes[i].percentile << " " <<
    std::setw(15) << all_prototypes[i].epsilon << " " <<
    std::setw(15) << all_prototypes[i].cluster_id << " " <<
    std::setw(15) << all_prototypes[i].n_members << " " <<
    std::setw(15) << all_prototypes[i].n_connected_sets_eliminated_protein << " " <<
    std::setw(15) << all_prototypes[i].n_connected_sets_retained_protein << " " <<
    std::setw(15) << all_prototypes[i].n_connected_sets_eliminated_solvent << " " <<
    std::setw(15) << all_prototypes[i].n_connected_sets_retained_solvent << " " <<
    std::setw(15) << all_prototypes[i].fraction_solvent_pre_edit  - settings.fraction_solvent_target << " " <<
    std::setw(15) << all_prototypes[i].fraction_solvent_post_edit - settings.fraction_solvent_target << " " <<
    std::setw(15) << all_prototypes[i].smc_with_known << " " <<
    std::setw(15) << all_prototypes[i].cc_smc_with_known << " " <<
    std::setw(15) << all_prototypes[i].mcc_with_known << " " <<
    std::setw(15) << all_prototypes[i].smc_inverted_with_known << " " <<
    std::setw(15) << all_prototypes[i].cc_smc_inverted_with_known << " " <<
    std::setw(15) << all_prototypes[i].mcc_inverted_with_known << " " <<
    std::setw(15) << all_prototypes[i].cc_with_measured_amplitudes << " " <<
    std::setw(15) << all_prototypes[i].rf_with_measured_amplitudes << std::endl;
  }
  
  // Hack to produce output required for papers
  
  // open file for run summary
  
  /*
   std::ofstream filewrite("cluster_summmary.csv", std::ofstream::out);
   
   for ( int i = 0; i < all_prototypes.size(); i++ )
   {
   filewrite <<
   std::setw(15) << all_prototypes[i].epsilon << "," <<
   std::setw(15) << all_prototypes[i].cluster_id << "," <<
   std::setw(15) << all_prototypes[i].n_members << "," <<
   std::setw(15) << all_prototypes[i].mean_mcc_of_cluster_members_with_known << "," <<
   std::setw(15) << all_prototypes[i].n_connected_sets_retained_protein << "," <<
   std::setw(15) << all_prototypes[i].fraction_solvent_post_edit - fraction_solvent_target << "," <<
   std::setw(15) << all_prototypes[i].mcc_with_known << "," <<
   std::setw(15) << all_prototypes[i].mcc_inverted_with_known << std::endl;
   }
   
   filewrite.close();
   */
  
  // Program is terminating normally
  outlog << "Normal termination." << std::endl;
  return true;
  
}

void ipa_functions::compute_phase_total_difference(const clipper::HKL_data<clipper::data64::F_phi> &fp1, const clipper::HKL_data <clipper::data64::Phi_fom> &apodization_weights1,
                                                   const clipper::HKL_data<clipper::data64::F_phi> &fp2, const clipper::HKL_data <clipper::data64::Phi_fom> &apodization_weights2,
                                                   float &mean_distance, double& mean_cc, float& mean_phi_difference, bool centrics_only, std::ostream &outlog, bool verbose)
{
  
  clipper::HKL_data_base::HKL_reference_index ih;
  
  // Used to calculate the mean absolute phase difference.
  double mean_phase_diff = 0.0;
  double summed_phase_differences = 0.0;
  double countp = 0;

  // Used to calculate the distance travelled in the complex plane 
  double summed_distances = 0.0;
  double xn = 0;
  double yn = 0;
  double countd = 0;

  // Used to calculate the map correlation coefficient
  double cc_num = 0;
  double cc_denom1 = 0;
  double cc_denom2 = 0;
  
  // Now loop over the asymmetric unit and compute the agreement measures
  // Because there is no guarantee that the first and second data sets have the same hkl list, we have to access the 2nd data set by explicit lookup of hkl
  
  for ( ih = fp1.first(); !ih.last(); ih.next() )
  {
    if ( !fp1[ih].missing() && !fp2[ih.hkl()].missing() && !(!ih.hkl_class().centric() && centrics_only) )
    {
      
      const double w = 1.0 / ih.hkl_class().epsilon();   // Get the standard statistical weights
      
      // work out the weights used for apodization
      const double wa1 = apodization_weights1[ih].fom(); // Get the apodization weights for the first data set
      const double wa2 = apodization_weights2[ih.hkl()].fom(); // Get the apodization weights for the second data set
      const double wa = std::min(wa1,wa2); // The apodization scheme could be different for the two data sets. So take the minimum as the apodization weight
      // Note, due to the program convention, wa1 will always be the lowest weight, therefore, no operational need to deal with wa2 really.

      // precompute these quantities for speed.
      const double c_phase_diff = cos( fp1[ih].phi() - fp2[ih.hkl()].phi() );
      const double f1 = fp1[ih].f();
      const double f2 = fp2[ih.hkl()].f();
      
      //  the absolute phase difference, accounting for periodicity 
      const double d = acos ( c_phase_diff );

      summed_phase_differences += (w * wa * d); // Augment the summation
      countp += (w * wa * 1); // Count should be incremented by the weights  ...

        
      // Augment the sums required for computation of the map correlation coefficient.
      // For details of the correlation coefficient calculation see 
      // Lunin VY, Woolfson MM. Mean phase error and the map-correlation coefficient. Acta Crystallogr D. 1993 Nov 1;49(Pt 6):530–3. 
      // also Bailey GD, Hyun J-KK, Mitra AK, Kingston RL. J Mol Biol. 2012 Mar 30;417(3):212–23. Supplementary Methods, Section 3.
      //
      // Note that the statistical weights must be included here as we are summing over the asymmetric unit only.
      
      if ( !((ih.hkl().h() == 0) && (ih.hkl().k() == 0) && (ih.hkl().l() == 0)) ) // F(000) must be excluded from these summations to ensure mean densities are zero
      {
        cc_num    += w*f1 * w*f2 * c_phase_diff;
        cc_denom1 += w*f1 * w*f1;
        cc_denom2 += w*f2 * w*f2;
      }    

      /// calculate the distance moved in the complex plane
      
      double x1, x2, y1, y2;
      x1 = cos(fp1[ih].phi()) * fp1[ih].f();
      y1 = sin(fp1[ih].phi()) * fp1[ih].f();      
      x2 = cos(fp2[ih.hkl()].phi()) * fp2[ih.hkl()].f();
      y2 = sin(fp2[ih.hkl()].phi()) * fp2[ih.hkl()].f();

      xn += (x2 - x1);
      yn += (y2 - y1);

      summed_distances += w * sqrt(pow(xn,2) + pow(yn,2)); // // Augment the summation
      countd += (w * 1); // Count should be incremented by the weights  ...  
          
    }
  }
  
  mean_phase_diff = summed_phase_differences/countp;    // Compute the mean absolute phase difference 
  mean_phi_difference = float(clipper::Util::rad2d(mean_phase_diff));   // Convert result to degrees.
  
  mean_cc = cc_num / std::sqrt(cc_denom1*cc_denom2);    // Compute the real space map correlation coefficient from the summations

  mean_distance  = summed_distances/countd;             // Compute the mean distance travelled in the complex plane 

                            
} 

double ipa_functions::calculate_circular_correlation_coeff(const clipper::HKL_data<clipper::data64::F_phi> &fp1, const clipper::HKL_data<clipper::data64::F_phi>& fp2)
{
  /*
  This function calculates the circular-circular association: T-Linear Association between two phase sets. As per Statistical Analysis of Circular Data, N. I. Fisher pg. 151.
  The resulting value is between -1 and 1, and measures the positive (or negative) association between two sets of angular data. Note there is no scalar association - angles only.
  The hypothesis of no T-Linear Association is rejected if the resultant differs too much from 0.

  We use the computationally efficient form of equation (6.35), shown as equation (6.37)

  */
  double circ_circ_association = 0.0; // Default is 0.0 precisely.

  clipper::HKL_data_base::HKL_reference_index ih;

  double fp1_mean = 0.0;
  double fp2_mean = 0.0;

  double A = 0.0;
  double B = 0.0;
  double C = 0.0;
  double D = 0.0;
  double E = 0.0;
  double F = 0.0;
  double G = 0.0;
  double H = 0.0;

  int n = 0; // How many reflections have been compared.

// Calculate all sum of operations: A, B, C, D, E, F, G, H, Equation 6.37, page 151.
  for (ih = fp1.hkl_info().first(); !ih.last(); ih.next())
  {
    if ( !fp1[ih].missing() && !fp2[ih.hkl()].missing())
    {
      A += cos(fp1[ih].phi() * cos(fp2[ih.hkl()].phi()));
      B += sin(fp1[ih].phi() * sin(fp2[ih.hkl()].phi()));
      C += cos(fp1[ih].phi() * sin(fp2[ih.hkl()].phi()));
      D += sin(fp1[ih].phi() * cos(fp2[ih.hkl()].phi()));
      E += cos(2.0*fp1[ih].phi());
      F += sin(2.0*fp1[ih].phi());
      G += cos(2.0*fp2[ih.hkl()].phi());
      H += sin(2.0*fp2[ih.hkl()].phi());

      n++; // Increment count.
    }
  }

  // Final Association calculation: Equation 6.36, page 151.
  circ_circ_association = (4*((A*B) - (C*D))) / sqrt(((n*n)-(E*E)-(F*F))*((n*n)-(G*G)-(H*H)));

  return circ_circ_association;
}



void ipa_functions::skeletonize_map(clipper::Xmap<float>& xn_input, const double& solvent_density, const double& density_cutoff, const int& sqrrad, const double& beta)
{
  /*
  This function takes a map, and some parameters, and imposes a skeletonisation constraint on the map. 
  cf.  "Baker et al (1993) PRISM: topologically constrained phased refinement for macromolecular crystallography"
  
  The idea is simple, generate the skeleton of the map, then regenerate a new map from the skeleton where the 
  density is replaced by some smoothsly decreasing function of the distance from the skeleton. 
  The distance is measured in Angstroms.
  The skeletonized density is calculated by the function: exp ( -(PI^2 / beta) * distance^2 )
  This function is obviously modulated by the value of beta (Angstroms^-2)
  
  Baker also describe a pruning and thinning process, which removes erroneous skeleton graphs, and minimises them based off their associated density values.

  It is likely this constraint will function best at intermediate resolution, where secondary structures are readily resolvable 
  The intended use of this function is as a real space constraint on the protein region of the map
  
  */

  clipper::Xmap_base::Map_reference_index ix;
  // Selection Sort algorithm to sort the density of the map in ascending order:

  clipper::Xmap<float> new_map;    // New map
  clipper::Xmap<int> skeleton_map; // 1 is skeleton node, 0 is non-skeleton node.
  new_map.init(xn_input.spacegroup(),xn_input.cell(), xn_input.grid_sampling());
  skeleton_map.init(xn_input.spacegroup(),xn_input.cell(), xn_input.grid_sampling());

  int count_skel = 0;
  for (ix = skeleton_map.first(); !ix.last(); ix.next()) {skeleton_map[ix] = 1; new_map[ix] = solvent_density; count_skel++;} // Initialise new_map and skeleton;
  std::cout << "Total map nodes: " << count_skel << std::endl;

  // Create a skeleton from the input map.  
  clipper::Skeleton_basic myskeletonobj;
  myskeletonobj(skeleton_map, xn_input); 

  //for (ix = xn_input.first(); !ix.last(); ix.next()) {xn_input[ix] = 0;} // Reinitialise the map to be solvent.

  count_skel = 0;

  clipper::Cell_descr rcd( xn_input.cell().descr() );
  clipper::Cell_descr vcd( 1.0,1.0,1.0, rcd.alpha(), rcd.beta(), rcd.gamma() );
  clipper::Cell vcell( vcd );

  clipper::Coord_grid g0(-sqrrad,-sqrrad,-sqrrad);
  clipper::Coord_grid g1( sqrrad, sqrrad, sqrrad);
  clipper::Grid_sampling vgrid( 1, 1, 1 );

  clipper::Xmap_base::Map_reference_coord ic(xn_input); // Set to the neighbouring point of interest.
  for (ix = xn_input.first(); !ix.last(); ix.next())
  {
    // For every point that is a skeletal node, and has a large enough density, iterate over its neighbours and set their values accordingly.
    if (skeleton_map[ix] == 1 && xn_input[ix] > density_cutoff) // Density cutoff determines if the skeleton node is large enough to be thinned or not.
    {
      count_skel++; // Count the number of skeleton nodes whilst we're here.
      clipper::Coord_grid iu, iv, iw;
      float thisd2;

      xn_input[ix] = 1.0;
      /*
      for ( iu = g0; iu.u() <= g1.u(); iu.u()++ ) {
        for ( iv = iu; iv.v() <= g1.v(); iv.v()++ ) {
          for ( iw = iv; iw.w() <= g1.w(); iw.w()++ ) {
            thisd2 = iw.coord_frac( vgrid ).lengthsq( vcell ); // Relative distance-squared (Angstroms^2).
            double pi2 = pow(clipper::Util::pi(),2); // π squared.
            double rho = exp( -(pi2 / beta ) * thisd2 ); // Electron density at distance from point.
            rho = std::max(rho, solvent_density); // Set the density to be either the solvent value, or our functioned value, whatever is larger.

            ic.set_coord(ix.coord() + iw); // This is the expensive part.

            xn_input[ic] = std::max((double)xn_input[ic], rho);
            //if (xn_input[ic] < rho)
            //{
            //  xn_input[ic] = rho;
            //}
          }
        }
      } 
      */
    } // Is a skeleton point
  }

  std::cout << "Total skeleton nodes: " << count_skel << std::endl;

}

void ipa_functions::atomize_map_region(clipper::Xmap<float>& map, const clipper::Xmap<int>& mask, const clipper::Atom_list& list_of_atoms, const clipper::HKL_info& hkl, const double threshold, const int region_flag, const double nearest_distance)
{
  /*
  INPUT SUMMARY:
  map               The real-valued map to be modified 
  mask              An integer-valued mask, with the same dimensions as the map, which defines the region to be modified. 
  list_of_atoms     List of clipper atom objects, input atomic coordinates are ignored.
  hkl               Reference to an hkl_info object, used when calculating the Fourier coefficients from the atom list
  threshold         Threshold electron density for atom placement
  region_flag       integer-valued flag, In association with the mask this defines the region considered for atom/pesudo-atom placement.
  nearest_distance  The nearest possible distance two atoms can be to each other. Map positions that are within this distance of placed atoms are ignored.

  This function will apply atomicity or globicity constraints. 
  
  We have modified atomsf.ccp of the clipper library to include the "globic" scattering factors from:

  Acta Cryst. (1995). A51,945-947. Use of globic scattering factors for protein structures at low resolution. By D. Y. GtJo, G. DAVID SMITH, JANE F. GRIFFIN and DAVID A. LANGS

  Globic scatterers represent the spherically-avaraged scattering of the backbone and the side chains
  
  Globic scatterers are identified as pseudo-atoms with the strings:
  "Xbb" (X = represents all amino acids, bb is code for backbone) 
  "Zsc" (Z = single letter representation for a specific amino acid, sc is code for side chain e.g. Vsc = the Valine side chain)

  Thus the protein backbone can be approximated by a number of Xbb pseudo-atoms. Zsc pseudo-atoms could also be included to represent the sidechains. 
  It might be possible to somehow treat these lists of atoms in pairs, so that each bb is associated with its sc in real space.

  As well as these globic scatters, regular atoms can also be used e.g. "C", "N", "O"... etc, see clipper/core/atomsf.ccp for the exhaustive list.
  For simplicity, variable names using "atom" are implied to include the globic pseudo-atoms as well
  
  TODO - Figure out how to properly augment atomsf.ccp without hacking the clipper library.  
  
  To impose an atomicity/globule constraint, we do the following:

  SETUP: Sort the Xmap density Largest -> smallest.
  
  1. Find the largest density within the affected region.
  2. Check it is not within some threshold distance to any previously placed points.
  3. Place an atom or pseudo-atom at this position 
  4. repeat 1-3 until we have iterated through all density points of sufficient density, or have exhausted the number of atoms/pseudo-atoms to place.
  5. Take the atom_list that has been generated and use Clipper::SFcalc_aniso_fft to generate a new Xmap.
  6. Return new Xmap.

  Note, the list_of_atoms is not fundementally changed. But the existing coordinates are ignored.

  TODO: is it worth sorting the list of atoms, so that they are largest to smallest? Who knows.

  */
  clipper::HKL_data<clipper::data64::F_phi> data_sf(hkl); // The Fourier coefficients to be calculated from the atom list.
  clipper::Atom_list placed_atoms; // Container to hold those atoms that were placed.

  int total_atoms = list_of_atoms.size(); // Pre-retrieved total number of atoms being passed into the function.
  int placed_atom_counter = 0; // Iterates through the atom list, anything before here has been placed, if it == total atoms, time to stop.

  // Order the density
  // We only sort those points which are to be considered, i.e. above threshold and within affected region
  std::vector<int> index;
  clipper::Xmap<float>::Map_reference_index ix;
  for (ix = map.first(); !ix.last(); ix.next())
    if (map[ix] > threshold && mask[ix] == region_flag)
      index.push_back(ix.index());

  clipper::Map_index_sort::sort_decreasing(map, index);


  // Iterate over the ordered indices of electron density.

  // This is how skeletonisation does it:
  //for ( int i = 0; i < index.size(); i++ )
  //    if ( !isInSkel( xskl, xskl.coord_of(index[i]), neigh, box_ ) )
  //      xskl.set_data( index[ i ], 0 );

  clipper::Xmap_base::Map_reference_coord ic(map); // Set to the neighbouring point of interest.

  for (int i = 0 ; i < index.size(); i ++ ) 
  {
    // Check against placed atoms.
    bool can_be_placed = true; // Check which remains true if an atom can be placed at this location.

    // Get our position in Orth coordinates (to compare to the atoms)
    clipper::Coord_grid new_point = map.coord_of(index[i]);
    clipper::Coord_orth new_orth = new_point.coord_frac(map.grid_sampling()).coord_orth(map.cell());

    // Check all pre-existing places atoms.
    for (int a = 0; a < placed_atoms.size(); a ++ )
    {
        clipper::Coord_orth existing_orth = placed_atoms[a].coord_orth();
        if (clipper::Coord_orth::length(new_orth, existing_orth) < nearest_distance)
        {
          // We are simply too close, ignore this map point!
          can_be_placed = false;
          break; // No need to check other already placed atoms.
        }
    }

    // Add an atom at this location!
    if (can_be_placed)
    {
      clipper::Atom new_atom(list_of_atoms[placed_atom_counter]);
      placed_atom_counter++;
      new_atom.set_coord_orth(new_orth); // Update the atoms position.
      placed_atoms.push_back(new_atom); // Add our new atom to the list of atoms to be fft'ed.
      //std::cout << "Placed atom at " << new_atom.coord_orth().format() << " with Grid Density " << map.get_data(index[i]) << std::endl;
    }

    if (placed_atom_counter == total_atoms) break; // We're done here.
  }

  clipper::SFcalc_aniso_fft<clipper::ftype64> calc(5.0, 1.5, 0.0); // Radius | Shannon Rate | u-add
  calc(data_sf, placed_atoms); // Calculate F_phi data for new map.
  map.fft_from(data_sf); // Calculate and update the provided map object with the new map SF's.
  return;
}

