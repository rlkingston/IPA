// IPA: Iterative Projection Algorithms for protein crystallography
// Phase consensus worker class  

#include "consensus_worker.h"

consensus_worker::consensus_worker()
{
  // Null constructor.
  is_idle = true;
}

consensus_worker::consensus_worker(consensus_phase_settings new_settings)
{
  // Initiate object with settings already in place.
  settings = new_settings;
  is_idle = true;
}

consensus_worker::~consensus_worker()
{
  // Deconstructor, simply closes the log file.
  outlog.close();
}

bool consensus_worker::add_phase_input_to_list(std::string solution_filename)
{
  /*
   This public function is the only safe way to add files to the phase agreement analysis. Upon completion of each phase determination run, the output phase file gets added to the list for analysis
   
   Importantly, the ordering of addition here needs to remain consistent throughout the experiment, that is the push_back of the flag into files_processed_flag and the filename into phase_file_list needs to be concurrent
  
   Due to this need, the function only adds unique files to the end of the list, and the list cannot be accessed otherwise ... once something is added, its part of the experiment forever. 
  
   If you want to change your list, then you make a new consensus_worker object.
  
   Because adding files to the list is controlled, the entire phase_list.txt file from the phase_determination step can simply be added, without worrying about duplication, as an internal check is made.
   
   // TODO: Note that the run_id is not carried over from the filenames in the phase_list.txt file, and the filenames are not strictly ordered by run_id since different threads finish at different times. 
   // Hence the file_id used during clustering, reflects the order of file addition, and is not the same as the run_id of the inputs ... 
   // We could attempt to fix this with a lookup table ... At the moment, we simply report the relevant filenames (which contain the run_id) at the conclusion of clustering.

   */
  
  
  if (is_idle) // We cannot add new phase sets to the system, if determination of a consensus phase set is in progress
  {
  
  for (int i = 0 ; i < phase_file_list.size(); i ++ )
  {
    if (phase_file_list[i] == solution_filename)
    {
      // We're trying to add something already in the phase_file_list, so ignore this entry
      //std::cout << "The filename was detected to be identical" << std::endl;
      return true;
    }
  }
  
  phase_file_list.push_back(solution_filename);
  
  files_processed_flag.push_back(false); // By default, files have NOT been processed.
  
  n_total_files = phase_file_list.size();
  
  return true;
  }
  return false;
}

std::ostream& consensus_worker::return_list_of_consensus_solutions(std::vector<cluster_data>& output_data, std::ostream& outlog)
{
  /*
   This function returns all the saved consensus phase sets from the experiment.
   It will output the list of files, and associated information into the requested log
   */
  
  //reportlog << "Retrieving consensus phase sets: " << std::endl;
  
  //Lets start by ensuring that the output is valid, i.e.  all the files HAVE been processed:
  for ( int i = 0 ; i < files_processed_flag.size(); i ++ )
  {
    if (!files_processed_flag[i])
    {
      outlog << util::parse_line(78," Not all the files added to the consensus worker have been processed!");
      outlog << util::parse_line(78," The following outputs do not represent the whole dataset:");
      break; // Was going to return false, but you never know, probably better to continue to give something.
    }
  }
  
  int n_solutions = all_prototypes.size();
  if (n_solutions == 0 )
  {
    // Find the most consistent pair of phase sets, to give the user something.
    double minimum_phase_difference = 360.0;
    int file1 = -1;
    int file2 = -1;
    
    int array_loc = -1; // This is actually unused, but it stores the array location of the minimum entry of the agreement matrix.
    // Simply find the minimum phase difference, and return the associated file id's 
    for (int i = 0; i < agreement_matrix.size(); i++ )
    {
      if (agreement_matrix[i] < minimum_phase_difference)
      {
        minimum_phase_difference = agreement_matrix[i];
        matrix_1d_to_filexy(i,file1,file2); // Convert coordinate in the 1D array to x and y coordinates, which are the required file ids.
        array_loc = i;
      }
    }
    
    util::print_section_line(78, outlog);
    outlog << util::parse_line(78, " NO HIGHLY CONSISTENT PHASE SETS WERE FOUND ");
    util::print_empty_line(78, outlog);
    outlog << util::parse_line(78, "The closest sets had a mean absolute phase difference of " + std::to_string(minimum_phase_difference));
    outlog << util::parse_line(78, "A pair of statistically independent phase sets will have a phase");
    outlog << util::parse_line(78, " difference of ~90°.", 1); // Split and fixed for special characters.
    outlog << util::parse_line(78, "A pair of true solutions will usually have a phase difference < 50°.", 1); // Fixed for special characters
    outlog << util::parse_line(78, "An intermediate value may indicate some progression toward location of the solution");
    outlog << util::parse_line(78, "The two closest phase sets are:");
    outlog << util::parse_line(78, phase_file_list[file1]); // return the filename associated with this file_id
    outlog << util::parse_line(78, phase_file_list[file2]); // return the filename associated with this file_id
    util::print_section_line(78, outlog);

    // TODO: Include reporting of the map correlation coefficients here, since we do now calculate them, and might be a more digestible metric than mean aboslute phase difference
    
     return outlog;
  }
  
  output_data.clear();
  
  // Copy all cluster prototypes into the output reference.

  // This line is somewhat redundant,
  //outlog << util::parse_line(78, "Retrieving " + std::to_string(all_prototypes.size()) + " consensus solutions from consensus worker.");
  
  for (int i = 0 ; i < all_prototypes.size(); i ++ )
  {
    output_data.push_back(all_prototypes[i]);
  }
  
  return outlog;
}

bool consensus_worker::generate_consensus_phases()
{
  /*
   We call this function after running clustering algorithm dbscan . 
   THis function generates the mean (consensus) phase set for each cluster, aka the cluster prototypes, and writes the outputs to the relevant folders.
   The function also ensures that the cluster_list is updated. The cluster_list can be accessed from other classes via the return_list_of_consensus_solutions method.
   */
  
  // TODO: This functon is too big, outputs too much, and isn't tidy. Clean it up sometime - especially the closing parts - which could be folded into final reporting
  
  // Clear and resize the cluster average array so it can contain an average for each cluster.
  
  cluster_averages.clear();
  cluster_averages.resize(number_of_clusters, input_fp[0]);
  
  clipper::HKL_data<clipper::data64::F_phi> work_fp(hkl); // We need this later for writing mtzs.
  
  outlog << "*------------------------------------------------------*" << std::endl;
  outlog << "Generating mean phases for each cluster" << std::endl;
  outlog << "*------------------------------------------------------*" << std::endl;
  
  for(int ic = 0; ic < number_of_clusters; ic++) // loop over the clusters
  {
    cluster_data current_prototype;
    
    // create an object (type  Phi + fom)  to hold both mean phase and sample circular variance
    
    clipper::HKL_data<clipper::data64::Phi_fom> mean_phase_and_scv(hkl);
    
    int im_ref;
    
    outlog << "\nCluster: " << ic << std::endl;
    outlog << "Members of this cluster: [";

    int unique_member_id = 0;
    int n_members_identified = 0;
    
    std::vector<int> list_of_members;
    std::vector<float> list_of_agreements;
    
    for(size_t im = 0; im < input_fp.size(); im++) // loop over all the input sets of structure factors
    {
      if (cluster_ids[im] == ic)
      {
        current_prototype.list_of_members.push_back(phase_file_list[im]); // Add the name of the file to the list of members
        n_members_identified += 1;
        unique_member_id += (im + im * 10); // just add together all the ID's in a way that'll be unique.
        if (n_members_identified == 1)  im_ref = im; // the first member of the cluster is the "reference", which is used to define the origin
        
        list_of_members.push_back(im);
        outlog << im << ", ";
      }
    }
    outlog << "]." << std::endl;

    
    // Get the mean mean agreement (collects all the pairwise mean agreements, and averages them). 
    float mean_agreement = 0;
    int mean_count = 0;
    for (int y = 0 ; y < list_of_members.size(); y ++ )
    {
      for (int x = y + 1 ; x < list_of_members.size(); x ++ )
      {
        if (x != y)
        {
          // x and y are a pair of cluster id's, they need to be averaged.
          int id = file_xy_get_1d(list_of_members[x],list_of_members[y]); // Be careful here, as j must be larger than i, else we break the diagonal matrix.
          outlog << "Whats connected: x: " << list_of_members[x] << " and y: " << list_of_members[y] << std::endl;
          outlog << "Furthermore, the mean absolute phase difference is: " << agreement_matrix[id] << std::endl;
          mean_agreement += agreement_matrix[id];
          mean_count++; // increment total
        }
      }
    }
    current_prototype.mean_phase_difference_across_cluster = mean_agreement / mean_count;
    
    outlog << std::endl;
    outlog << "Number of members of this cluster: " << n_members_identified << std::endl;
    outlog << "Origin is defined by member : " << im_ref << std::endl;
    outlog << "Average of the pairwise mean absolute phase differences for this cluster: " << current_prototype.mean_phase_difference_across_cluster << std::endl;
    outlog << std::endl;
    
    current_prototype.epsilon = settings.epsilon_input;
    current_prototype.cluster_id = ic;
    current_prototype.n_members = n_members_identified;
    
    // loop over the asymmetric unit and average the origin-corrected phases for this cluster
    // These  are variables on the unit circle ...
    // see Statistical Analysis of Circular data. Fisher, Section 2.3.1
    
    for ( clipper::HKL_data_base::HKL_reference_index ih = hkl.first(); !ih.last(); ih.next() )
    {
      double sum_sin = 0.0;
      double sum_cos = 0.0;
      double sum_amp = 0.0;
      
      for(int im = 0; im < input_fp.size(); im++) // loop over all the input sets of structure factors
      {
        if (cluster_ids[im] == ic) // if true, the phase set is a member of the current cluster
        {
          if (!input_fp[im][ih].missing()) // Check data isn't missing
          {
            double delta_phi = 0.0; // By default, no origin shift necessary.
            bool inverse_required = false;
            
            if (settings.allow_origin_shifts)
            {
              // Get the lookup x and y for the 1D origin shifts array between im and im_ref
              int x = std::max(im_ref,im);
              int y = std::min(im_ref,im);
              
              if (x != y) // Ensure I'm not the im_ref.
              {
                int matrix_lookup = file_xy_get_1d(x,y); // Get the postion in the 1D array from the x and y coordinates (file_ids)

                if (im > im_ref) // We don't neet to invert or shift the reference !!
                {
                  inverse_required = inverse_flags[matrix_lookup]; // Do we need to invert the phase set im to maximize agreement?. The matrix inverse_flags tells us if a phase set is inverted, relative to im_ref
                                    
                  clipper::Coord_frac origin_shift = shift_coords[matrix_lookup]; // here is the real space origin shift required to maximize the agreement
                  
                  clipper::Coord_reci_frac h( ih.hkl() ); // fractional coordinates of current observation ih
                  const double hx = h * origin_shift;
                  delta_phi = clipper::Util::twopi() * ( hx - floor(hx) );    // Convert the real space origin shift into a phase difference
                }
              }
            }
            
            //
            if (inverse_required)  // We need to invert the data to match the reference
            {
              // We could make the decision to permanently invert the inputs, and record if this occurred.
              // Then we wouldnt have to do this here.
              input_fp[im][ih].friedel();
              sum_sin += sin( input_fp[im][ih].phi() + delta_phi );
              sum_cos += cos( input_fp[im][ih].phi() + delta_phi );
              input_fp[im][ih].friedel();
            }
            else // data has the same hand as the reference
            {
              sum_sin += sin( input_fp[im][ih].phi() + delta_phi );
              sum_cos += cos( input_fp[im][ih].phi() + delta_phi );
            }
            sum_amp += input_fp[im][ih].f();
          }
        }
      }
      
      // compute the mean phase
      
      cluster_averages[ic][ih].phi() = atan2(sum_sin,sum_cos);
      
      // atan2 returns the result over the interval -pi to pi. Correct so the phase is specified over the interval 0 to 2*pi.
      if (cluster_averages[ic][ih].phi() < 0)
      {
        cluster_averages[ic][ih].phi() += clipper::Util::twopi();
      }
      
      // compute the mean structure factor amplitude 
      
      // For the missing data the amplitudes are not subject to any Fourier space constraint, and may differ between the input data sets  ...
      // For all other data, the terms are likely to be identical, and the averaging has no effect. 
      // This is because the ER algorithm will almost always be used at the end of the phase determination step
      // Hence Fourier space projection will be the last operation performed and the amplitudes will have been reset their measured values (multipled by any apodization function)

      // TODO: Compute the standard deviation of the structure factor amplitudes. Might be useful to see how much variation there is in the amplitude estimates for the missing data, between the differing solutions
           
      cluster_averages[ic][ih].f() = sum_amp / double(n_members_identified);
            
      // compute the sample circular variance ... see Fisher page 32
      
      float resultant_length = sqrt( pow(sum_cos, 2.0 ) + pow(sum_sin, 2.0 ) );
      float mean_resultant_length =  resultant_length / float(n_members_identified);
      float sample_circular_variance = 1.0 - mean_resultant_length;
      
      // as of right now,  the sample circular variance is not retained beyond this loop
      // and we don't use the duplicated mean phase estimate for anything
      // however let's keep these two related quantities together here for completeness, allowing for future development ...
      
      mean_phase_and_scv[ih].phi() = cluster_averages[ic][ih].phi();
      mean_phase_and_scv[ih].fom() = sample_circular_variance;
      
    } // end the averaging loop over the asymmetric unit
    
    // Report the sample circular variance for this cluster
    
    clipper::Resolution_ordinal resord; // Resolution Ordinal - so far as I can figure out, once initiated this returns the approximate fractional ordinal within a dataset (in the range 0...1) for any specified value of 1 / s^2
    resord.init( hkl, 1.0 );
    
    // initialize vectors with the appropriate number of elements, all set to zero
    // bins 0,1, ...,n-1  store the individual bin statistics. bin n stores statistics for the entire data set.
    
    std::vector<double> observation_counts_all(settings.res_bins + 1, 0.0);
    std::vector<double> scv_all(settings.res_bins + 1, 0.0);
    std::vector<double> weights_all(settings.res_bins + 1, 0.0);
    
    std::vector<double> observation_counts_centric(settings.res_bins + 1, 0.0);
    std::vector<double> scv_centric(settings.res_bins + 1, 0.0);
    std::vector<double> weights_centric(settings.res_bins + 1, 0.0);
    
    std::vector<double> observation_counts_acentric(settings.res_bins + 1, 0.0);
    std::vector<double> scv_acentric(settings.res_bins + 1, 0.0);
    std::vector<double> weights_acentric(settings.res_bins + 1, 0.0);
    
    for (clipper::HKL_data_base::HKL_reference_index ih = hkl.first(); !ih.last(); ih.next() )
    {
      if ( !mean_phase_and_scv[ih].missing() )
      {
        const int bin = std::min ( int(double(settings.res_bins) * resord.ordinal(ih.invresolsq())) , settings.res_bins-1 ); // Which bin does this observation belong to ?
        const double w = 1.0 / ih.hkl_class().epsilon();  // Get the standard statistical weights
        const double wa = zero_phase_weight[ih].fom(); // Get the weights associated with the apodization function
        
        //  outlog << ih.hkl().h() << " " << ih.hkl().k() << " " << ih.hkl().l() << " " << std::endl;
        //  outlog <<  mean_phase_and_scv[ih].phi() << " "  <<  mean_phase_and_scv[ih].fom() << std::endl;
        
        observation_counts_all[bin]      += 1.0;
        scv_all[bin]                     += w * wa * mean_phase_and_scv[ih].fom();
        weights_all[bin]                 += w * wa;
        
        observation_counts_all[settings.res_bins] += 1.0;
        scv_all[settings.res_bins]                += w * wa * mean_phase_and_scv[ih].fom();
        weights_all[settings.res_bins]            += w * wa;
        
        if ( ih.hkl_class().centric() )
        {
          observation_counts_centric[bin]       += 1.0;
          scv_centric[bin]                      += w * wa * mean_phase_and_scv[ih].fom();
          weights_centric[bin]                  += w * wa;
          
          observation_counts_centric[settings.res_bins]  += 1.0;
          scv_centric[settings.res_bins]                 += w * wa * mean_phase_and_scv[ih].fom();
          weights_centric[settings.res_bins]             += w * wa;
        }
        else
        {
          observation_counts_acentric[bin]      += 1.0;
          scv_acentric[bin]                     += w * wa * mean_phase_and_scv[ih].fom();
          weights_acentric[bin]                 += w * wa;
          
          observation_counts_acentric[settings.res_bins] += 1.0;
          scv_acentric[settings.res_bins]                += w * wa * mean_phase_and_scv[ih].fom();
          weights_acentric[settings.res_bins]            += w * wa;
        }
      }
      
    }
    
    
    // Compute the mean scv
    
    for ( int bin = 0; bin < (settings.res_bins + 1); bin++ )
    {
      scv_all[bin]      /= std::max( weights_all[bin],      1.0 );
      scv_centric[bin]  /= std::max( weights_centric[bin],  1.0 );
      scv_acentric[bin] /= std::max( weights_acentric[bin], 1.0 );
    }
    
    // statistics as a function of resolution
    
    resord.invert(); // invert the distribution so that now we specify the ordinal value and it returns 1/|s|**2
    
    double lower_limit;
    double upper_limit;
    
    outlog << "\nWeighted sample circular variance as a function of resolution" << std::endl;
    outlog <<   "------------------------------------------------------------" << std::endl;
    
    outlog << "                                        all data              centric data             acentric data          " << std::endl;
    outlog << " |s| lo   |s| hi  1/|s| hi       #   mean w  mean scv      #   mean w  mean scv      #   mean w  mean scv     " << std::endl;
    
    for ( int bin = 0; bin < settings.res_bins; bin++ )
    {
      if (bin == 0)
      {
        lower_limit = 0.0;
        upper_limit = sqrt( resord.ordinal( double(bin+1)/double(settings.res_bins) ) );
      }
      else
      {
        lower_limit = sqrt( resord.ordinal( double(bin)/double(settings.res_bins) ) );
        upper_limit = sqrt( resord.ordinal( double(bin+1)/double(settings.res_bins) ) );
      }
      outlog << std::fixed <<
      std::setw(8) << lower_limit << " " <<
      std::setw(8) << upper_limit << " " <<
      std::setw(8) << 1/upper_limit <<  " "  <<
      std::setw(8) << int(observation_counts_all[bin]) <<  " "  <<
      std::setw(8) << weights_all[bin]/observation_counts_all[bin] <<  " "  <<
      std::setw(8) << scv_all[bin] <<
      std::setw(8) << int(observation_counts_centric[bin]) <<  " "  <<
      std::setw(8) << weights_centric[bin]/observation_counts_centric[bin] <<  " "  <<
      std::setw(8) << scv_centric[bin] <<
      std::setw(8) << int(observation_counts_acentric[bin]) <<  " "  <<
      std::setw(8) << weights_acentric[bin]/observation_counts_acentric[bin] <<  " "  <<
      std::setw(8) << scv_acentric[bin] << std::endl;
    }
    
    
    // overall statistics
    
    lower_limit = 0.0;
    upper_limit = sqrt(resord.ordinal( 1.0 ));
    outlog << "Overall statistics:" << std::endl;
    outlog << std::fixed <<
    std::setw(8) << lower_limit << " " <<
    std::setw(8) << upper_limit << " " <<
    std::setw(8) << 1/upper_limit <<  " "  <<
    std::setw(8) << int(observation_counts_all[settings.res_bins]) <<  " "  <<
    std::setw(8) << weights_all[settings.res_bins]/observation_counts_all[settings.res_bins] <<  " "  <<
    std::setw(8) << scv_all[settings.res_bins] <<
    std::setw(8) << int(observation_counts_centric[settings.res_bins]) <<  " "  <<
    std::setw(8) << weights_centric[settings.res_bins]/observation_counts_centric[settings.res_bins] <<  " "  <<
    std::setw(8) << scv_centric[settings.res_bins] <<
    std::setw(8) << int(observation_counts_acentric[settings.res_bins]) <<  " "  <<
    std::setw(8) << weights_acentric[settings.res_bins]/observation_counts_acentric[settings.res_bins] <<  " "  <<
    std::setw(8) << scv_acentric[settings.res_bins] << std::endl;
    
    current_prototype.scv = scv_all[settings.res_bins];
        
    // TODO: here we are overwriting the consensus phase estimates everytime the function is called
    // It would require a careful analysis of the DBSCAN algorithm, but I'm not certain this is safe. The number of clusters may shrink as more solutions were added to the analysis. 
    // Then we would end up with an "orphaned" consensus phase file that was no longer relevant to the analysis. 
    // Seems like this could only happen due to some unfortunate arrangement of boundary points.
    // Avoiding the need for hard thinking, it might be simpler and cleaner to write to a different subdirectory every time phase clustering is invoked.
    // CURERNT FIX: All consensus phase solutions have a unique ID calculated by the mixture of file ids that went into them, thus no over-writes can occur, the cluster membership can always be tracked down in the log.
    
    char output_filename[256] = "";
    clipper::String output_col_fp         = "/*/*/[F,PHIC]";
    clipper::CCP4MTZfile mtzout;
    
    // Write out the mean phase estimates and correspondent map for this cluster
    //--------------------------------------------------------------------------
    
    //Altered below to ensure it goes to the out directory.
    snprintf(output_filename, sizeof(output_filename),"%sconsensus_phase_%d_%i.mtz", settings.working_dir.c_str(), ic, unique_member_id);
    outlog << "Writing Amplitudes and Phases to file: " << output_filename << std::endl;
    current_prototype.consensus_filename = output_filename; // this reports the mtz filename only, but it should be clear enough. We can tidy up later if needed

    mtzout.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
    mtzout.open_write( output_filename );
    mtzout.export_hkl_info( hkl );
    mtzout.export_hkl_data(cluster_averages[ic], output_col_fp);
    mtzout.close_write();

    // For convenience of the user, we also output a finely gridded map
    snprintf(output_filename, sizeof(output_filename),"%sconsensus_phase_%d_%i.ccp4", settings.working_dir.c_str(), ic, unique_member_id);
    outlog << "Writing correspondent map to file: " << output_filename << std::endl;
    
    clipper::Grid_sampling output_grid;
    clipper::Xmap<float> output_map;
    clipper::CCP4MAPfile map_out;
     
    output_grid.init(hkl.spacegroup(),hkl.cell(), hkl.resolution(), 3.0); // create grid object with Linear Shannon rate of 3.0 
    output_map.init(hkl.spacegroup(), hkl.cell(), output_grid );
    output_map.fft_from( cluster_averages[ic] );  // Fourier transform to generate the map
  
    map_out.open_write(output_filename);
    map_out.export_xmap(output_map);
    map_out.close_write();
  
    // Write out the inverted mean phase estimates and correspondent map for this cluster
    //-----------------------------------------------------------------------------------

    //Altered below to ensure it goes to the out directory.
    snprintf(output_filename, sizeof(output_filename), "%sconsensus_phase_%d_%i_inverse.mtz", settings.working_dir.c_str(),  ic, unique_member_id);
    outlog << "Writing Amplitudes and Phases to file: " << output_filename << std::endl;
    current_prototype.consensus_filename_inverse = output_filename; // this reports the mtz filename only, but it should be clear enough. We can tidy up later if needed
       
    for ( clipper::HKL_data_base::HKL_reference_index ih = cluster_averages[ic].first(); !ih.last(); ih.next() )
    {
      work_fp[ih] = cluster_averages[ic][ih];
      work_fp[ih].friedel(); // invert the phases
    }
     
    if (!space_group_is_chiral) // space group is achiral, and is invariant under inversion
    {

      mtzout.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
      mtzout.open_write( output_filename );
      mtzout.export_hkl_info( hkl );
      mtzout.export_hkl_data(work_fp, output_col_fp);
      mtzout.close_write();

      // For convenience of the user, we also output a finely gridded map
      snprintf(output_filename, sizeof(output_filename),"%sconsensus_phase_%d_%i_inverse.ccp4", settings.working_dir.c_str(), ic, unique_member_id);
      outlog << "Writing correspondent map to file: " << output_filename << std::endl;
        
      output_grid.init(hkl.spacegroup(),hkl.cell(), hkl.resolution(), 3.0); // create grid object with Linear Shannon rate of 3.0 
      output_map.init(hkl.spacegroup(), hkl.cell(), output_grid );
      output_map.fft_from( work_fp );  // Fourier transform to generate the map
    
      map_out.open_write(output_filename);
      map_out.export_xmap(output_map);
      map_out.close_write();
      
    }
    else  // space group is chiral, and changes under inversion
    {  
      // There may well be a simpler and neater way to change the hand of the space group, but the following works.
      
      // create a new HKL_info ... the same as the original but switch to the enantiomorphic space group
      int sg_enantiomorph = util::enantiomorphic_space_group(hkl.spacegroup().spacegroup_number());
      clipper::HKL_info hkl_enantiomorph( clipper::Spacegroup( clipper::Spgr_descr( sg_enantiomorph ) ), hkl.cell(), hkl.resolution() ); 
      hkl_enantiomorph.generate_hkl_list(); // initially there are no observations hkl in the list - better make some
            
      clipper::HKL_data<clipper::data64::F_phi> work_fp_enantiomorph(hkl_enantiomorph);

      for ( clipper::HKL_data_base::HKL_reference_index ih = hkl_enantiomorph.first(); !ih.last(); ih.next() ) {
        work_fp_enantiomorph[ih] = work_fp[ih.hkl()];
      }

      mtzout.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
      mtzout.open_write( output_filename );
      mtzout.export_hkl_info(hkl_enantiomorph );
      mtzout.export_hkl_data(work_fp_enantiomorph, output_col_fp);
      mtzout.close_write();

      // For convenience of the user, we also output a finely gridded map
      snprintf(output_filename, sizeof(output_filename),"%sconsensus_phase_%d_%i_inverse.ccp4", settings.working_dir.c_str(), ic, unique_member_id);
      outlog << "Writing correspondent map to file: " << output_filename << std::endl;
        
      output_grid.init(hkl_enantiomorph.spacegroup(),hkl_enantiomorph.cell(), hkl_enantiomorph.resolution(), 3.0); // create grid object with Linear Shannon rate of 3.0 
      output_map.init(hkl_enantiomorph.spacegroup(), hkl_enantiomorph.cell(), output_grid );
      output_map.fft_from( work_fp_enantiomorph );  // Fourier transform to generate the map
    
      map_out.open_write(output_filename);
      map_out.export_xmap(output_map);
      map_out.close_write();
      
    }

      
    // If a known phase set is input, then compare with the cluster prototypes ...
    
    // TODO: ... there is another issue to be resolved here ...
    // If the amplitude data in the input phase sets has been apodized, then strictly speaking the known amplitude data in fp_known would be equivalently apodized
    // (This currently does not happen). The place this is important is in the call to the phased translation function, from within
    // compute_phase_agreement, which determines the origin shift between the two phases sets. Probably we will still get the same origin shift if one set of data
    // is unapodized, and the other is not. However ideally, we should ensure both data sets are apodized equivalently. 
    // To make this work robustly, we need to know what apodization function had been applied to the input data sets without reference to the settings etc 
    // The best solution long term solution would be to specify the apodization function via some metadata added to the map or mtz header
    // Then the apodization scheme would be specified together with the data, eliminating many possibilities for error.
       
    current_prototype.mean_phase_diff_all_with_known = 0.0;
    current_prototype.mean_phase_diff_all_with_known_inverse = 0.0;
    if (settings.input_known_phases)
    {
      outlog << "\nAgreement with Known Phase Set, including inversion, if permitted by the space group" << std::endl;
      outlog << "-------------------------------------------------------------------------------------\n" << std::endl;
      
      
      float phase_diff_all = 0.0;
      float phase_diff_centric = 0.0;
      float phase_diff_acentric = 0.0;
      float map_cc = 0.0;
            
      bool allow_origin_shifts = true; // This the intended behavior. In this setting we always allow origin shifts.
      clipper::Coord_frac origin_shift( 0.0, 0.0, 0.0 );
      
      bool invert;
      
      invert = false;
           
      outlog << "\n Agreement between known phases and the mean phases " << "\n" << std::endl;
      outlog << std::setw(10);
      
      ipa_functions::compute_phase_agreement(fp_known, cluster_averages[ic], zero_phase_weight,
                                             settings.res_bins, allow_origin_shifts, invert, phase_diff_all, phase_diff_centric,  phase_diff_acentric, map_cc, origin_shift, outlog, true);
      
      current_prototype.mean_phase_diff_all_with_known = phase_diff_all;
      
      // If requested, apply the origin shift to maximize agreement with the known phase set, and write out the modified phases
      
      if (settings.shift_to_known_origin)
      {
        clipper::HKL_data<clipper::data64::F_phi> fp_shifted(hkl);
        ipa_functions::shift_phases(cluster_averages[ic], fp_shifted, origin_shift);
        
        char output_filename[128] = "";
        // TODO - change the output filename here so that it matches naming above ...
        // Also write out a map for convenience, as above
        snprintf(output_filename, sizeof(output_filename), "min_points_%d_epsilon_%d_prototype%d_shifted.mtz",     settings.min_pts, static_cast<int>(100*settings.epsilon_input), ic);
        
        mtzout.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
        mtzout.open_write( output_filename );
        mtzout.export_hkl_info( hkl );
        mtzout.export_hkl_data(fp_shifted, output_col_fp);
        mtzout.close_write();
      }
      
      
      if (!space_group_is_chiral) // space group is invariant under change of hand
      {
        invert = true;
        //outlog << "\nAgreement between known phases and the *inverted* mean phases for cluster: " << i << "\n" << std::endl;
        outlog << "\nAgreement between known phases and the *inverted* mean phases" << "\n" << std::endl;
        
        ipa_functions::compute_phase_agreement(fp_known, cluster_averages[ic], zero_phase_weight,
                                               settings.res_bins, allow_origin_shifts, invert, phase_diff_all, phase_diff_centric,  phase_diff_acentric, map_cc, origin_shift, outlog, true);
        
        // If requested, apply the origin shifts to match with the known phase set, and write out the modified phases
        current_prototype.mean_phase_diff_all_with_known_inverse = phase_diff_all;
        
        if (settings.shift_to_known_origin)
        {
          for (clipper::HKL_data_base::HKL_reference_index ih = cluster_averages[ic].first(); !ih.last(); ih.next() )
          {
            work_fp[ih] = cluster_averages[ic][ih];
            work_fp[ih].friedel(); // invert the phases
          }
          
          clipper::HKL_data<clipper::data64::F_phi> fp_shifted(hkl);
          ipa_functions::shift_phases(work_fp, fp_shifted, origin_shift);
          
          char output_filename[128] = "";
          // TODO - change the output filename here so that it matches naming above ...
          // Also write out a map for convenience, as above
          snprintf(output_filename, sizeof(output_filename), "min_points_%d_epsilon_%d_prototype%d_inverted_shifted.mtz",     settings.min_pts, static_cast<int>(100*settings.epsilon_input), ic);
          
          mtzout.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
          mtzout.open_write( output_filename );
          mtzout.export_hkl_info( hkl );
          mtzout.export_hkl_data(fp_shifted, output_col_fp);
          mtzout.close_write();
          
        }
      } // End of if Chiral
      
      
    } // End of if InputPhaseknown
    
    
    // Store all the information about this prototype
    
    all_prototypes.push_back(current_prototype);
    
  }
  
  return true;

}

bool consensus_worker::check_for_consensus_phases(std::atomic<float>& MyProgress)
{
  is_idle = false; // Lock down the object for additional inputs whilst it calculates.
  /*
   Here we perform everything that is needed to complete the phase consensus, however, this is done without recomputing anything that has been calculated in previous steps. 
   
   The steps to perform phase consensus are the following:
      
   - Import data from the file list into the vector input_fp.
   
   - Augment the matrices containing distances and origin shifts using the new data
   
   - execute DBSCAN to identify clusters.
   
   - Generate averages for each cluster 
   
   - Output the cluster averages (and their inverse) as Fourier coefficients and maps.  These are the potential solutions to the phase retrieval problem.
   
   Note the MyProgress atomic float needs to be 9.0 - 9.99 to indicate the progress for phase consensus to the threader.
   
   */
  
  // Quick check to see if we even have a phase file to do this...
  if (phase_file_list.size() == 0)
  {
    return false;
  }
  // if we haven't already, initialise our constants, and start a log.
  if (!is_init)
  {
    is_init = initialise();
    if (is_init != true)
    {
      outlog << "Something has gone horribly wrong... or minorly so... check your inputs." << std::endl;
      return false;
    }
  }
  
  // ------------------------------------------------------- //
  // --------- READ FILES, POPULATE MATRIX  ---------------- //
  // ------------------------------------------------------- //
  
  for (int i = 0; i < n_total_files; i++) // For each phase file.
  {
    MyProgress = 9 + 0.80*(float(i) / float(n_total_files));
    if (files_processed_flag[i] == false) // For each file not yet processed.
    {
      // Import file into data input_fp vectors.
      outlog << "*------------------------------------------------------*" << std::endl;
      input_data_from_filename(i);
      
      
      
      if (settings.input_known_phases)
      {
        // TODO: Calculate agreement with known phase set, if desired. May be useful to have the individual agreements reported, as well as agreement with the consensus
      }
      
      outlog << "Calculating agreement matrix for new input_fp[" << i << "]." << std::endl;
      outlog << "*------------------------------------------------------*" << std::endl;

      // Append the agreement data for this file to the 1d matrix
      calculate_and_append_data_to_matrix(i);
      
      files_processed_flag[i] = true; // This file has now been processed.
    }
  }
  outlog << "*------------------------------------------------------*" << std::endl;
  outlog << "Matrix representation for input phase files:" << std::endl;
  print_agreement_matrix_to_log(agreement_matrix,outlog);
  if (space_group_is_chiral)
  {
    outlog << "No inversion was performed as the space group is chiral." << std::endl;
    
  } else // space group is achiral
  {
    print_inverse_matrix_to_log(inverse_flags,outlog);
  }
  outlog << "*------------------------------------------------------*" << std::endl;

  // ------------------------------------------------------- //
  // --------- DBSCAN INTO CLUSTERS & AVERAGE -------------- //
  // ------------------------------------------------------- //
  
  // At this point, all files have have been processed and the matrices containing distances, origin shifts, and inversion indicators have been generated
  outlog << "Performing DBSCAN to identify clusters: " << std::endl;
  outlog << "Threshold mean absolute phase difference for clustering: " << settings.epsilon_input << std::endl;
  // Time to perform DBSCAN to generate clusters.
  
  // Assumes and reads the 1D agreement (distance) matrix.
  number_of_clusters = 0;
  cluster_ids.clear();
  
  dbscan(settings.epsilon_input, settings.min_pts, number_of_clusters, cluster_ids); // execute DBSCAN to generate clusters
  
  MyProgress = 9 + 0.95; // Just an abitrary amount to indicate progress is almost done.
  
  log_cluster_info(); // Log the current cluster information.
  
  outlog << "Averaging across clusters to generate consensus phases" << std::endl;
  generate_consensus_phases();
  
  outlog << "Finished averaging !" << std::endl;
  outlog << "*------------------------------------------------------*" << std::endl;
  
  
  std::vector<cluster_data> temp;
  return_list_of_consensus_solutions(temp, outlog); // Write this to the outlog as well.

  is_idle = true; // Release the atomic bool check, we can safely add more inputs to the consensus worker.
  
  return true;
}

bool consensus_worker::input_data_from_filename(int fileid)
{
  /*
   Looks for filename found at list_of_files[fileid] and adds the dataset to vector input_fp in the same position.
   */
  if (input_fp.size() == fileid) // Check to ensure we're adding to the end, which is also the correct place
  {
    clipper::CCP4MTZfile mtzin;
    outlog << "Reading input_fp[" << fileid << "] from file \"" << phase_file_list[fileid] << "\"" << std::endl;
    mtzin.open_read( phase_file_list[fileid] ); // open the mtz file for reading
    clipper::HKL_data<clipper::data64::F_phi> work_fp(hkl);
    mtzin.import_hkl_data(work_fp, settings.col_fp);
    mtzin.close_read();
    input_fp.push_back(work_fp); // augment the vector
    return true;
  }
  outlog << "ERROR: The input filename [" << fileid << "] is trying to input itself at the wrong location [" << input_fp.size() << "]..." << std::endl;
  
  return false;
}

bool consensus_worker::calculate_and_append_data_to_matrix(int fileid)
{
  /*
   Here we process the addition of input_fp[fileid] to the matrices, it is important that this is done sequentially, and in the order of the file ids... we cannot add a file to a position without over-writing. 
   As we should never "need" to over-write, we simply add to the end, and check at the start whether our starting place is valid.
   
   Note that this function will increase in overhead by 1/2(N*N)-N complexity, where N is fileid.
   i.e. for each entry added, we must calculate the agreements (distances) with all previous entries.
   
   */


  for (int i = 0; i < fileid; i ++) // For every file before this one...
  {
    // Important variables to populate:
    float phase_diff_agreement = 0.0;
    float phase_diff_agreement_with_inversion = 0.0;
    clipper::Coord_frac origin_shift( 0.0, 0.0, 0.0 );
    clipper::Coord_frac origin_shift_with_inversion( 0.0, 0.0, 0.0 );
    
    // Made simply for references to pass:
    float phase_diff_centric = 0.0;
    float phase_diff_acentric = 0.0;
    float map_cc = 0.0;
    

    // First we assess agreement without inverting
    bool invert = false;
    ipa_functions::compute_phase_agreement(input_fp[i],
                                           input_fp[fileid],
                                           zero_phase_weight,
                                           settings.res_bins,
                                           settings.allow_origin_shifts,
                                           invert,
                                           phase_diff_agreement,
                                           phase_diff_centric,
                                           phase_diff_acentric,
                                           map_cc,
                                           origin_shift,
                                           outlog);
    outlog << "        Agreement with [" << fileid << "] and [" << i << "]: " << std::setprecision(3) << std::setw(10) << phase_diff_agreement << std::setw(17) << "Origin shift: " << origin_shift.format() << std::endl;

    phase_diff_agreement_with_inversion = phase_diff_agreement;
    // Now assess agreement with inversion, but only if the spacegroup is achiral
    if (!space_group_is_chiral)
    {
    invert = true;
    ipa_functions::compute_phase_agreement(input_fp[i],
                                           input_fp[fileid], // This is the phase set that is inverted.
                                           zero_phase_weight,
                                           settings.res_bins,
                                           settings.allow_origin_shifts,
                                           invert,
                                           phase_diff_agreement_with_inversion,
                                           phase_diff_centric,
                                           phase_diff_acentric,
                                           map_cc,
                                           origin_shift_with_inversion,
                                           outlog);
    
    
    // Record the smallest mean phase difference, and note if inversion was required to generate this difference
    outlog << "Inverse agreement with [" << fileid << "] and [" << i << "]: " << std::setprecision(3) << std::setw(10) << phase_diff_agreement_with_inversion << "Origin shift: " << origin_shift_with_inversion.format() << std::endl;
    }
    if (phase_diff_agreement <= phase_diff_agreement_with_inversion) // Smallest phase difference was generated without inversion (or inversion wasn't performed).
    {
      // Write all relevant values to the matrices:
      agreement_matrix.push_back(phase_diff_agreement);
      shift_coords.push_back(origin_shift);
      inverse_flags.push_back(false);
    }
    else // Smallest phase difference was generated with inversion
    {
      // Write all relevant values to the matrices:
      agreement_matrix.push_back(phase_diff_agreement_with_inversion);
      shift_coords.push_back(origin_shift_with_inversion);
      inverse_flags.push_back(true);
    }
  }
  return true;
}

std::vector<int> consensus_worker::return_neighbours(const double& epsilon, const int fileid, const std::vector<float>& matrix)
{
  /*
   A phase set is a neighbour to the input set, when the distance between then falls below a threshold value (epsilon)
  
   This function takes the fileid of a phase set, and generates a list of its neighbours, in the form of integer lookups into the 1D agreement matrix. 
   
   To identify neighors we need to look at all the *potential* neighbors in the 1D agreement matrix, and check which of those pairings have a distance below the threshold.
  
   To do this we:
   - Skip all entries before the column defined by the fileid 
   - Examine all entries in the column defined by the fileid
   - Examine all entries in the row defined by the fileid 
   Note that the diagonal entries (representing self-agreement) are not stored in the 1D agreement matrix 
   
     0    1    2    3     4  ... n   <- FileIds
0        [0]  [1]  [3]   [6]   [n1]
1             [2]  [4]   [7]   [n2]  <- Integer lookups for triangular agreement matrix.
2                  [5]   [8]   [n3]
3                        [9]   [n4]
4                              [n5]
n
 Practically this means:
  
 Potential Neighbours of 1 = [0,   2,4,7]
 Potential Neighbours of 2 = [1,2,   5,8]
 Potential Neighbours of 3 = [3,4,5,   9]
 Potential Neighbours of 4 = [6,7,8,9   ]
 ...
   */
  
  // Declare new array to hold the neighbour id's
  std::vector<int> neighbours;
  
  //outlog << "Finding neighbours for file " << fileid << std::endl;

  for (int i = 0 ; i < n_total_files; i++)
  {
    if (i < fileid) continue; // skip further tests, these entries precede the relevant column
    
    int startloc = file_xy_get_1d(i); // y = 0 by default, which returns the start location.
    
    // test all potential neighbours in the column associated with fileid
    if (i == fileid)
    {
      //outlog << "Checking my column " << i << "." << std::endl;
      int t = 0;
      while (t < fileid) // Columns always have a length of fileid.
      {
        //eps to test = startloc + t
        if (matrix[startloc+t] <= epsilon)
        {
          //outlog << "Neighbour idenfied within own column, " << t << std::endl;
          neighbours.push_back(t); // The neighbour is the filename id for this row
        }
        t++;
      }
      continue; // Skip further tests.
    }
    
    // test all potential neighbours in the row associated with fileid
    // Only columns > fileid are relevant here.
    if (i > fileid)
    {
      //outlog << "Checking column " << i << "." << std::endl;

      if (matrix[startloc+fileid] <= epsilon) // Test our row for file 'i' column.
      {
        //outlog << "Neighbour identified from other files (" << i << ") column." << std::endl;
        neighbours.push_back(i); // the neighbour is the filename id for this column
      }
      continue;
    }
  }
  
  //outlog << "Found " << neighbours.size() << " neighbours for file " << fileid << std::endl;
  
  return neighbours;
}

void consensus_worker::dbscan(const double &epsilon, const double &minPoints, int &n_clusters, std::vector<int> &cluster_id)
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
  
  int n_elements = n_total_files;//distance[0].size();
  
  cluster_id.clear();
  cluster_id.resize(n_elements, -1); // initialize vector of size n_elements, with integer entries -1, to allow rapid lookup of the cluster assignment
  
  int cluster_counter = 0;
  
  // Loop over all elements
  for(int i = 0; i < n_elements; i++)
  {
    //outlog << "Processing Element: " << i << std::endl;
    
    if ( cluster_id[i] == -1 ) // only execute if this element is currently unclassified
    {
      std::vector<int> neighbors; // create a integer valued vector to hold the neighbors
      
      neighbors = return_neighbours(epsilon, i, agreement_matrix); // Populate neighbours, not including self.
      
      if ( (neighbors.size() + 1)  < minPoints)  // Fails the core point condition
      {
        cluster_id[i] = -2;                   // Label this element as a noise point
      }
      else                                    // Satisfies the Core point condition, so now we expand the cluster
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
            neighbors = return_neighbours(epsilon, seed_set[CountSeeds], agreement_matrix);
            neighbors.push_back(seed_set[CountSeeds]); // include J (self) as a neighbour.
            // This list now *includes* the query point.
            
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

void consensus_worker::matrix_1d_to_filexy(const int i, int& x, int& y)
{
  /*
   This functions take a position (i) in the 1D agreement matrix, and evaluates the x (column) and y (row) positions, which represent the file ID's that give rise to the agreement at position (i) 
   */
  x = 0;
  y = 0;
  int startloc;
  for (int n = 0; n < n_total_files; n++)
  {
    startloc = file_xy_get_1d(n); // Find the startloc of file n,
    if (i < startloc) // i is less than startloc, :. n = the file after i;
    {
      x = n - 1; // our first file is the file before n
      break;
    }
  }
  
  // we can get the new startloc of x, knowing the startloc of x+1 with the following formula:
  int x_startloc = startloc - x; // all columns are of length of fileid.
  
  y = i - x_startloc; // our second file is the difference between n and the startloc of our file.
}

int consensus_worker::file_xy_get_1d(const int x, const int y)
{
  /*
   x maps to columns, y maps to rows
   This function takes a fileid (x), and a fileid (y), and returns the corresponding position in the 1D agreement matrix 
   
   Matrices will look like this: for files 0,1,2,3,4, through to n...
        0    1    2    3     4  ... n
   0        [0]  [1]  [3]   [6]   [n1]
   1             [2]  [4]   [7]   [n2]
   2                  [5]   [8]   [n3]
   3                        [9]   [n4]
   4                              [n5]
   n
  
   The diagonals are technically irrelevent, as maps are not compared with themselves... so therefore only need to populate a matrix from item 1 (x) for 0 (y) and onwards...
   
   Our matrix should just be a 1d array... where the lookup is:
   file n will read n values, However starts n starts at 0, so nth file is n-1th if we use an adjusted "nth triangular" forumla:  ( ( n ( n + 1 )) / 2 ) - n  ... therefore
   file 1, reads 1, starting from pos 0 ((( 1 ( 1+1 )) / 2) - 1)
   file 2, reads 2, starting from pos 1
   file 3, reads 3, starting from pos 3
   file 4, reads 4, starting from pos 6
  
   seems to be working, and i'm not going to do a proof to infinity.
   
   */
  
  int n;
  
  n = ((x*(x + 1))/2)-x; // start location
  n += y; // Append row number.
  
  return n;
}

bool consensus_worker::initialise()
{
  /*
   This function reads in some crystallographic information (cell, spacegroup, resolution) from the first file processed, this should only need to be called once. 
   The first phase set encountered will dictate whether subsequent phase sets are accepted, as they should all have the same base crystallographic information
   */
  // Create custom log file for this phasing_consensus run.
  char log_filename[256];
  snprintf(log_filename, sizeof(log_filename), "%sphase_consensus_log.txt", settings.working_dir.c_str());
  outlog.open(log_filename, std::ios_base::app); // we do this so we never overwrite previous logs. Probably a good idea for consensus phase only, as its run multiple times.
  if (!outlog.good())
  {
    std::cerr << "Bad phase consensus log file generation!";
    return false;
  }

  util::version_fingerprint(outlog);
  
  outlog << "Phase Consensus Log file" << std::endl;
  outlog << "Starting phase file list count = " << phase_file_list.size() << std::endl;
  
  outlog << "Reading in space group information from first mtz file." << std::endl;
  // Read the first phases from mtz into hkl (member variable)
  clipper::CCP4MTZfile mtzin;
  mtzin.open_read( phase_file_list[0] );
  mtzin.import_hkl_info(hkl); // read sg, cell, reso
  mtzin.close_read();
  
  // Record chirality of space group.
  space_group_is_chiral = util::space_group_is_chiral(hkl.spacegroup().spacegroup_number());
  
  if (space_group_is_chiral) {outlog << "Space group is chiral" << std::endl;}
  else                       {outlog << "Space group is achiral" << std::endl;}
  
  if (settings.apodization) // apodization is being applied and weights are non-unitary 
  {
  outlog << "Calculating the weights associated with the apodization function" << std::endl;
  outlog << "These weights will be used when calculating mean absolute phase differences" << std::endl;  
  }
  
  // Initialise the object holding the weights.
  zero_phase_weight.init(hkl,hkl.cell());
  float sigma = 1.0/(2.4477*settings.effective_resolution_limit);
  float sigma_squared = std::pow(sigma , 2);
  clipper::HKL_data_base::HKL_reference_index ih;
  for ( ih = hkl.first(); !ih.last(); ih.next() )
  {
    if (settings.apodization)
    {
      float s_squared = ih.invresolsq();
      float weight = exp(-s_squared/(2*sigma_squared));
      
      zero_phase_weight[ih].phi() = 0.0;
      zero_phase_weight[ih].fom() = weight;
    }
    else
    {
      zero_phase_weight[ih].phi() = 0.0;
      zero_phase_weight[ih].fom() = 1.0;
    }
  }
  
  outlog << "Performing check on settings" << std::endl;
  // Some basic checks for the settings.
  if ( settings.threshold_specified && ((settings.epsilon_input >= 90.0) || (settings.epsilon_input <= 0.0)) )
  {
    outlog << "Value for epsilon (°) should be > 0 and < 90\n";
    return false;
  }
  
  if ( settings.min_pts < 1 )
  {
    outlog << "Minimum number of points cannot be < 1\n";
    return false;
  }
  
  // if an known phase set is input, read it here
  
  if ( settings.input_known_phases && settings.input_filename_known_phase_set != "NONE" && settings.known_phase_set_fp_column_labels != "NONE")
  {
    // Currently commented, pending investigation of bug which is causing incorrect reporting of statistics in this setting
    // outlog << "An known phase set was specified in the experiment parameter file... importing the known phase data for testing purposes ..." << std::endl;
    
    // Read the known phases
    mtzin.open_read( settings.input_filename_known_phase_set ); // open the mtz file for reading
    mtzin.import_hkl_info( hkl_known ); // read sg, cell, reso
    
    if (hkl_known.spacegroup().spacegroup_number() != hkl.spacegroup().spacegroup_number())
    {
      outlog << "The space group of the known phase set does not match the space group of the target structure. This is likely a user input error, and the incorrect file has been passed in the settings." << std::endl;
      return false;
    }
    /*
    // TODO - reenable this check with some reasonable tolerance. As it is, insignificant differences in the resolution  limit will trigger this message, which is probably not good ?
    // There is some clipper functionality for this ...
    if (hkl_known.resolution().limit() != hkl.resolution().limit() )
    {
      double limitd = hkl.resolution().limit();
      double limitk = hkl_known.resolution().limit();
      outlog << "The resolution limit of the known phase set, and the target data differs: " << limitk << " (known) is not equal to " << limitd  << " (target)." << std::endl;
    }
    */
    //outlog << "Importing known data with labels: " << settings.col_fp_known << " " << std::endl;
    fp_known.init(hkl_known,hkl_known.cell());
    mtzin.import_hkl_data( fp_known, settings.known_phase_set_fp_column_labels ); // Here are the known amplitudes and phases
    mtzin.close_read();
    
  }
  else
  {
    // Set this to false so things aren't confused, really this shouldn't be much of an issue.
    settings.input_known_phases = false;
  }
  
  is_init = true; // No need to do any of this ever again!
  
  return true;
}

void consensus_worker::print_agreement_matrix_to_log(const std::vector<float>& matrix, std::ofstream& outlog)
{
  /*
   Helper function that will print the current matrix to the log file in a nice format.
   Should work for bools and floats..... hopefully.
   */
  outlog << "\n Agreement Matrix for currently processed files:" << std::endl;
  outlog << "   FILE:";
  for (int i = 0; i < n_total_files; i++)
  {
    outlog << " [" << std::setw(3) << i << "]" << std::fixed;
  }
  outlog << std::endl;
  
  for (int y = 0; y < n_total_files; y++)
  {
    outlog << std::setprecision(1) << "   [" <<
    std::setw(3) << y << "]" << std::fixed;
    
    for (int x = 0; x < n_total_files; x++)
    {
      if (x > y)
      {
        //Inside triangle, lookup matrix value.
        int j = file_xy_get_1d(x,y);
        
        if (matrix[j] < settings.epsilon_input)
        {
          //std::cout<<"\033[31;1;4mHello\033[0m"
          outlog << std::setprecision(1) << "*" <<
          std::setw(5) << matrix[j] << std::fixed;
        }
        else
        {
          outlog << std::setprecision(1) << " " <<
          std::setw(5) << matrix[j] << std::fixed;
        }
        
        
      }
      if (x <= y)
      {
        // Off triangle fill with blank space.
        outlog << " " << std::setw(5) << " " << std::fixed;
      }
    }
    outlog << std::endl;
  }
}

void consensus_worker::print_inverse_matrix_to_log(const std::vector<bool>& matrix, std::ofstream& outlog)
{
  /*
   Helper function that will print the current matrix to the log file in a nice format.
   Should work for bools and floats..... hopefully.
   */
  outlog << "\n Inverse Flag Matrix for currently processed files:" << std::endl;
  outlog << "   FILE:";
  for (int i = 0; i < n_total_files; i++)
  {
    outlog << " [" << std::setw(3) << i << "]" << std::fixed;
  }
  outlog << std::endl;
  
  for (int y = 0; y < n_total_files; y++)
  {
    outlog << std::setprecision(1) << "   [" <<
    std::setw(3) << y << "]" << std::fixed;
    
    for (int x = 0; x < n_total_files; x++)
    {
      if (x > y)
      {
        //Inside triangle, lookup matrix value.
        int j = file_xy_get_1d(x,y);
        std::string invertedbool = "  ";
        if (matrix[j] == true) invertedbool = "-1";
        
        outlog << " " <<
        std::setw(4) << invertedbool << " " << std::fixed;
      }
      if (x <= y)
      {
        // Off triangle fill with blank space.
        outlog << " " << std::setw(5) << " " << std::fixed;
      }
    }
    outlog << std::endl;
  }
}

void consensus_worker::log_cluster_info()
{
  outlog << number_of_clusters << " clusters are Identified. \n (cluster values of -2 is equivalent to noise)" << std::endl;
  outlog << "Cluster : FileId";
  if (settings.input_known_phases) {outlog << "(Known)";}
  outlog << std::endl;
  
  // Print out the cluster info real quick:
  for (int i = number_of_clusters; i >= -2; i --)
  {
    for (int fileid = 0; fileid < n_total_files; fileid++)
    {
      if (cluster_ids[fileid] == i)
      {
        //if (settings.input_known_phases) // Was causing segmentation fault when not using an known.
        //{
        //  outlog << std::setw(7) << i << " : " << fileid << " " << phase_file_list[fileid] << " (" << agreement_with_known[fileid] << ")" << std::endl;
        //}
        //else
        //{
        //TODO: Update this to include the known again once that is re-added and fixed.
          outlog << std::setw(7) << i << " : " << fileid << " " << phase_file_list[fileid] << " (" << ")" << std::endl;
        //}
      }
    }
  }
}


float consensus_worker::calculate_agreement_with_known(int file_id)
{
  /*
   Take the file_id of an input to the consensus phase, and populates the array of known phase agreements for reporting later.
   */
  if (!settings.input_known_phases)
  {
    outlog << "Calculating the agreement with a known phase set was requested, however no known phase set was specified" << std::endl;
    return -1.0;
  }
  
  float known_agreement = 0.0;
  
  bool invert = false;
  float phase_diff_agreement = 0.0;
  float phase_diff_centric = 0.0;
  float phase_diff_acentric = 0.0;
  float map_cc = 0.0;
  clipper::Coord_frac origin_shift( 0.0, 0.0, 0.0 );
  
  ipa_functions::compute_phase_agreement(fp_known,
                                         input_fp[file_id],
                                         zero_phase_weight,
                                         settings.res_bins,
                                         settings.allow_origin_shifts,
                                         invert,
                                         phase_diff_agreement,
                                         phase_diff_centric,
                                         phase_diff_acentric,
                                         map_cc,
                                         origin_shift,
                                         outlog, true);
  
  known_agreement = phase_diff_agreement;
  
  if (!space_group_is_chiral)
  {
    invert = true;
    ipa_functions::compute_phase_agreement(fp_known,
                                           input_fp[file_id],
                                           zero_phase_weight,
                                           settings.res_bins,
                                           settings.allow_origin_shifts,
                                           invert,
                                           phase_diff_agreement,
                                           phase_diff_centric,
                                           phase_diff_acentric,
                                           map_cc,
                                           origin_shift,
                                           outlog, false);
    
    if (known_agreement > phase_diff_agreement)
    {
      known_agreement = -1*phase_diff_agreement; // We make the whole thing negative to indicate inversion
    }
  }
  
  if (agreement_with_known.size() == file_id)
  {
    // Perfect! We can push back.
    agreement_with_known.push_back(known_agreement);
  }
  else if (agreement_with_known.size() < file_id)
  {
    agreement_with_known.resize(file_id+1,0.0); // Ensure there is enough room
    agreement_with_known[file_id] = known_agreement; // Before setting the location outright.
  }
  
  return agreement_with_known[file_id];
}

