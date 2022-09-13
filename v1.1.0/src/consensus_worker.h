// IPA: Iterative Projection Algorithms for protein crystallography
// Header for the phase concensus worker class  

/*
 
 This class determines if a cluster of consistent phase sets can be identified among the available "solutions" (i.e. do we have a true solution to the phase problem) 

 It is designed to run in parallel with the phase determination process, such that after each phase determination job completes, the new solution is compared with pre-existing solutions, to check if a cluster exists
 
 In the achiral space groups handedness also needs to be taken into account - ab initio phase determination can result in generation of the true solution, or its inverse.
 
 A single instance of this class will exist in its own thread, as solutions are found, they will be safely appended through an update function one by one. Solutions are written to .mtz files at the conclusion of phase determination, and are communicated between objects via the file name and path.
 
 The single instance of consensus_worker is linked to the lifespan of the ipa_manager which is performing phase determination. Once the ipa_manager completes all phase runs, it will wait for the final consensus jobs to complete before its routines return to the experiment manager.
 
 Importantly, the consensus_worker object itself will not occupy its own thread, but a function to perform the consensus with the added phase solution will be run in parallel. As items are appended in a case by case bases, we should be able to control for any non-thread safe behaviour this way.
 */
#pragma once

#include <iostream>
#include <fstream> // file streams
#include <iomanip> // std::setw & std::setprecison & std::boolalpha
#include <algorithm> // std::sort & std::min
#include <cmath> // lround, lroundf, std::abs, std::pow
#include <time.h> // clock
#include <random> // random number generation

#include <thread> // THIS MIGHT BE REDUNDANT HERE...
#include <atomic> // Needed for updating results of things.
#include <chrono> // Useful for testing threads and for reporting system time.

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

#include "ipa-lib.h"
#include "svn-version.h"

#ifndef CONSENSUS_WORKER
#define CONSENSUS_WORKER

struct member_data
{
  float mean_phase_diff_all_with_known;
  float mean_phase_diff_all_with_known_inverse;
};

struct cluster_data
{
  float  epsilon;
  int cluster_id;
  int n_members; // Redundant for a call to list_of_members.size()?
  std::vector<std::string> list_of_members; // We remember which files went into this cluster, now.
  std::string consensus_filename; // The output average of the members.
  //std::string consensus_filename_inverse = "NONE"; // if space group allows, this will be populated. else its None.
  std::string consensus_filename_inverse; // Inverse is always generated now ... if the space group is chiral, we switch to the enantiomorphic space group when we invert
  
  float scv;
  float mean_phase_difference_across_cluster;
  
  float mean_phase_diff_all_with_known;
  float mean_phase_diff_all_with_known_inverse;
};

class consensus_worker
{
public:
  // ------------------------------------------------------- //
  // ---------------- PUBLIC Methods ----------------------- //
  // ------------------------------------------------------- //
  
  bool add_phase_input_to_list(std::string solution_filename); // Add solution "dir/filename.mtz" to be included in experiment.
  
  std::ostream& return_list_of_consensus_solutions(std::vector<cluster_data>& solutions, std::ostream& outlog); // Return information for both the filenames of the solutions, and the output of the current state.
  
  bool check_for_consensus_phases(std::atomic<float>& MyProgress); // Run the consensus experiment (requires at least 1 filename to be added) Best run in its own thread.
  
  consensus_worker();
  consensus_worker(consensus_phase_settings new_settings); // Constructor alternative: settings passed in directly.
  ~consensus_worker(); // Deconstructor.
  
  // ------------------------------------------------------- //
  // ---------------- PUBLIC Variables --------------------- //
  // ------------------------------------------------------- //
  
  std::atomic<bool> is_idle; // Simple check for thread safety, ensures only one thread of calculating phases ever occurs at a time. This is particularly important as we only add to the list of files to process when the worker is idle. Else issues.
  
  consensus_phase_settings settings; // Settings.
  
private:
  
  // ------------------------------------------------------- //
  // ---------------- PRIVATE Methods ---------------------- //
  // ------------------------------------------------------- //
  
  //void matrix_2d_to_1d(int& i,const int x, const int y = 0); // y = 0 as, not including the 7 lookup simply gives you the startlocation of the column. Columns always go for the length of the file id.
  
  void matrix_1d_to_filexy(const int i, int& x, int& y); //Converts a point in the 1D distance array to the file id's of the two phase sets that contributed to that point.
  
  int file_xy_get_1d(const int x, const int y = 0); // Get 1d matrix location of the comparison between file x and y. The reverse of matrix_1d_to_filexy().
  
  std::vector<int> return_neighbours(const double& epsilon, const int fileid, const std::vector<float>& matrix);
  
  // Possibly rename this dbscan_1D to avoid confusion with the "2d" version in the main library ?
  void dbscan(const double &epsilon, const double &minPoints, int &n_clusters, std::vector<int> &cluster_id); // Executes the DBSCAN algorithm using a 1D distance matrix. Returns cluster_id (a list detailing which cluster each file belongs to), and the total number of clusters.
  
  void log_cluster_info(); // Helper function to print to log.
  
  bool initialise(); // Called once, initialises the object properly after the first phase is given.
  
  void print_agreement_matrix_to_log(const std::vector<float>& matrix, std::ofstream& outlog);
  
  void print_inverse_matrix_to_log(const std::vector<bool>& matrix, std::ofstream& outlog);
  
  bool input_data_from_filename(int fileid);
  
  bool calculate_and_append_data_to_matrix(int fileid);
  
  bool generate_consensus_phases();
  
  float calculate_agreement_with_known(int file_id);
  
  // ------------------------------------------------------- //
  // ---------------- PRIVATE Variables -------------------- //
  // ------------------------------------------------------- //
  
  bool is_init = false; // flag for initialisation.
  std::ofstream outlog; // Log file var.
  
  // ----------- Crystallographic data constants ------------//
  
  clipper::HKL_info hkl;

  bool space_group_is_chiral; 
    
  // space group chirality affects how we test for phase agreement, and the role of inversion
  //
  // If the space group is achiral ... 
  //   We need to allow for inversion when testing agreement of phase sets, since the hand of the reconstruction is arbitrary
  //   We need export both the solution and its inverse in the original space group (space group doesn't change under inversion)
  // If the space group is chiral ... 
  //    We don't need to allow for inversion when testing agreement of phase sets, since the hand of the reconstruction is fixed by the space group choice.
  //    However unless we are working with a known structure, the space group choice is arbitrary - it could be either of a pair of enantiomorphs.
  //    We need to export both the solution and its inverse, switching to the enantiomorphic space group when we invert (space group changes under inversion)
  
  clipper::HKL_data<clipper::data64::Phi_fom> zero_phase_weight; // phase weights.
  
  // Known Phases (for testing purposes)
  clipper::HKL_info hkl_known;
  clipper::HKL_data<clipper::data64::F_phi> fp_known;
  
  // ----------- Phase file management ----------------------//
  int n_total_files; //  Quick Ref: Current total number of files i.e. phase_file_list.size()
  std::vector<std::string> phase_file_list; // Vector list of Files to be processed.
  std::vector<std::string> solutions_file_list; // Vector list of solutions found and written.
  std::vector<float> agreement_with_known; // Vector list of each files agreement to the known phase set (is negative if inversion was required)

  std::vector<bool> files_processed_flag; // Flags each file to know if its been calculated within the concensus or not.
  std::vector< clipper::HKL_data<clipper::data64::F_phi> > input_fp; // Vector of input phases, should be the same size as n_total_files.
  
  // ----------- Matrices -----------------------------------//
  int matrix_size = 0; // We cannot add a file to the matrix unless we want to overwrite columns.
  
  std::vector<float> agreement_matrix; // A vector describing the matrix containing the pairwise phase differences, initialized to zero. The matrix is symmetrical about the diagonal, and no diagonal or other side is stored.
  std::vector<float> agreement_matrix_sorted; // The same matrix, however we have sorted the distances in order, we can no longer determine which files came from which values however. We could make a map if needed to retain that info?
  std::vector<bool> inverse_flags; // flags for whether the dataset required to be inverted or not for maximum agreement.
  std::vector<clipper::Coord_frac> shift_coords; // Coordinates to align best phase agreements. We can have a 1D version of this.
  
  // ----------- Consensus clustering -----------------------//
  int number_of_clusters;          // How many clusters have been identified.
  std::vector<int> cluster_ids;    // Which cluster each file belongs to.
  
  
  // ----------- Consensus Averaging ------------------------//
  std::vector< clipper::HKL_data<clipper::data64::F_phi> > cluster_averages;
  
  std::vector<cluster_data> all_prototypes;
  
  // ------------------------------------------------------- //
  // ------------------------- NOTES ----------------------- //
  // ------------------------------------------------------- //

  
  /*
   Matrices will look like this: for files 0,1,2,3,4, through to n...
        0    1    2    3     4  ... n
   0        [0]  [1]  [3]   [6]   [n1]
   1             [2]  [4]   [7]   [n2]
   2                  [5]   [8]   [n3]
   3                        [9]   [n4]
   4                              [n5]
   n
  
   The diagonals are technically irrelevent, as maps are not compared with themselves... so therefore only need to populate a matrix from item 1 (x) for 0 (y) and onwards...
   
   // Our matrix should just be a 1d array... where the lookup is:
   file n will read n values, However starts n starts at 0, so nth file is n-1th if we use an adjusted "nth triangular" forumla:  ( ( n ( n + 1 )) / 2 ) - n  ... therefore
   file 1, reads 1, starting from pos 0 ((( 1 ( 1+1 )) / 2) - 1)
   file 2, reads 2, starting from pos 1
   file 3, reads 3, starting from pos 3
   file 4, reads 4, starting from pos 6
  
   seems to be working, and i'm not going to do a proof to infinity.
   
   The reason we may want to switch to a 1D triangle based matrix is because a) its symmetrical around the diagonal anyway, and b) this means we don't need to know the total size of the number of files, and can simply append to the triangle as each new phase input file becomes apparent.
   i.e. for 4 files: i.e. 0,1,2,3, total size of 4, to add the 1d matrix values of the 4th file:
   
   for (int i = 0; i < (n_total_size-1); i ++)
   {
      value = compare_files(n_total_size-1, i)
      matrix_1d.pushback(value);
   
      // appends into positions [6,7,8,9]
   }
   
   Annoyingly, its convenient for the DB_SCAN algorithm to recieve a full 2D matrix, as this makes it easier to quickly lookup a single files distance with all others by iterating over columns, we will write a simple function to populate a 2D matrix from the 1D one, so as to pass this to DB_SCAN with no issues.
   .... ooorrrrr we could rewrite db_scan to work on the 1d matrix... the only "tricky" thing there will be identifing neighbours appropriately...
   
         0    1    2    3     4  ... n
    0        [0]  [1]  [3]   [6]   [n1]
    1             [2]  [4]   [7]   [n2]
    2                  [5]   [8]   [n3]
    3                        [9]   [n4]
    4                              [n5]
    n
   
   Neighbours of 1 = [0,   2,4,7]
   Neighbours of 2 = [1,2,   5,8]
   Neighbours of 3 = [3,4,5,   9]
   Neighbours of 4 = [6,7,8,9   ]
   
   The pattern here is take your own column n, then add the nth term from everyone elses column...
   (iff that file column is greater than n)
   we can do this using the rules from above:
   i.e. ...
   
   int n = target;
   std::vector<int> neighbours;
   for (int i = 0 ; i < n_total_files; i++)
   {
      int startloc = function_get_1d_file_position(i)
      if (i == n)
        t = 0
        while t < n
        neighbours.add(startloc + t)
        continue
      if (i > n)
      
      neighbours.add(startloc + n)
   }
   
   this will populate a neighbours[] array with the integer values of where the neighbours are in the 1d table.
   */
  
};

#endif
