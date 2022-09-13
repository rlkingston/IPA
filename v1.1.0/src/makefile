# Simple makefile for compiling program IPA 

# Compiler Specification:
CC = g++ # C++ Compiler 
FC = gfortran # Fortran compiler 

# Update these paths to reflect the local environment:
FORTRAN_LIB_DIR:=/usr/local/lib/# Directory containing gfortran library
CLIPPER_LIB_DIR:=/Applications/devtools/install/lib/# Directory containing clipper libraries: libclipper-ccp4, libclipper-core, and libclipper-contrib
CLIPPER_HEADER_DIR:=../clipper_src/# Location of the clipper c++ header files clipper/clipper-ccp4.h, clipper/clipper-ccp4.h, etc
EIGEN_DIR:=../eigen_src/# Directory containing the Eigen3 package

# Nothing below should need updating:
B:=../build/# Directory for building common .o objects.
F:=../bin/#Final output Directory
C:=../src/# Directory of local src, redundant for now.

build: $(F)ipa
	echo "IPA build successful!"

$(F)ipa: $(F) $(B) $(B)ipa.o $(B)ipa-lib.o $(B)util-lib.o $(B)ipa_worker.o $(B)consensus_worker.o $(B)ipa_manager.o $(B)exp_manager.o $(B)local_regression.o $(B)gamma_distribution.o $(B)apodization.o $(B)rogers_analysis.o $(B)atom_count.o
	$(CC) -o ../bin/ipa $(B)ipa.o $(B)ipa-lib.o $(B)util-lib.o $(B)ipa_worker.o $(B)consensus_worker.o $(B)ipa_manager.o $(B)exp_manager.o $(B)local_regression.o $(B)gamma_distribution.o $(B)apodization.o $(B)rogers_analysis.o $(B)atom_count.o -std=c++11 -lstdc++ -L $(FORTRAN_LIB_DIR) -lgfortran -L $(CLIPPER_LIB_DIR) -lclipper-ccp4 -lclipper-core -lclipper-contrib -lccp4c -lccp4f -lfftw -lrfftw -pthread
$(B)ipa.o: $(C)ipa.cc
	$(CC) -c $(C)ipa.cc -o $(B)ipa.o -std=c++11 -I $(CLIPPER_HEADER_DIR) -I $(EIGEN_DIR)
$(B)ipa-lib.o: $(C)ipa-lib.cc
	$(CC) -c $(C)ipa-lib.cc -o $(B)ipa-lib.o -std=c++11 -I $(CLIPPER_HEADER_DIR) -I $(EIGEN_DIR)
$(B)util-lib.o: $(C)util-lib.cc $(C)svn-version.h
	$(CC) -c $(C)util-lib.cc -o $(B)util-lib.o -std=c++11 -I $(CLIPPER_HEADER_DIR) -I $(EIGEN_DIR)
$(B)ipa_worker.o: $(C)ipa_worker.cc
	$(CC) -c $(C)ipa_worker.cc -o $(B)ipa_worker.o -std=c++11 -I $(CLIPPER_HEADER_DIR) -I $(EIGEN_DIR)
$(B)consensus_worker.o: $(C)consensus_worker.cc
	$(CC) -c $(C)consensus_worker.cc -o $(B)consensus_worker.o -std=c++11 -I $(CLIPPER_HEADER_DIR) -I $(EIGEN_DIR)
$(B)ipa_manager.o: $(C)ipa_manager.cc
	$(CC) -c $(C)ipa_manager.cc -o $(B)ipa_manager.o -std=c++11 -I $(CLIPPER_HEADER_DIR) -I $(EIGEN_DIR)
$(B)exp_manager.o: $(C)exp_manager.cc
	$(CC) -c $(C)exp_manager.cc -o $(B)exp_manager.o -std=c++11 -I $(CLIPPER_HEADER_DIR) -I $(EIGEN_DIR)
$(B)local_regression.o: $(C)local_regression.f90
	$(FC) -c $(C)local_regression.f90 -o $(B)local_regression.o -J $(B)
$(B)gamma_distribution.o: $(C)gamma_distribution.f90
	$(FC) -c $(C)gamma_distribution.f90 -o $(B)gamma_distribution.o -J $(B)
$(B)apodization.o: $(C)apodization.f90
	$(FC)  -ffree-form -ffree-line-length-0 -std=f2008 -c $(C)apodization.f90 -o $(B)apodization.o -J $(B)
$(B)rogers_analysis.o: $(C)rogers_analysis.f90
	$(FC) -c $(C)rogers_analysis.f90 -o $(B)rogers_analysis.o -J $(B)
$(B)atom_count.o: $(C)atom_count.f90
	$(FC) -c $(C)atom_count.f90 -o $(B)atom_count.o -J $(B)
# The following are just used to ensure the correct directories exist.
$(B):
	mkdir $(B)
$(F):
	mkdir $(F)

clean:
	rm -f $(B)ipa.o $(B)ipa-lib.o $(B)util-lib.o $(B)ipa_worker.o $(B)consensus_worker.o $(B)ipa_manager.o $(B)exp_manager.o $(B)local_regression.o $(B)apodization.o $(B)rogers_analysis.o $(B)atom_count.o $(B)gamma_distribution.o $(F)ipa
	rm -rf $(B)