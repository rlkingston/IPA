module atom_count

  use,intrinsic :: iso_c_binding
  implicit none
  
  
  integer  (c_int), parameter, private      :: n_elements = 98 ! The number of elements over which we count. Controls the dimensioning of atom_counts. A change to this number must be reflected in any C/C++ driver
  
  contains

    subroutine count_protein_atoms(c_string_protein_sequence_filename, n_chains, &
                                  atom_counts, c_string_output_filename) bind(c)
      
      ! Calculates the total number of atoms of each element in a protein, given a user supplied protein sequence
      ! Based very loosly on the Fortran77 program atmcnt, authored by Steven Sheriff
      
      ! Protein sequence (one letter code) must be in plain text format. 
      ! Maximum line length is set by the parameter max_len
      ! White spaces are ignored, any other symbols encountered that are not those of the canonical 20 amino acids will cause the subroutine to terminate
      
      use,intrinsic :: iso_c_binding
      implicit none
      type       (c_ptr),  intent(in), target                     :: c_string_protein_sequence_filename  ! the name of the file containing the protein sequence information. Note that we are defining a pointer to a C Char type variable here
      integer    (c_int),  intent(in)                             :: n_chains ! The number of individual polypeptide chains present in the sequence file
      integer    (c_int),  intent(out),   dimension(n_elements)   :: atom_counts  ! Atom counts for each element
      !integer    (c_int),  intent(in)                             :: output_unit ! Fortran unit for output 
      type       (c_ptr),  intent(in), target                     :: c_string_output_filename  ! the filename for the output. Note that we are defining a pointer to a C Char type variable here
      
      integer :: i
            
      character(len=128), pointer :: f_string_protein_sequence_filename
      integer :: string_length_protein_sequence_filename
      
      character(len=128), pointer :: f_string_output_filename
      integer :: string_length_output_filename

      integer :: input_unit = 15
      integer :: output_unit = 20

      integer, parameter :: max_len = 256 
      character(len=max_len) :: line

      character(len=1), allocatable :: aa_sequence(:) 
      integer :: sequence_length
      integer :: n_peptide_bonds

      integer :: ia, il, ie
      integer :: istat

      integer, parameter :: n_elements_protein = 5 
      integer, parameter, dimension(n_elements_protein) :: atomic_number =  [ 1, 6, 7, 8, 16] 
      character(len=2), parameter, dimension(n_elements_protein) :: element_symbol =  [ "H ", "C ", "N ", "O ", "S "] 
      !real, parameter, dimension(n_elements_protein) :: atomic_mass = [1.00784, 12.0107, 14.0067, 15.999, 32.065] 

      integer, parameter :: n_amino_acid_types = 20   
      integer :: aa_composition(n_amino_acid_types)
      
          
      type  amino_acid_data
        character (len=1) :: single_letter_code
        character (len=3) :: triple_letter_code
        integer, dimension(n_elements_protein) :: n_atoms
      end type

      type (amino_acid_data), dimension(n_amino_acid_types) :: amino_acid
      
      amino_acid(1)%triple_letter_code = "Ala"
      amino_acid(1)%single_letter_code = "A"
      amino_acid(1)%n_atoms = [7,3,1,2,0]
      
      amino_acid(2)%triple_letter_code =  "Arg"
      amino_acid(2)%single_letter_code = "R"
      amino_acid(2)%n_atoms = [14,6,4,2,0]
      
      amino_acid(3)%triple_letter_code =  "Asn"
      amino_acid(3)%single_letter_code = "N" 
      amino_acid(3)%n_atoms = [8,4,2,3,0]

      amino_acid(4)%triple_letter_code =  "Asp" 
      amino_acid(4)%single_letter_code = "D" 
      amino_acid(4)%n_atoms = [7,4,1,4,0]

      amino_acid(5)%triple_letter_code =  "Cys"
      amino_acid(5)%single_letter_code = "C" 
      amino_acid(5)%n_atoms = [7,3,1,2,1]

      amino_acid(6)%triple_letter_code =  "Gln" 
      amino_acid(6)%single_letter_code = "Q" 
      amino_acid(6)%n_atoms = [10,5,2,3,0]

      amino_acid(7)%triple_letter_code =  "Glu"  
      amino_acid(7)%single_letter_code = "E" 
      amino_acid(7)%n_atoms = [9,5,1,4,0]

      amino_acid(8)%triple_letter_code =  "Gly"  
      amino_acid(8)%single_letter_code = "G" 
      amino_acid(8)%n_atoms = [5,2,1,2,0]

      amino_acid(9)%triple_letter_code =  "His"  
      amino_acid(9)%single_letter_code = "H" 
      amino_acid(9)%n_atoms = [9,6,3,2,0]

      amino_acid(10)%triple_letter_code =  "Ile" 
      amino_acid(10)%single_letter_code = "I" 
      amino_acid(10)%n_atoms = [13,6,1,2,0]

      amino_acid(11)%triple_letter_code =  "Leu"  
      amino_acid(11)%single_letter_code = "L" 
      amino_acid(11)%n_atoms = [13,6,1,2,0]

      amino_acid(12)%triple_letter_code =  "Lys" 
      amino_acid(12)%single_letter_code = "K" 
      amino_acid(12)%n_atoms = [14,6,2,2,0]

      amino_acid(13)%triple_letter_code =  "Met" 
      amino_acid(13)%single_letter_code = "M" 
      amino_acid(13)%n_atoms = [11,5,1,2,1]

      amino_acid(14)%triple_letter_code =  "Phe" 
      amino_acid(14)%single_letter_code = "F" 
      amino_acid(14)%n_atoms = [11,9,1,2,0]

      amino_acid(15)%triple_letter_code =  "Pro"  
      amino_acid(15)%single_letter_code = "P" 
      amino_acid(15)%n_atoms = [9,5,1,2,0]

      amino_acid(16)%triple_letter_code =  "Ser"  
      amino_acid(16)%single_letter_code = "S" 
      amino_acid(16)%n_atoms = [7,3,1,3,0]

      amino_acid(17)%triple_letter_code =  "Thr"  
      amino_acid(17)%single_letter_code = "T" 
      amino_acid(17)%n_atoms = [9,4,1,3,0]

      amino_acid(18)%triple_letter_code =  "Trp"  
      amino_acid(18)%single_letter_code = "W" 
      amino_acid(18)%n_atoms = [12,11,2,2,0]

      amino_acid(19)%triple_letter_code =  "Tyr"  
      amino_acid(19)%single_letter_code = "Y" 
      amino_acid(19)%n_atoms = [11,9,1,3,0]

      amino_acid(20)%triple_letter_code =  "Val" 
      amino_acid(20)%single_letter_code = "V" 
      amino_acid(20)%n_atoms = [11,5,1,2,0]
      
      ! Passing a character string from C to Fortran is only straighforward when the the string is of length 1 !!
      ! There are some fixes to this problem in the Fortran 2018 standard, but they have not yet been widely implemented
      ! Here we adopt a solution based on pointers which I dislike, but at least it works. 
      
      ! Do some work to get some filenames we can use  ...
      call c_f_pointer(c_loc( c_string_protein_sequence_filename ), f_string_protein_sequence_filename) ! convert c pointer to fortran pointer ... c_loc obtains the c address of the argument.
      string_length_protein_sequence_filename = index(f_string_protein_sequence_filename, c_null_char) - 1 ! The text that we need can be found in f_string_protein_sequence_filename(1:string_length_protein_sequence_filename). Hurrah.      
 
      call c_f_pointer(c_loc( c_string_output_filename ), f_string_output_filename) ! convert c pointer to fortran pointer ... c_loc obtains the c address of the argument.
      string_length_output_filename = index(f_string_output_filename, c_null_char) - 1 ! The text that we need can be found in f_string_output_filename(1:string_length_ouput_filename). Hurrah.      
 
      open(output_unit, file=f_string_output_filename(1:string_length_output_filename), &
           status="old", position="append", action="write") ! Open the output file
          
      write(output_unit,*) 
      write(output_unit,*) "Calculating atomic composition, based on the protein sequence provided"
      write(output_unit,*)

      open(input_unit,status="OLD",file=f_string_protein_sequence_filename(1:string_length_protein_sequence_filename))    ! open the file containing the protein sequence

      ! run through the sequence and compute the amino acid composition

      sequence_length = 0
      aa_composition(:) = 0
      line(:) = ""


      sequence_read: do 

        read(input_unit,iostat=istat,fmt="(A)") line
        if (istat /= 0) exit

        line_loop: do il = 1, max_len-1
          
          aa_loop: do ia =1,n_amino_acid_types
            if  ( line(il:il) == amino_acid(ia)%single_letter_code ) then
              sequence_length = sequence_length + 1
              aa_composition(ia) = aa_composition(ia) + 1
              exit
            end if
          end do aa_loop 

!         if we have run through all the amino acid types and what we have is not a blank space, that's bad
!         Note that at the end of a complete run through the aa_loop,  ia == n_amino_acid_types+1.

          if ((ia == (n_amino_acid_types+1)) .and. (line(il:il) /= " ")) then
            write(output_unit,*) "I don't know what to do with this character : ", line(il:il)
            write(output_unit,*) "Stopping"
            stop
          end if

        end do line_loop

      end do sequence_read

      close(input_unit)    ! close the file containing the protein sequence
      
      write(output_unit,*) 
      write(output_unit,*) "Total Sequence Length :", sequence_length
      write(output_unit,*) 
      
      write(output_unit,*) 
      write(output_unit,*) "Amino Acid Composition: "
      write(output_unit,*) 

      do ia=1,n_amino_acid_types
        write(output_unit,*) amino_acid(ia)%single_letter_code," ", amino_acid(ia)%triple_letter_code," ", aa_composition(ia)
      end do
    
      ! Calculate total atom counts from amino acid composition
      
      do ia=1,n_amino_acid_types
        do ie=1,n_elements_protein
          atom_counts( atomic_number(ie) ) = atom_counts( atomic_number(ie) ) + aa_composition(ia) * amino_acid(ia)%n_atoms(ie)
        end do
      end do

      write(output_unit,*) 
      write(output_unit,*) "Raw atom counts, computed by summing the contributions of the individual amino acids." 
      write(output_unit,*)
      
      do ie=1,n_elements_protein
        write(output_unit,*) atomic_number(ie), element_symbol(ie), atom_counts( atomic_number(ie) )
      end do

      ! Correct for peptide bond formation - two Hydrogen atoms and one oxygen atom eliminated for each peptide bond formed
      
      n_peptide_bonds = sequence_length - n_chains
      
      atom_counts( 1 ) = atom_counts( 1 ) - 2*n_peptide_bonds ! Correction for Hydrogen, atomic number 1
      atom_counts( 8 ) = atom_counts( 8 ) - 1*n_peptide_bonds ! Correction for Oxygen,   atomic number 8
       
      write(output_unit,*) 
      write(output_unit,*) "Atom Counts, Corrected For Peptide Bond Formation"
      write(output_unit,*) "Assuming ", n_chains, "polypeptide chains are present in the sequence file"  
      write(output_unit,*)
      
      do ie=1,n_elements_protein
        write(output_unit,*) atomic_number(ie), element_symbol(ie), atom_counts( atomic_number(ie) )
      end do

      close(output_unit) ! close the output file
      
    end subroutine count_protein_atoms

    subroutine count_na_atoms(c_string_na_sequence_filename, is_dna, n_chains, &
                              atom_counts, c_string_output_filename) bind(c)
      
      ! Calculates the total number of atoms of each element in dna or rna, given a user supplied nucleic acid sequence
      
      ! Nucleic acid sequence (one letter code) must be in plain text format, case is unimportant
      ! Maximum line length is set by the parameter max_len
      ! White spaces are ignored, any other symbols encountered that are not those of the canonical 5 nucleosides used in DNA or RNA will cause the subroutine to terminate
      
      use,intrinsic :: iso_c_binding
      implicit none
      type       (c_ptr),  intent(in), target                     :: c_string_na_sequence_filename  ! the name of the file containing the nucleic acid sequence information. Note that we are defining a pointer to a C Char type variable here
      integer    (c_int),  intent(in)                             :: n_chains ! The number of individual nucleic acid chains present in the sequence file
      logical    (c_bool), intent(in)                             :: is_dna ! If true, the sequence supplied is DNA. If false, the sequence supplied is RNA
      integer    (c_int),  intent(out),   dimension(n_elements)   :: atom_counts  ! Atom counts for each element
      !integer    (c_int),  intent(in)                             :: output_unit ! Fortran unit for output 
      type       (c_ptr),  intent(in), target                     :: c_string_output_filename  ! the filename for the output. Note that we are defining a pointer to a C Char type variable here
      
      integer :: i
            
      character(len=128), pointer :: f_string_na_sequence_filename
      integer :: string_length_na_sequence_filename

      character(len=128), pointer :: f_string_output_filename
      integer :: string_length_output_filename

      integer :: input_unit = 15
      integer :: output_unit = 20

      integer, parameter :: max_len = 256 
      character(len=max_len) :: line

      character(len=1), allocatable :: na_sequence(:) 
      integer :: sequence_length
      integer :: n_phosphodiester_bonds

      integer :: ib, il, ie
      integer :: istat
      character(len=1) :: current_symbol
      
      integer, parameter :: n_elements_nucleic_acid = 5 
      integer, parameter, dimension(n_elements_nucleic_acid) :: atomic_number =  [ 1, 6, 7, 8, 15] 
      character(len=2), parameter, dimension(n_elements_nucleic_acid) :: element_symbol =  [ "H ", "C ", "N ", "O ", "P "] 
      !real, parameter, dimension(n_elements_nucleic_acid) :: atomic_mass = [1.00784, 12.0107, 14.0067, 15.999, 30.9738] 

      integer, parameter :: n_nucleosides = 5   
      integer :: na_composition(n_nucleosides)
      
          
      type  na_data
        character (len=1) :: single_letter_code
        integer, dimension(n_elements_nucleic_acid) :: n_atoms_dna
        integer, dimension(n_elements_nucleic_acid) :: n_atoms_rna
      end type

      type (na_data), dimension(n_nucleosides) :: nucleoside
      
      nucleoside(1)%single_letter_code = "A" ! deoxyAdenosine / Adenosine
      nucleoside(1)%n_atoms_dna = [13,10,5,3,0] !  Adenine linked to 2'-deoxyribofuranose 
      nucleoside(1)%n_atoms_rna = [13,10,5,4,0] !  Adenine linked to        ribofuranose 

      nucleoside(2)%single_letter_code = "C" ! deoxyCytidine / Cytidine
      nucleoside(2)%n_atoms_dna = [13,9,3,4,0] !  Cytosine linked to 2'-deoxyribofuranose 
      nucleoside(2)%n_atoms_rna = [13,9,3,5,0] !  Cytosine linked to        ribofuranose 

      nucleoside(3)%single_letter_code = "G" ! deoxyGuanosine / Guanosine 
      nucleoside(3)%n_atoms_dna = [13,10,5,4,0]  !  Guanine linked to 2'-deoxyribofuranose 
      nucleoside(3)%n_atoms_rna = [13,10,5,5,0]  !  Guanine linked to         ribofuranose 

      nucleoside(4)%single_letter_code = "T"  ! deoxyThymidine / Thymidine
      nucleoside(4)%n_atoms_dna = [14,10,2,5,0]  !  Thymine linked to 2'-deoxyribofuranose  
      nucleoside(4)%n_atoms_rna = [14,10,2,6,0]  !  Thymine linked to         ribofuranose - not generally found in RNA of course

      nucleoside(5)%single_letter_code = "U"  ! dexoyUridine / Uridine
      nucleoside(5)%n_atoms_dna = [12,9,2,5,0] ! Uracil linked to 2'-deoxyribofuranose  - not generally found in DNA of course
      nucleoside(5)%n_atoms_rna = [12,9,2,6,0] ! Uracil linked to         ribofuranose 


      ! Passing a character string from C to Fortran is only straighforward when the the string is of length 1 !!
      ! There are some fixes to this problem in the Fortran 2018 standard, but they have not yet been widely implemented
      ! Here we adopt a solution based on pointers which I dislike, but at least it works. 
      
      ! Do some work to get some filenames we can use  ...

      call c_f_pointer(c_loc( c_string_na_sequence_filename ), f_string_na_sequence_filename) ! convert c pointer to fortran pointer ... c_loc obtains the c address of the argument.
      string_length_na_sequence_filename = index(f_string_na_sequence_filename, c_null_char) - 1 ! The text that we need can be found in f_string_na_sequence_filename(1:string_length_na_sequence_filename). Hurrah.      

      call c_f_pointer(c_loc( c_string_output_filename ), f_string_output_filename) ! convert c pointer to fortran pointer ... c_loc obtains the c address of the argument.
      string_length_output_filename = index(f_string_output_filename, c_null_char) - 1 ! The text that we need can be found in f_string_output_filename(1:string_length_ouput_filename). Hurrah.      

      open(output_unit, file=f_string_output_filename(1:string_length_output_filename), &
           status="old", position="append", action="write") ! Open the output file
      
      if (is_dna) then      
        write(output_unit,*) 
        write(output_unit,*) "Calculating atomic composition, based on the DNA sequence provided"
        write(output_unit,*)
      else
        write(output_unit,*) 
        write(output_unit,*) "Calculating atomic composition, based on the RNA sequence provided"
        write(output_unit,*)
      end if        
      
      open(input_unit,status="OLD",file=f_string_na_sequence_filename(1:string_length_na_sequence_filename))    ! open the file containing the nucleic acid sequence

      ! run through the nucleic acid sequence and compute the nucleoside composition
      ! Since nucleic acid sequences are commonly given in lower case, we need to do a case conversion to catch this possibility
      
      ! TODO Should we throw a warning if supplied DNA sequence has symbol "U", or supplied RNA sequence has symbol "T"? Probably
      
      sequence_length = 0
      na_composition(:) = 0
      line(:) = ""

      sequence_read: do 

        read(input_unit,iostat=istat,fmt="(A)") line
        if (istat /= 0) exit

        line_loop: do il = 1, max_len-1
          
          na_loop: do ib =1,n_nucleosides

            current_symbol = to_upper( line(il:il) ) ! convert to upper case representation
            
            if  ( current_symbol == nucleoside(ib)%single_letter_code ) then
              sequence_length = sequence_length + 1
              na_composition(ib) = na_composition(ib) + 1
              exit
            end if
          end do na_loop 

!         if we have run through all the nucleosides and what we have is not a blank space, that's bad
!         Note that at the end of a complete run through the na_loop,  ib == n_nucleosides+1.

          if ((ib == (n_nucleosides+1)) .and. (current_symbol /= " ")) then
            write(output_unit,*) "I don't know what to do with this character : ", line(il:il)
            write(output_unit,*) "Stopping"
            stop
          end if

        end do line_loop

      end do sequence_read

      close(input_unit)    ! close the file containing the nucleic acid sequence
      
      write(output_unit,*) 
      write(output_unit,*) "Total Sequence Length :", sequence_length
      write(output_unit,*) 
      
      write(output_unit,*) 
      write(output_unit,*) "Nucleoside Composition: "
      write(output_unit,*) 

      do ib=1,n_nucleosides
        if ( .not.(is_dna) .and. (ib == 4) ) cycle ! No Thymidine    in RNA
        if (      (is_dna) .and. (ib == 5) ) cycle ! No deoxyUridine in DNA
        write(output_unit,*) nucleoside(ib)%single_letter_code," ", na_composition(ib)
      end do

      ! Calculate total atom counts from the nucleoside composition
      
      do ib=1,n_nucleosides
        do ie=1,n_elements_nucleic_acid
          if (is_dna) then
            atom_counts( atomic_number(ie) ) = &
            atom_counts( atomic_number(ie) ) + na_composition(ib) * nucleoside(ib)%n_atoms_dna(ie)
          else
            atom_counts( atomic_number(ie) ) = &
            atom_counts( atomic_number(ie) ) + na_composition(ib) * nucleoside(ib)%n_atoms_rna(ie) 
          end if           
        end do
      end do

      write(output_unit,*) 
      write(output_unit,*) "Raw atom counts, computed by summing the contributions of the individual nucleosides" 
      write(output_unit,*)
      
      do ie=1,n_elements_nucleic_acid
        write(output_unit,*) atomic_number(ie), element_symbol(ie), atom_counts( atomic_number(ie) )
      end do

      ! Correct for phosphodiester bond formation ... net change is two hydrogen atoms eliminated / two oxygen atoms and one phosphorus atom gained for each phosphodiester formed
      ! This assumes that the DNA is made synthetically and has a Hydroxyl group (not phosphate) attached to the 5' end of the chain  

      n_phosphodiester_bonds = sequence_length - n_chains
      
      atom_counts( 1 )  = atom_counts( 1 )  - 2*n_phosphodiester_bonds ! Correction for Hydrogen,   atomic number 1
      atom_counts( 8 )  = atom_counts( 8 )  + 2*n_phosphodiester_bonds ! Correction for Oxygen,     atomic number 8
      atom_counts( 15 ) = atom_counts( 15 ) + 1*n_phosphodiester_bonds ! Correction for Phosphorus, atomic number 15
       
      write(output_unit,*) 
      write(output_unit,*) "Atom Counts, Corrected For Phosphodiester Bond Formation"
      write(output_unit,*) "Assuming ", n_chains, "nucleic acid chains are present in the sequence file"  
      write(output_unit,*) "And the 5' end of each chain carries a Hydroxyl group"  

      write(output_unit,*)
      
      do ie=1,n_elements_nucleic_acid
        write(output_unit,*) atomic_number(ie), element_symbol(ie), atom_counts( atomic_number(ie) )
      end do

      close(output_unit) ! close the output file

    end subroutine count_na_atoms

    pure function to_upper (str) Result (string)
    
    ! Simple conversion function courtesy of Seth Morton
    
    !   ==============================
    !   Changes a string to upper case
    !   ==============================

        Implicit None
        Character(*), Intent(In) :: str
        Character(LEN(str))      :: string

        Integer :: ic, i

        Character(26), Parameter :: upp = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    !   Capitalize each letter if it is lowercase
        string = str
        do i = 1, LEN_TRIM(str)
            ic = INDEX(low, str(i:i))
            if (ic > 0) string(i:i) = upp(ic:ic)
        end do

    end function to_upper
    
end module atom_count
