module rogers_analysis

  ! These are the core subroutines from Robert Blessing's program Rogers (written in Fortran 77) which was obtained from the author in the mid 90's
  ! Some minimal changes have been made to achieve a standardized treatment of dummy arguments (Fortran 90 standard), C/C++ interoperability (Fortran2003 standard) 
  ! and some elimination of features that were deleted in the Fortran2018 standard (principally in subroutine MATINV)
  ! The grid dimensions for the Patterson calculation, and a few constants are now also defined in the module declaration section immediately below, for clarity and ease of modification

  ! The original program did not apply the exact constraints which exist on the elements of the tensor B/U due to symmetry 
  ! Modifications have been made in the function SOLVE to apply the symmetry constraints to B/U, following least squares analysis
  ! The deviations from the symmetry constraints following least squares analysis are *very very* slight so this is a purely cosmetic issue. 
  
  
  ! The main Fortran 77 program has been replaced with a C++ calling function - defined inline in exp_manager.cc ...
  ! The central purpose of the calling function is to create and transform the intensity data (I -> I/ε), and expand it to a full hemisphere (space group P1).
  ! The c++ driver thereby assimilates the functions of former subroutines DECODE, DATA, READ1 and FCALC
  ! All subsequent calculations are handled by the Fortran subroutines from Robert Blessing, modified as above. 
  !  
  ! As per the original documentation the routines below are used to estimate 
  ! "THE ABSOLUTE SCALE FACTOR AND OVERALL ANISOTROPIC MEAN-SQUARE ATOMIC DISPLACEMENT PARAMETERS 
  !  BY MEANS OF AN ANALYSIS OF THE PATTERSON ORIGIN PEAK.  
  !  THE METHOD IS BASED ON AN IDEA DUE TO DONALD ROGERS (1965) 
  !  AS DEVELOPED AND IMPLEMENTED BY R. H. BLESSING AND D. A. LANGS, ACTA CRYST. A44, 729-735 (1988)."


  use,intrinsic :: iso_c_binding
  implicit none
  
  ! Any changes made to the grid dimensions *must* be mirrored in the c++ driving function
   
  ! number of grid points for the calculation of the Patterson along the axial directions a, b, c
  integer  (c_int), parameter, private       :: n_ax=200 
  
  ! Polar coordinate grid dimensions
  integer  (c_int), parameter, private       :: NR=25 ! radial grid spacing
  integer  (c_int), parameter, private       :: NA=12 ! angular grid spacing, this must be some multiple of 4. 
  integer  (c_int), parameter, private       :: n_pc = 1 + NR + (NR*NA*NA)/4 ! total number of points in the polar coordinate grid
  
  ! Some angular constants, used throughout
  
  real  (c_double), parameter, private      :: PI=ACOS(-1.0)
  real  (c_double), parameter, private      :: TWOPI=2*PI
  real  (c_double), parameter, private      :: RAD=PI/180
  real  (c_double), parameter, private      :: DEG=180/PI
    
  contains
!-----------------------------------------------------------------------
    SUBROUTINE FCALC (N,X,A1,B1,A2,B2,A3,B3,A4,B4,C0,FP,FPP,S,SUMFSQ) bind(c, name='FCALC')

       ! deprecated - this calculation is currently done in the main c++ calling routine.
             
      integer  (c_int), intent(in)                :: N
      integer  (c_int), intent(in),  dimension(N) :: X
      real  (c_double), intent(in),  dimension(N) :: A1
      real  (c_double), intent(in),  dimension(N) :: B1
      real  (c_double), intent(in),  dimension(N) :: A2
      real  (c_double), intent(in),  dimension(N) :: B2
      real  (c_double), intent(in),  dimension(N) :: A3
      real  (c_double), intent(in),  dimension(N) :: B3
      real  (c_double), intent(in),  dimension(N) :: A4
      real  (c_double), intent(in),  dimension(N) :: B4
      real  (c_double), intent(in),  dimension(N) :: C0
      real  (c_double), intent(in),  dimension(N) :: FP
      real  (c_double), intent(in),  dimension(N) :: FPP
      real  (c_double), intent(in)                :: S     !  = sin(theta)/lambda *not* 2sin(theta)/lambda
      real  (c_double), intent(out)               :: SUMFSQ
      
      integer  (c_int) :: i
            
      SUMFSQ=0
      DO I=1,N
        SUMFSQ=SUMFSQ+DBLE(X(I))*((A1(I)*EXP(-B1(I)*S**2) &
                                  +A2(I)*EXP(-B2(I)*S**2) &
                                  +A3(I)*EXP(-B3(I)*S**2) &
                                  +A4(I)*EXP(-B4(I)*S**2)+C0(I)+FP(I))**2 + FPP(I)**2)
      END DO
      RETURN   
    END SUBROUTINE FCALC
    
!-----------------------------------------------------------------------
    SUBROUTINE FTABLE (IZ,RM,A1,B1,A2,B2,A3,B3,A4,B4,C) bind(c, name='FTABLE')
        
! This version assumes provision of the Atomic Number Z, and returns relative atomic masses, and the coefficients required to compute the normal scattering of the associated element 

!
! NEUTRAL ATOM SCATTERING FACTOR COEFFICIENTS FROM CROMER & WEBER
! (INTERNATIONAL TABLES, VOL. IV, 1974)
!
! F(S) = A1*EXP(-B1*S**2) + A2*EXP(-B2*S**2) + A3*EXP(-B3*S**2) 
!          + A4*EXP(-B4*S**2) + C
!

      use,intrinsic :: iso_c_binding
      implicit none

      integer  (c_int), intent(in)  :: IZ
      real  (c_double), intent(out) :: RM
      real  (c_double), intent(out) :: A1
      real  (c_double), intent(out) :: B1
      real  (c_double), intent(out) :: A2
      real  (c_double), intent(out) :: B2
      real  (c_double), intent(out) :: A3
      real  (c_double), intent(out) :: B3
      real  (c_double), intent(out) :: A4
      real  (c_double), intent(out) :: B4
      real  (c_double), intent(out) :: C
      
      character*2,      dimension(98) :: SYMBOL
      real  (c_double), dimension(98) :: RMASS
      real  (c_double), dimension(98) :: AA1
      real  (c_double), dimension(98) :: BB1
      real  (c_double), dimension(98) :: AA2
      real  (c_double), dimension(98) :: BB2
      real  (c_double), dimension(98) :: AA3
      real  (c_double), dimension(98) :: BB3
      real  (c_double), dimension(98) :: AA4
      real  (c_double), dimension(98) :: BB4
      real  (c_double), dimension(98) :: CC
      real  (c_double), dimension(98) :: FPMO
      real  (c_double), dimension(98) :: FPPMO
      real  (c_double), dimension(98) :: FPCU
      real  (c_double), dimension(98) :: FPPCU
      real  (c_double), dimension(98) :: FPKISS
      
      integer  (c_int) :: I

!
! ELEMENT SYMBOLS
!
!rlk Element Symbols are currently unuused

      DATA SYMBOL / &
      'H ',                              'HE', &
      'LI','BE','B ','C ','N ','O ','F ','NE', &
      'NA','MG','AL','SI','P ','S ','CL','AR', & 
      'K ','CA', &
                'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN', &
                'GA','GE','AS','SE','BR','KR', &
      'RB','SR', &
                'Y ','ZR','NB','MO','TC','RU','RH','PD','AG','CD', &
                'IN','SN','SB','TE','I ','XE', &
      'CS','BA', &
                'LA', &
                     'CE','PR','ND','PM','SM','EU','GD', &
                       'TB','DY','HO','ER','TM','YB','LU', &
                     'HF','TA','W ','RE','OS','IR','PT','AU','HG', &
                'TL','PB','BI','PO','AT','RN', &
      'FR','RA', &
                'AC', &
                     'TH','PA','U ','NP','PU','AM','CM','BK','CF'/
!
! RELATIVE ATOMIC MASSES (NATURAL ISOTOPIC ABUNDANCES)
!
      DATA RMASS / 1.00797, 4.0026, 6.939, 9.0122, 10.811, 12.01115, &
       14.0067, 15.9994, 18.9984, 20.183, 22.9898, 24.312, 26.9815, &
       28.086, 30.9738, 32.064, 35.453, 39.948, 39.102, 40.08, 44.956, &
       47.90, 50.942, 51.996, 54.938, 55.847, 58.933, 58.71, 63.54, &
       65.37, 69.72, 72.59, 74.922, 78.96, 79.909, 83.80, 85.47, 87.62, &
       88.905, 91.22, 92.906, 95.94, 98, 101.07, 102.905, 106.4, &
       107.870, 112.40, 114.82, 118.69, 121.75, 127.60, 126.904, 131.30, &
       132.905, 137.34, 138.91, 140.12, 140.907, 144.24, 147, 150.35, &
       151.96, 157.25, 158.924, 162.50, 164.930, 167.26, 168.934, &
       173.04, 174.97, 178.49, 180.948, 183.85, 186.2, 190.2, 192.2, &
       195.09, 196.967, 200.59, 204.37, 207.19, 208.980, 210, 210, 222, &
       223, 226, 227, 232.038, 231, 238.03, 237, 242, 243, 247, 247, &
       249/
!
! SCATTERING FACTOR COEFFICIENTS
!
      DATA AA1 /   0.493002,  0.873400,  1.128200,  1.591900,  2.054500, &
        2.310000, 12.212590,  3.048500,  3.539200,  3.955300,  4.762600, &
        5.420400,  6.420200,  6.291500,  6.434500,  6.905300, 11.460390, &
        7.484500,  8.218600,  8.626600,  9.189000,  9.759500, 10.297100, &
       10.640600, 11.281900, 11.769490, 12.284090, 12.837590, 13.338000, &
       14.074290, 15.235400, 16.081600, 16.672300, 17.000590, 17.178890, &
       17.355490, 17.178400, 17.566290, 17.776000, 17.876490, 17.614190, &
        3.702500, 19.130090, 19.267390, 19.295700, 19.331890, 19.280800, &
       19.221400, 19.162390, 19.188900, 19.641790, 19.964400, 20.147200, &
       20.293300, 20.389200, 20.336100, 20.578000, 21.167090, 22.044000, &
       22.684490, 23.340490, 24.004190, 24.627390, 25.070900, 25.897590, &
       26.507000, 26.904900, 27.656290, 28.181900, 28.664090, 28.947600, &
       29.143990, 29.202390, 29.081800, 28.762100, 28.189400, 27.304900, &
       27.005900, 16.881890, 20.680890, 27.544600, 31.061700, 33.368890, &
       34.672600, 35.316290, 35.563090, 35.929900, 35.763000, 35.659690, &
       35.564490, 35.884700, 36.022790, 36.187390, 36.525400, 36.670590, &
       36.648800, 36.788100, 36.918500/
      DATA BB1 /  10.510890,  9.103700,  3.954600, 43.642700, 23.218500, &
       20.843900,  0.005700, 13.277090, 10.282500,  8.404200,  3.285000, &
        2.827500,  3.038700,  2.438600,  1.906700,  1.467900,  0.010400, &
        0.907200, 12.794890, 10.442090,  9.021300,  7.850800,  6.865700, &
        6.103800,  5.340900,  4.761100,  4.279100,  3.878500,  3.582800, &
        3.265500,  3.066900,  2.850900,  2.634500,  2.409800,  2.172300, &
        1.938400,  1.788800,  1.556400,  1.402900,  1.276180,  1.188650, &
        0.277200,  0.864132,  0.808520,  0.751536,  0.698655,  0.644600, &
        0.594600,  0.547600,  5.830300,  5.303400,  4.817420,  4.347000, &
        3.928200,  3.569000,  3.216000,  2.948170,  2.812190,  2.773930, &
        2.662480,  2.562700,  2.472740,  2.387900,  2.253410,  2.242560, &
        2.180200,  2.070510,  2.073560,  2.028590,  1.988900,  1.901820, &
        1.832620,  1.773330,  1.720290,  1.671910,  1.629030,  1.592790, &
        1.512930,  0.461100,  0.545000,  0.655150,  0.690200,  0.704000, &
        0.700999,  0.685870,  0.663100,  0.646453,  0.616341,  0.589092, &
        0.563359,  0.547751,  0.529300,  0.511929,  0.499384,  0.483629, &
        0.465154,  0.451018,  0.437533/
      DATA AA2 /   0.322912,  0.630900,  0.750800,  1.127800,  1.332600, &
        1.020000,  3.132200,  2.286800,  2.641200,  3.112500,  3.173600, &
        2.173500,  1.900200,  3.035300,  4.179100,  5.203400,  7.196400, &
        6.772300,  7.439800,  7.387300,  7.367900,  7.355800,  7.351100, &
        7.353700,  7.357300,  7.357300,  7.340900,  7.292000,  7.167600, &
        7.031800,  6.700600,  6.374700,  6.070100,  5.819600,  5.235800, &
        6.728600,  9.643500,  9.818400, 10.294590, 10.947990, 12.014390, &
       17.235590, 11.094790, 12.918190, 14.350090, 15.501700, 16.688500, &
       17.644390, 18.559600, 19.100490, 19.045500, 19.013790, 18.994900, &
       19.029800, 19.106200, 19.296990, 19.598990, 19.769500, 19.669690, &
       19.684690, 19.609490, 19.425790, 19.088590, 19.079800, 18.218500, &
       17.638300, 17.294000, 16.428490, 15.885100, 15.434490, 15.220800, &
       15.172590, 15.229290, 15.430000, 15.718890, 16.154990, 16.729590, &
       17.763900, 18.591290, 19.041700, 19.158400, 13.063690, 12.951000, &
       15.473290, 19.021100, 21.281600, 23.054700, 22.906400, 23.103190, &
       23.421900, 23.294790, 23.412790, 23.596400, 23.808300, 24.099190, &
       24.409600, 24.773600, 25.199490/
      DATA BB2 /  26.125700,  3.356800,  1.052400,  1.862300,  1.021000, &
       10.207500,  9.893300,  5.701100,  4.294400,  3.426200,  8.842200, &
       79.261090,  0.742600, 32.333690, 27.156990, 22.215100,  1.166200, &
       14.840700,  0.774800,  0.659900,  0.572900,  0.500000,  0.438500, &
        0.392000,  0.343200,  0.307200,  0.278400,  0.256500,  0.247000, &
        0.233300,  0.241200,  0.251600,  0.264700,  0.272600, 16.579600, &
       16.562300, 17.315090, 14.098790, 12.800600, 11.916000, 11.765990, &
        1.095800,  8.144870,  8.434670,  8.217580,  7.989290,  7.472600, &
        6.908900,  6.377600,  0.503100,  0.460700,  0.420885,  0.381400, &
        0.344000,  0.310700,  0.275600,  0.244475,  0.226836,  0.222087, &
        0.210628,  0.202088,  0.196451,  0.194200,  0.181951,  0.196143, &
        0.202172,  0.197940,  0.223545,  0.238849,  0.257119,  9.985190, &
        9.599900,  9.370460,  9.225900,  9.092270,  8.979480,  8.865530, &
        8.811740,  8.621600,  8.448400,  8.707510,  2.357600,  2.923800, &
        3.550780,  3.974580,  4.069100,  4.176190,  3.871350,  3.651550, &
        3.462040,  3.415190,  3.325300,  3.253960,  3.263710,  3.206470, &
        3.089970,  3.046190,  3.007750/
      DATA AA3 /   0.140191,  0.311200,  0.617500,  0.539100,  1.097900, &
        1.588600,  2.012500,  1.546300,  1.517000,  1.454600,  1.267400, &
        1.226900,  1.593600,  1.989100,  1.780000,  1.437900,  6.255600, &
        0.653900,  1.051900,  1.589900,  1.640900,  1.699100,  2.070300, &
        3.324000,  3.019300,  3.522200,  4.003400,  4.443800,  5.615800, &
        5.165200,  4.359100,  3.706800,  3.431300,  3.973100,  5.637700, &
        5.549300,  5.139900,  5.422000,  5.726290,  5.417320,  4.041830, &
       12.887590,  4.649010,  4.863370,  4.734250,  5.295370,  4.804500, &
        4.461000,  4.294800,  4.458500,  5.037100,  6.144870,  7.513800, &
        8.976700, 10.661990, 10.887990, 11.372690, 11.851300, 12.385600, &
       12.774000, 13.123490, 13.439590, 13.760290, 13.851790, 14.316690, &
       14.559590, 14.558300, 14.977890, 15.154190, 15.308690, 15.100000, &
       14.758600, 14.513500, 14.432700, 14.556400, 14.930500, 15.611490, &
       15.713100, 25.558190, 21.657500, 15.538000, 18.442000, 16.587700, &
       13.113800,  9.498870,  8.003700, 12.143890, 12.473890, 12.597700, &
       12.747300, 14.189100, 14.949090, 15.640190, 16.770700, 17.341500, &
       17.399000, 17.891900, 18.331690/
      DATA BB3 /   3.142360, 22.927590, 85.390500,103.483000, 60.349790, &
        0.568700, 28.997490,  0.323900,  0.261500,  0.230600,  0.313600, &
        0.380800, 31.547190,  0.678500,  0.526000,  0.253600, 18.519390, &
       43.898300,213.186900, 85.748390,136.108000, 35.633800, 26.893790, &
       20.262600, 17.867400, 15.353500, 13.535900, 12.176300, 11.396590, &
       10.316300, 10.780500, 11.446800, 12.947890, 15.237190,  0.260900, &
        0.226100,  0.274800,  0.166400,  0.125599,  0.117622,  0.204785, &
       11.003990, 21.570690, 24.799690, 25.874890, 25.205200, 24.660500, &
       24.700800, 25.849890, 26.890890, 27.907390, 28.528390, 27.766000, &
       26.465890, 24.387890, 20.207300, 18.772590, 17.608300, 16.766900, &
       15.885000, 15.100890, 14.399600, 13.754590, 12.933090, 12.664790, &
       12.189900, 11.440690, 11.360400, 10.997500, 10.664690,  0.261033, &
        0.275116,  0.295977,  0.321703,  0.350500,  0.382661,  0.417916, &
        0.424593,  1.482600,  1.572900,  1.963470,  8.618000,  8.793700, &
        9.556420, 11.382390, 14.042200, 23.105190, 19.988690, 18.598990, &
       17.830900, 16.923490, 16.092690, 15.362190, 14.945500, 14.313590, &
       13.434590, 12.894590, 12.404390/
      DATA AA4 /   0.040810,  0.178000,  0.465300,  0.702900,  0.706800, &
        0.865000,  1.166300,  0.867000,  1.024300,  1.125100,  1.112800, &
        2.307300,  1.964600,  1.541000,  1.490800,  1.586300,  1.645500, &
        1.644200,  0.865900,  1.021100,  1.468000,  1.902100,  2.057100, &
        1.492200,  2.244100,  2.304500,  2.348800,  2.380000,  1.673500, &
        2.410000,  2.962300,  3.683000,  4.277900,  4.354300,  3.985100, &
        3.537500,  1.529200,  2.669400,  3.265880,  3.657210,  3.533460, &
        3.742900,  2.712630,  1.567560,  1.289180,  0.605844,  1.046300, &
        1.602900,  2.039600,  2.466300,  2.682700,  2.523900,  2.273500, &
        1.990000,  1.495300,  2.695900,  3.287190,  3.330490,  2.824280, &
        2.851370,  2.875160,  2.896040,  2.922700,  3.545450,  2.953540, &
        2.965770,  3.638370,  2.982330,  2.987060,  2.989630,  3.716010, &
        4.300130,  4.764920,  5.119820,  5.441740,  5.675890,  5.833770, &
        5.783700,  5.860000,  5.967600,  5.525930,  5.969600,  6.469200, &
        7.025880,  7.425180,  7.443300,  2.112530,  3.210970,  4.086550, &
        4.807030,  4.172870,  4.188000,  4.185500,  3.479470,  3.493310, &
        4.216650,  4.232840,  4.243910/
      DATA BB4 /  57.799690,  0.982100,168.261000,  0.542000,  0.140300, &
       51.651190,  0.582600, 32.908900, 26.147590, 21.718390,129.423900, &
        7.193700, 85.088590, 81.693690, 68.164500, 56.171990, 47.778390, &
       33.392890, 41.684090,178.436900, 51.353100,116.104900,102.477900, &
       98.739890, 83.754300, 76.880490, 71.169200, 66.342100, 64.812600, &
       58.709700, 61.413490, 54.762490, 47.797190, 43.816290, 41.432800, &
       39.397200,164.934000,132.376000,104.354000, 87.662700, 69.795700, &
       61.658400, 86.847190, 94.292800, 98.606200, 76.898600, 99.815590, &
       87.482490, 92.802900, 83.957100, 75.282500, 70.840300, 66.877590, &
       64.265790,213.904000,167.201900,133.123900,127.113000,143.643900, &
      137.902900,132.720900,128.007000,123.173900,101.397900,115.361900, &
      111.873900, 92.656600,105.703000,102.960900,100.417000, 84.329800, &
       72.029000, 63.364390, 57.055990, 52.086100, 48.164700, 45.001090, &
       38.610300, 36.395590, 38.324600, 45.814890, 47.257900, 48.009290, &
       47.004500, 45.471490, 44.247290,150.645000,142.324900,117.020000, &
       99.172190,105.251000,100.613000, 97.490790,105.979900,102.272900, &
       88.483390, 86.003000, 83.788100/
      DATA CC /    0.003038,  0.006400,  0.037700,  0.038500, -0.193200, &
        0.215600,-11.529000,  0.250800,  0.277600,  0.351500,  0.676000, &
        0.858400,  1.115100,  1.140700,  1.114900,  0.866900, -9.557400, &
        1.444500,  1.422800,  1.375100,  1.332900,  1.280700,  1.219900, &
        1.183200,  1.089600,  1.036900,  1.011800,  1.034100,  1.191000, &
        1.304100,  1.718900,  2.131300,  2.531000,  2.840900,  2.955700, &
        2.825000,  3.487300,  2.506400,  1.912130,  2.069290,  3.755910, &
        4.387500,  5.404280,  5.378740,  5.328000,  5.265930,  5.179000, &
        5.069400,  4.939100,  4.782100,  4.590900,  4.352000,  4.071200, &
        3.711800,  3.335200,  2.773100,  2.146780,  1.862640,  2.058300, &
        1.984860,  2.028760,  2.209630,  2.574500,  2.419600,  3.583240, &
        4.297280,  4.567960,  5.920460,  6.756210,  7.566720,  7.976280, &
        8.581540,  9.243540,  9.887500, 10.472000, 11.000490, 11.472200, &
       11.688300, 12.065790, 12.608900, 13.174590, 13.411800, 13.578200, &
       13.677000, 13.710800, 13.690500, 13.724690, 13.621100, 13.526590, &
       13.431400, 13.428700, 13.396590, 13.357290, 13.381190, 13.359190, &
       13.288700, 13.275400, 13.267390/

! These anomalous dispersion corrections are not currently employed ...  

!
! ANOMALOUS DISPERSION CORRECTIONS FOR MO AND CU X-RAYS FROM CROMER  
! LIBERMAN (INTERNATIONAL TABLES, VOL. IV, 1974)
!
      DATA FPMO /  &
        0.000,  0.000,  0.000,  0.000,  0.000,  0.002,  0.004,  0.008, &
        0.014,  0.021,  0.030,  0.042,  0.056,  0.072,  0.090,  0.110, &
        0.132,  0.155,  0.179,  0.203,  0.226,  0.248,  0.267,  0.284, &
        0.295,  0.301,  0.299,  0.285,  0.263,  0.222,  0.163,  0.081, &
       -0.030, -0.178, -0.374, -0.652, -1.044, -1.657, -2.951, -2.965, &
       -2.197, -1.825, -1.590, -1.420, -1.287, -1.177, -1.085, -1.005, &
       -0.936, -0.873, -0.816, -0.772, -0.726, -0.684, -0.644, -0.613, &
       -0.588, -0.564, -0.530, -0.535, -0.530, -0.533, -0.542, -0.564, &
       -0.591, -0.619, -0.666, -0.723, -0.795, -0.884, -0.988, -1.118, &
       -1.258, -1.421, -1.598, -1.816, -2.066, -2.352, -2.688, -3.084, &
       -3.556, -4.133, -4.861, -5.924, -7.444, -8.862, -7.912, -7.620, &
       -7.725, -8.127, -8.960,-10.673,-11.158, -9.725, -8.926, -8.416, &
       -7.990, -7.683/
      DATA FPPMO / &
        0.000,  0.000,  0.000,  0.000,  0.001,  0.002,  0.003,  0.006, &
        0.010,  0.016,  0.025,  0.036,  0.052,  0.071,  0.095,  0.124, &
        0.159,  0.201,  0.250,  0.306,  0.372,  0.446,  0.530,  0.624, &
        0.729,  0.845,  0.973,  1.113,  1.266,  1.431,  1.609,  1.801, &
        2.007,  2.223,  2.456,  2.713,  2.973,  3.264,  3.542,  0.560, &
        0.621,  0.688,  0.759,  0.836,  0.919,  1.007,  1.101,  1.202, &
        1.310,  1.424,  1.546,  1.675,  1.812,  1.958,  2.119,  2.282, &
        2.452,  2.632,  2.845,  3.018,  3.225,  3.442,  3.669,  3.904, &
        4.151,  4.410,  4.678,  4.958,  5.248,  5.548,  5.858,  6.185, &
        6.523,  6.872,  7.232,  7.605,  7.990,  8.388,  8.798,  9.223, &
        9.659, 10.102, 10.559, 11.042,  9.961, 10.403,  7.754,  8.105, &
        8.472,  8.870,  9.284,  9.654,  4.148,  4.330,  4.511,  4.697, &
        4.908,  5.107/
      DATA FPCU / &
        0.000,  0.000,  0.001,  0.003,  0.008,  0.017,  0.029,  0.047, &
        0.069,  0.097,  0.129,  0.165,  0.204,  0.244,  0.283,  0.319, &
        0.348,  0.366,  0.365,  0.341,  0.285,  0.189,  0.035, -0.198, &
       -0.568, -1.179, -2.464, -2.956, -2.019, -1.612, -1.354, -1.163, &
       -1.011, -0.879, -0.767, -0.665, -0.574, -0.465, -0.386, -0.314, &
       -0.248, -0.191, -0.145, -0.105, -0.077, -0.059, -0.060, -0.079, &
       -0.126, -0.194, -0.287, -0.418, -0.579, -0.783, -1.022, -1.334, &
       -1.716, -2.170, -2.939, -3.431, -4.357, -5.696, -7.718, -9.242, &
       -9.498,-10.423,-12.255, -9.733, -8.488, -7.701, -7.713, -6.715, &
       -6.351, -6.048, -5.790, -5.581, -5.391, -5.233, -5.096, -4.990, &
       -4.883, -4.818, -4.776, -4.756, -4.772, -4.787, -4.833, -4.898, &
       -4.994, -5.091, -5.216, -5.359, -5.529, -5.712, -5.930, -6.176, &
       -6.498, -6.798/
      DATA FPPCU / &
        0.000,  0.000,  0.000,  0.001,  0.004,  0.009,  0.018,  0.032, &
        0.053,  0.083,  0.124,  0.177,  0.246,  0.330,  0.434,  0.557, &
        0.702,  0.872,  1.066,  1.286,  1.533,  1.807,  2.110,  2.443, &
        2.808,  3.204,  3.608,  0.509,  0.589,  0.678,  0.777,  0.886, &
        1.006,  1.139,  1.283,  1.439,  1.608,  1.820,  2.025,  2.245, &
        2.482,  2.735,  3.005,  3.296,  3.605,  3.934,  4.282,  4.653, &
        5.045,  5.459,  5.894,  6.352,  6.835,  7.348,  7.904,  8.460, &
        9.036,  9.648, 10.535, 10.933, 11.614, 12.320, 11.276, 11.946, &
        9.242,  9.748,  3.704,  3.937,  4.181,  4.432,  4.693,  4.977, &
        5.271,  5.577,  5.891,  6.221,  6.566,  6.925,  7.297,  7.686, &
        8.089,  8.505,  8.930,  9.383,  9.843, 10.317, 10.803, 11.296, &
       11.799, 12.330, 12.868, 13.409, 13.967, 14.536, 15.087, 15.634, &
       16.317, 16.930/
!
! WAVELENGTH-INDEPENDENT KISSEL AND PRATT CORRECTIONS TO THE CROMER AND
! LIBERMAN FP VALUES [L. KISSEL AND R.H. PRATT (1990).  ACTA CRYST. A46,
! 170-175].
!
      DATA FPKISS / &
        0.000,  0.000,  0.001,  0.000,  0.001,  0.001,  0.002,  0.003, &
        0.004,  0.004,  0.006,  0.008,  0.008,  0.011,  0.012,  0.014, &
        0.017,  0.020,  0.022,  0.025,  0.028,  0.031,  0.035,  0.039, &
        0.042,  0.048,  0.052,  0.057,  0.061,  0.067,  0.073,  0.079, &
        0.085,  0.092,  0.099,  0.106,  0.114,  0.122,  0.130,  0.138, &
        0.147,  0.156,  0.166,  0.175,  0.186,  0.196,  0.207,  0.219, &
        0.230,  0.242,  0.255,  0.267,  0.281,  0.294,  0.308,  0.323, &
        0.338,  0.354,  0.369,  0.386,  0.402,  0.419,  0.437,  0.455, &
        0.474,  0.493,  0.512,  0.532,  0.553,  0.574,  0.596,  0.617, &
        0.640,  0.663,  0.687,  0.711,  0.736,  0.762,  0.788,  0.814, &
        0.842,  0.870,  0.899,  0.928,  0.957,  0.988,  1.018,  1.050, &
        1.083,  1.115,  1.149,  1.184,  1.219,  1.255,  1.292,  1.330, &
        1.368,  1.407/

      RM=RMASS(IZ)
      A1=AA1(IZ)
      B1=BB1(IZ)
      A2=AA2(IZ)
      B2=BB2(IZ)
      A3=AA3(IZ)
      B3=BB3(IZ)
      A4=AA4(IZ)
      B4=BB4(IZ)
      C=CC(IZ)
!rlk      IF (XRAY.EQ.'MO') THEN
!rlk        FP=FPMO(I)+FPKISS(I)
!rlk        FPP=FPPMO(I)
!rlk      ELSE IF (XRAY.EQ.'CU') THEN
!rlk        FP=FPCU(I)+FPKISS(I)
!rlk        FPP=FPPCU(I)
!rlk      ELSE
!rlk        FP=0
!rlk        FPP=0
!rlk      END IF
      RETURN
    END SUBROUTINE FTABLE
!-----------------------------------------------------------------------
    SUBROUTINE GRID (R,GINV,U,V,W) bind(c, name='GRID')

!
! CALCULATE FRACTIONAL CRYSTAL COORDINATES OF A HEMISPHERICAL POLAR
! COORDINATE GRID.
!
!     W|   R        U = R SIN(TH) COS(PH)
!      |   /|       V = R SIN(TH) SIN(PH)
!      |TH/ |       W = R COS(TH)
!      | /  |
!      |/   |
!      _____|__V    
!     /\    |  
!    /  \   |
!   / PH \  |
! U/      \ |
!          \|
!
!
! A GRID WITH 25 POINTS ALONG RADII AT 30 DEGREE INTERVALS OF THETA AND
! PHI GIVES 1 + 25 + 25*12*3 = 926 GRID POINTS.
!
 
      use,intrinsic :: iso_c_binding
      implicit none
 
      real  (c_double), intent(inout)                  :: R
      real  (c_double), intent(in),  dimension(3,3)    :: GINV
      real  (c_double), intent(out),  dimension(n_pc)  :: U
      real  (c_double), intent(out),  dimension(n_pc)  :: V
      real  (c_double), intent(out),  dimension(n_pc)  :: W

      real  (c_double),   dimension(3,3) :: G
      real  (c_double),   dimension(3,3) :: T
      real  (c_double),   dimension(3)   :: X
      real  (c_double),   dimension(3)   :: Y
      real  (c_double)                   :: PH
      real  (c_double)                   :: SINPH
      real  (c_double)                   :: COSPH
      real  (c_double)                   :: TH
      real  (c_double)                   :: SINTH
      real  (c_double)                   :: COSTH
      real  (c_double)                   :: DA
      real  (c_double)                   :: DR
      real  (c_double)                   :: DETG
      real  (c_double)                   :: DETGI
      integer  (c_int) :: I
      integer  (c_int) :: J
      integer  (c_int) :: K

      integer  (c_int) :: N     
            
      DR=R/NR
      DA=2*PI/NA
      N=0
      DO I=0,NR
        N=N+1
        U(N)=0
        V(N)=0
        W(N)=I*DR
      END DO
      DO I=0,NA-1
        PH=I*DA
        SINPH=SIN(PH)
        COSPH=COS(PH)
      DO J=1,NA/4
        TH=J*DA
        SINTH=SIN(TH)
        COSTH=COS(TH)
      DO K=1,NR
        R=K*DR
        N=N+1
        U(N)=R*SINTH*COSPH
        V(N)=R*SINTH*SINPH
        W(N)=R*COSTH
      END DO
      END DO
      END DO
!
! TRANSFORM TO FRACTIONAL CRYSTAL COORDINATES FROM ORTHONORMAL
! COORDINATES WITH X ALONG A, Y PERPENDICULAR TO A IN THE AB PLANE, AND
! Z PERPENDICULAR TO THE AB PLANE.
!
      CALL MATINV3 (GINV,G,DETGI)
      DETG=1/DETGI
      T(1,1)=1/SQRT(G(1,1))
      T(1,2)=-G(1,2)/SQRT(DETG*G(1,1)*GINV(3,3))
      T(1,3)=GINV(3,1)/SQRT(GINV(3,3))
      T(2,1)=0
      T(2,2)=SQRT(G(1,1)/(DETG*GINV(3,3)))
      T(2,3)=GINV(3,2)/SQRT(GINV(3,3))
      T(3,1)=0
      T(3,2)=0
      T(3,3)=SQRT(GINV(3,3))
      DO I=1,N
        X(1)=U(I)
        X(2)=V(I)
        X(3)=W(I)
        DO J=1,3
          Y(J)=0
        DO K=1,3
          Y(J)=Y(J)+T(J,K)*X(K)
        END DO
        END DO
        U(I)=Y(1)
        V(I)=Y(2)
        W(I)=Y(3)
      END DO
      RETURN
    END SUBROUTINE GRID
!-----------------------------------------------------------------------
    SUBROUTINE MATINV (N,A,D) bind(c, name='MATINV')
!
!         PURPOSE
!           INVERT A SYMMETRIC MATRIX AND CALCULATE ITS DETERMINANT.
!
!         USAGE
!           CALL MATINV (N,A,D,DLOG)
!
!         DESCRIPTION OF PARAMETERS
!           A    - INPUT MATRIX, WHICH IS REPLACED BY ITS INVERSE
!           N    - DEGREE OF MATRIX (ORDER OF DETERMINANT)
!           D    - DETERMINANT OF EITHER THE INPUT MATRIX OR THE SCALED
!                  INPUT MATRIX
!           DLOG - THE NATURAL LOGARITHM OF THE DETERMINANT OF THE INPUT
!                  MATRIX
!
!           DLOG IS RETURNED IN CASE THE MAGNITUDE OF D = EXP(DLOG) IS
!           TOO LARGE.  THE CALLING PROGRAM SHOULD COMPARE LOG(D) AND
!           DLOG TO DETERMINE WHICH VALUE OF D WAS RETURNED.
!
!         SUBROUTINES AND FUNCTION SUBPROGRAM REQUIRED
!           NONE
!

      use,intrinsic :: iso_c_binding
      implicit none
  
      integer  (c_int), intent(in)                      :: N
      real  (c_double), intent(inout),   dimension(N,N) :: A
      real  (c_double), intent(out)                     :: D

      real  (c_double), dimension(N) :: DIAG
      integer  (c_int), dimension(N) :: IK
      integer  (c_int), dimension(N) :: JK

      
      integer  (c_int) :: I
      integer  (c_int) :: J
      integer  (c_int) :: K
      integer  (c_int) :: L
      real  (c_double) :: AMAX
      real  (c_double) :: DLOG
      real  (c_double) :: QMAX
      real  (c_double) :: QMIN
      real  (c_double) :: T

!
!         SAVE SQUARE-ROOT DIAGONAL ELEMENTS OF MATRIX.
!
      DO I=1,N
 10     DIAG(I)=SQRT(ABS(A(I,I)))
      END DO
!
!         SCALE MATRIX TO AVOID ARITHMETIC OVERFLOW OR UNDERFLOW DURING
!         INVERSION.
!
      DO I=1,N
        DO J=I,N
          A(I,J)=A(I,J)/(DIAG(I)*DIAG(J))
 11       A(J,I)=A(I,J)
        END DO
      END DO
      
!      
!         INVERT SCALED MATRIX AND CALCULATE ITS DETERMINANT.
!
      D=1
      DO K=1,N
!
!         FIND LARGEST ELEMENT A(I,J) IN REST OF MATRIX.
!
        AMAX=0
 21     DO I=K,N
          DO J=K,N
            IF (ABS(AMAX).GT.ABS(A(I,J))) GO TO 30
            AMAX=A(I,J)
            IK(K)=I
            JK(K)=J
 30         CONTINUE
          END DO
        END DO
!
!         INTERCHANGE ROWS AND COLUMNS TO PUT AMAX IN A(K,K).
!
        IF (AMAX.NE.0) GO TO 41
        D=0
        RETURN
 41     I=IK(K)
 
!        IF (I-K) 21,51,43
!       Replacement for the arithmetic if statement - a deleted feature in Fortran 2018
        IF      ((I-K) .lt. 0) THEN
          GOTO 21
        ELSE IF ((I-K) .eq. 0) THEN
          GOTO 51
        ELSE IF ((I-K) .gt. 0) THEN
          GOTO 43
        END IF
           
 43     DO J=1,N
          T=A(K,J)
          A(K,J)=A(I,J)
 50       A(I,J)=-T
        END DO
 51     J=JK(K)
!        IF (J-K) 21,61,53
!       Replacement for the arithmetic if statement - a deleted feature in Fortran 2018
        IF      ((J-K) .lt. 0) THEN
          GOTO 21
        ELSE IF ((J-K) .eq. 0) THEN
          GOTO 61
        ELSE IF ((J-K) .gt. 0) THEN
          GOTO 53
        END IF

 53     DO I=1,N
          T=A(I,K)
          A(I,K)=A(I,J)
 60       A(I,J)=-T
        END DO
!
!         ACCUMULATE ELEMENTS OF INVERSE MATRIX.
!
 61     DO I=1,N
          IF (I.EQ.K) GO TO 70
          A(I,K)=-A(I,K)/AMAX
 70       CONTINUE
        END DO
        DO I=1,N
          DO J=1,N
            IF (I.EQ.K) GO TO 80
            IF (J.EQ.K) GO TO 80
            A(I,J)=A(I,J)+A(I,K)*A(K,J)
 80         CONTINUE
          END DO
        END DO
        DO J=1,N
          IF (J.EQ.K) GO TO 90
          A(K,J)=A(K,J)/AMAX
 90       CONTINUE
        END DO
        A(K,K)=1/AMAX
 100    D=D*AMAX
      END DO
!
!         RESTORE ORDERING OF MATRIX.
!
      DO L=1,N
        K=N-L+1
        J=IK(K)
        IF (J.LE.K) GO TO 111
        DO I=1,N
          T=A(I,K)
          A(I,K)=-A(I,J)
 110      A(I,J)=T
        END DO
 111    I=JK(K)
        IF (I.LE.K) GO TO 130
        DO J=1,N
          T=A(K,J)
          A(K,J)=-A(I,J)
 120      A(I,J)=T
        END DO
 130    CONTINUE
      END DO
!
!         SCALE INVERSE MATRIX.
!
      DO I=1,N
        DO J=I,N
          A(I,J)=A(I,J)/(DIAG(I)*DIAG(J))
 140      A(J,I)=A(I,J)
        END DO
      END DO
!
!         SCALE DETERMINANT (IF POSSIBLE).
!
      DLOG=0
      DO I=1,N
 150    DLOG=DLOG+LOG(DIAG(I))
      END DO
      DLOG=LOG(ABS(D))+2*DLOG
!
!         TEST AGAINST RANGE OF MACHINE-ALLOWED MAGNITUDES.
!
!         (QMIN = 2**(-128) = 0.29E-38 AND QMAX = 2**(+128)/2 = 1.7E+38
!         FOR VAX COMPUTERS.)
!
      QMIN=1E-37
      QMAX=1E+37
      IF (LOG(QMIN).LT.ABS(DLOG).AND.ABS(DLOG).LT.LOG(QMAX)) D=(D/ABS(D))*EXP(DLOG)
      RETURN
    END SUBROUTINE MATINV

!-----------------------------------------------------------------------
    SUBROUTINE MATINV3 (A,B,D) bind(c, name='MATINV3')
      
          !
          ! FOR A SYMMETRIC 3 X 3 MATRIX A, B = A**(-1) AND D = DET(A).
          !
          
          use,intrinsic :: iso_c_binding
          implicit none
            
          real  (c_double), intent(in),  dimension(3,3) :: A
          real  (c_double), intent(out), dimension(3,3) :: B
          real  (c_double), intent(out)                 :: D
 
          integer  (c_int)                          :: I
          integer  (c_int)                          :: J
          integer  (c_int)                          :: K
          integer  (c_int)                          :: L
          integer  (c_int)                          :: M
          integer  (c_int)                          :: N
            
          DO I=1,3
            K=1+MOD(I+1,3)
            M=1+MOD(I+3,3)
          DO J=1,3
            L=1+MOD(J+1,3)
            N=1+MOD(J+3,3)
            B(I,J)=A(L,K)*A(N,M)-A(L,M)*A(N,K)
          END DO
          END DO
          D=0
          DO I=1,3
            D=D+A(I,1)*B(1,I)
          END DO
          DO I=1,3
          DO J=1,3
            B(I,J)=B(I,J)/D
          END DO
          END DO
          RETURN
    END SUBROUTINE MATINV3
!-----------------------------------------------------------------------
    SUBROUTINE METRIC (A,V,GA,GB) bind(c, name='METRIC')
!
! A  = CRYSTAL LATTICE PARAMETERS
! B  = RECIPROCAL LATTICE PARAMETERS
! GA = CRYSTAL SPACE METRIC MATRIX
! GB = RECIPROCAL SCACE METRIC MATRIX
!

      use,intrinsic :: iso_c_binding
      implicit none

      real  (c_double), intent(in),  dimension(6)   :: A
      real  (c_double), intent(out)                 :: V
      real  (c_double), intent(out), dimension(3,3) :: GA
      real  (c_double), intent(out), dimension(3,3) :: GB
      

      real  (c_double) :: CA
      real  (c_double) :: CB
      real  (c_double) :: CG
      real  (c_double) :: SA
      real  (c_double) :: SB
      real  (c_double) :: SG
      real  (c_double), dimension(6) :: B
      
      CA=COS(A(4)*RAD)
      CB=COS(A(5)*RAD)
      CG=COS(A(6)*RAD)
      SA=SIN(A(4)*RAD)
      SB=SIN(A(5)*RAD)
      SG=SIN(A(6)*RAD)
      V=A(1)*A(2)*A(3)*SQRT(1-CA**2-CB**2-CG**2+2*CA*CB*CG)
      B(1)=A(2)*A(3)*SA/V
      B(2)=A(1)*A(3)*SB/V
      B(3)=A(1)*A(2)*SG/V
      B(4)=(CB*CG-CA)/(SB*SG)
      B(5)=(CA*CG-CB)/(SA*SG)
      B(6)=(CA*CB-CG)/(SA*SB)
      GA(1,1)=A(1)*A(1)
      GA(1,2)=A(1)*A(2)*CG
      GA(1,3)=A(1)*A(3)*CB
      GA(2,1)=GA(1,2)
      GA(2,2)=A(2)*A(2)
      GA(2,3)=A(2)*A(3)*CA
      GA(3,1)=GA(1,3)
      GA(3,2)=GA(2,3)
      GA(3,3)=A(3)*A(3)
      GB(1,1)=B(1)*B(1)
      GB(1,2)=B(1)*B(2)*B(6)
      GB(1,3)=B(1)*B(3)*B(5)
      GB(2,1)=GB(1,2)
      GB(2,2)=B(2)*B(2)
      GB(2,3)=B(2)*B(3)*B(4)
      GB(3,1)=GB(1,3)
      GB(3,2)=GB(2,3)
      GB(3,3)=B(3)*B(3)
      B(4)=ACOS(B(4))*DEG
      B(5)=ACOS(B(5))*DEG
      B(6)=ACOS(B(6))*DEG
      RETURN
    END SUBROUTINE METRIC
!-----------------------------------------------------------------------
    SUBROUTINE PADJ (P,PMIN,P0) bind(c, name='PADJ')
    
!
! SCULPT A MINIMUM PROFILE ORIGIN PEAK.
!

      use,intrinsic :: iso_c_binding
      implicit none

      real  (c_double), intent(inout),    dimension(n_pc) :: P
      real  (c_double), intent(inout)                     :: PMIN
      real  (c_double), intent(out)                       :: P0

      integer  (c_int) :: I
      integer  (c_int) :: J
      integer  (c_int) :: K
      real  (c_double) :: Q
      
      integer  (c_int) :: N
            
      N=1
      Q=P(1)
      DO I=1,NR
        N=N+1
        P(N)=MIN(P(N),Q)
        Q=P(N)
      END DO
      DO I=0,NA-1
      DO J=1,NA/4
        Q=P(1)
      DO K=1,NR
        N=N+1
        P(N)=MIN(P(N),Q)
        Q=P(N)
      END DO
      END DO
      END DO
      DO I=1,N
        PMIN=MIN(PMIN,P(I))
      END DO
      P0=P(1)-PMIN
      RETURN
    END SUBROUTINE PADJ
!-----------------------------------------------------------------------
    SUBROUTINE PFIT (U,V,W,P,PMIN,P0,PP,VCELL,F000,SUMZSQ,F,NCYCLE, R) bind(c, name='PFIT')

      !
      ! REFINE THE FIT OF THE TRIVARIATE GAUSSIAN DENSITY FUNCTION TO THE
      ! PATTERSON ORIGIN PEAK.
      !
      !   P(U,V,W) - PMIN = P0*EXP(-(P11*U**2 + P22*V**2 + P33*W**2
      !                              + 2*P12*U*V + 2*P13*U*W + 2*P23*V*W))
      !
      ! WHERE
      !
      !   PMIN = -F000**2/(SCALEK**2*VCELL*SUM(J) Z(J)**2)
      !
      ! AND
      !
      !   SCALEK**2 = SQRT(DET(PIJ))/(SQRT(PI**3)*VCELL*P0)
      !

      use,intrinsic :: iso_c_binding
      implicit none

      real  (c_double), intent(in),   dimension(n_pc)   :: U
      real  (c_double), intent(in),   dimension(n_pc)   :: V
      real  (c_double), intent(in),   dimension(n_pc)   :: W
      real  (c_double), intent(in),   dimension(n_pc)   :: P
      real  (c_double), intent(inout)                  :: PMIN
      real  (c_double), intent(inout)                  :: P0
      real  (c_double), intent(inout),  dimension(3,3) :: PP
      real  (c_double), intent(in)                   :: VCELL
      real  (c_double), intent(in)                   :: F000
      real  (c_double), intent(in)                   :: SUMZSQ
      real  (c_double), intent(in)                   :: F
      integer  (c_int), intent(out)                  :: NCYCLE
      real  (c_double), intent(out)                  :: R
      
      real  (c_double), dimension(n_pc) :: WT
      real  (c_double), dimension(8,8)  :: AA
      real  (c_double), dimension(8)    :: BB
      real  (c_double), dimension(8)    :: XX
      real  (c_double), dimension(8,7)  :: QQ
      real  (c_double), dimension(7,8)  :: QQT
      real  (c_double), dimension(7,7)  :: CC
      real  (c_double), dimension(7)    :: DD
      real  (c_double), dimension(7)    :: UU
      integer  (c_int) :: I
      integer  (c_int) :: J
      integer  (c_int) :: K
      integer  (c_int) :: L
      integer  (c_int) :: NOBS
      integer  (c_int) :: NPAR
      integer  (c_int) :: NCON
      real  (c_double) :: CHISQ
      real  (c_double) :: CONST1
      real  (c_double) :: CONST2
      real  (c_double) :: SUMSQ
      real  (c_double) :: D1
      real  (c_double) :: D2
      real  (c_double) :: D3
      real  (c_double) :: D4
      real  (c_double) :: D5
      real  (c_double) :: D6
      real  (c_double) :: DET
      real  (c_double) :: DETP
      real  (c_double) :: R1
      real  (c_double) :: SCLSQ
      real  (c_double) :: SUMWT
      real  (c_double) :: T
      real  (c_double) :: UI
      real  (c_double) :: VI
      real  (c_double) :: WI
      real  (c_double) :: YI
      real  (c_double) :: Z
      real  (c_double) :: Z1
            
      integer  (c_int) :: N
      N = n_pc
      
      NOBS=N
      NPAR=8
      NCON=1
      CONST1=VCELL*SQRT(PI**3)
      CONST2=VCELL*SUMZSQ
      !
      ! CALCULATE INITIAL GOODNESS-OF-FIT STATISTICS.
      !
      CHISQ=0
      SUMSQ=0
      DO I=1,N
        CHISQ=CHISQ+(P(I)-PMIN-PUVW(U(I),V(I),W(I),P0,PP))**2
        SUMSQ=SUMSQ+P(I)**2
      END DO
      Z=SQRT(CHISQ/(NOBS-NPAR))
      R=SQRT(CHISQ/SUMSQ)
      !
      ! ITERATE THROUGH LEAST-SQUARES LOOP.
      !
      NCYCLE=0
    1 NCYCLE=NCYCLE+1
      !
      ! EVALUATE THE DERIVATIVES, AND BUILD THE NORMAL MATRIX AND VECTOR.
      !
      NOBS=0
      DO J=1,NPAR
        DO K=1,NPAR
          AA(J,K)=0
        END DO
        BB(J)=0
      END DO
      DO I=1,N
        CALL WEIGHT (P(I)-PMIN,P0,PP,U(I),V(I),W(I),Z,WT(I))
        IF (WT(I).NE.0) THEN
          UI=U(I)
          VI=V(I)
          WI=W(I)
          YI=PUVW(UI,VI,WI,P0,PP)
          XX(1)=-UI**2*YI
          XX(2)=-VI**2*YI
          XX(3)=-WI**2*YI
          XX(4)=-2*UI*VI*YI
          XX(5)=-2*UI*WI*YI
          XX(6)=-2*VI*WI*YI
          XX(7)=YI/P0
          XX(8)=1
          YI=P(I)-PMIN-YI
          NOBS=NOBS+1
          DO J=1,NPAR
            DO K=J,NPAR
              AA(J,K)=AA(J,K)+WT(I)*XX(J)*XX(K)
              AA(K,J)=AA(J,K)
            END DO
            BB(J)=BB(J)+WT(I)*XX(J)*YI
          END DO
        END IF
      END DO
      !
      ! BUILD PMIN CONSTRAINT MATRIX QQ(8,7).
      !
      ! CONSTRAINT MATRIX ELEMENT (I,J) IS THE DERIVATIVE OF UNCONSTRAINED
      ! PARAMETER (I) WITH RESPECT TO CONSTRAINED PARAMETER (J).
      !
      DO I=1,7
        DO J=1,7
          QQ(I,J)=0
        END DO
        QQ(I,I)=1
      END DO
      DETP=  PP(1,1)*PP(2,2)*PP(3,3)+PP(1,2)*PP(2,3)*PP(3,1) &
            +PP(1,3)*PP(2,1)*PP(3,2)-PP(3,1)*PP(2,2)*PP(1,3) &
            -PP(3,2)*PP(2,3)*PP(1,1)-PP(3,3)*PP(2,1)*PP(1,2)
      D1=(PP(2,2)*PP(3,3)-PP(2,3)**2)/(2*DETP)
      D2=(PP(1,1)*PP(3,3)-PP(1,3)**2)/(2*DETP)
      D3=(PP(1,1)*PP(2,2)-PP(1,2)**2)/(2*DETP)
      D4=(2*PP(1,3)*PP(2,3)-PP(1,2)*PP(3,3))/(2*DETP)
      D5=(2*PP(1,2)*PP(3,1)-PP(1,3)*PP(2,2))/(2*DETP)
      D6=(2*PP(2,1)*PP(3,1)-PP(2,3)*PP(1,1))/(2*DETP)
      T=PMIN
      SCLSQ=SQRT(DETP)/(CONST1*P0)
      PMIN=-F000**2/(SCLSQ*CONST2)
      QQ(8,1)=-0.5*PMIN*D1/DETP
      QQ(8,2)=-0.5*PMIN*D2/DETP
      QQ(8,3)=-0.5*PMIN*D3/DETP
      QQ(8,4)=-0.5*PMIN*D4/DETP
      QQ(8,5)=-0.5*PMIN*D5/DETP
      QQ(8,6)=-0.5*PMIN*D6/DETP
      QQ(8,7)=PMIN/P0
      PMIN=T
      !
      ! APPLY CONSTRAINTS TO THE NORMAL MATRIX AND VECTOR.
      !
      DO I=1,NPAR-NCON
        DO J=1,NPAR
          QQT(I,J)=QQ(J,I)
        END DO
      END DO
      DO I=1,NPAR-NCON
        DO J=1,NPAR-NCON
          CC(I,J)=0
          DO K=1,NPAR
            DO L=1,NPAR
              CC(I,J)=CC(I,J)+QQT(I,K)*AA(K,L)*QQ(L,J)
            END DO
          END DO
        END DO
      END DO
      DO I=1,NPAR-NCON
        DD(I)=0
        DO J=1,NPAR
          DD(I)=DD(I)+QQT(I,J)*BB(J)
        END DO
      END DO
      !
      ! INVERT THE CONSTRAINED NORMAL MATRIX.
      !
      CALL MATINV (NPAR-NCON,CC,DET)
      IF (DET.EQ.0) STOP 'SINGULAR NORMAL MATRIX'
      !
      ! CALCULATE THE CONSTRAINED PARAMETER SHIFTS.
      !
      DO I=1,NPAR-NCON
        UU(I)=0
        DO J=1,NPAR-NCON
          UU(I)=UU(I)+CC(I,J)*DD(J)
        END DO
      END DO
      !
      ! CALCULATE AND APPLY THE UNCONSTRAINED PARAMETER SHIFTS.
      !
      DO I=1,NPAR
        XX(I)=0
        DO J=1,NPAR-NCON
          XX(I)=XX(I)+QQ(I,J)*UU(J)
        END DO
      END DO
      PP(1,1)=PP(1,1)+F*XX(1)
      PP(2,2)=PP(2,2)+F*XX(2)
      PP(3,3)=PP(3,3)+F*XX(3)
      PP(1,2)=PP(1,2)+F*XX(4)
      PP(1,3)=PP(1,3)+F*XX(5)
      PP(2,3)=PP(2,3)+F*XX(6)
      PP(2,1)=PP(1,2)
      PP(3,1)=PP(1,3)
      PP(3,2)=PP(2,3)
      P0=P0+F*XX(7)
      PMIN=PMIN+F*XX(8)
      !
      ! CALCULATE GOODNESS-OF-FIT STATISTICS.
      !
      CHISQ=0
      SUMSQ=0
      SUMWT=0
      DO I=1,N
        CHISQ=CHISQ+WT(I)*(P(I)-PMIN-PUVW(U(I),V(I),W(I),P0,PP))**2
        SUMSQ=SUMSQ+WT(I)*P(I)**2
        SUMWT=SUMWT+WT(I)
      END DO
      Z1=Z
      R1=R
      Z=SQRT((CHISQ/SUMWT)*NOBS/(NOBS-(NPAR-NCON)))
      R=SQRT(CHISQ/SUMSQ)
      !
      ! TEST FOR CONVERGENCE.
      !
      T=(R-R1)/R1
      IF (T.GT.0) GO TO 9
      IF (T.GT.-1E-6) RETURN
      IF (NCYCLE.LE.50) GO TO 1
    9 CONTINUE
      !
      ! DIVERGENT REFINEMENT.  BACK-SHIFT THE PARAMETERS.
      !
      NCYCLE=NCYCLE-1
      Z=Z1
      R=R1
      PP(1,1)=PP(1,1)-F*XX(1)
      PP(2,2)=PP(2,2)-F*XX(2)
      PP(3,3)=PP(3,3)-F*XX(3)
      PP(1,2)=PP(1,2)-F*XX(4)
      PP(1,3)=PP(1,3)-F*XX(5)
      PP(2,3)=PP(2,3)-F*XX(6)
      PP(2,1)=PP(1,2)
      PP(3,1)=PP(1,3)
      PP(3,2)=PP(2,3)
      P0=P0-F*XX(7)
      PMIN=PMIN-F*XX(8)
      RETURN
    END SUBROUTINE PFIT
!-----------------------------------------------------------------------
    SUBROUTINE PLOT1 (P,c_string_output_filename) bind(c, name='PLOT1')

        use,intrinsic :: iso_c_binding
        implicit none
        
        real  (c_double), intent(in),  dimension(3*n_ax) :: P
        !integer  (c_int), intent(in)                  :: ILP
        type     (c_ptr), intent(in), target          :: c_string_output_filename  ! the filename for the output. Note that we are defining a pointer to a C Char type variable here

        integer  (c_int) :: ILP = 15

        integer  (c_int), parameter :: NX=100
        integer  (c_int), parameter :: NY=50
      
        character*1, dimension(0:NX,0:NY) :: AA
        integer  (c_int) :: IX
        integer  (c_int) :: IY
        real  (c_double) :: XMIN
        real  (c_double) :: XMAX
        real  (c_double) :: YMIN
        real  (c_double) :: YMAX
        real  (c_double), dimension(3) :: Y
        real  (c_double) :: XI
        real  (c_double) :: YI
        real  (c_double) :: DX
        real  (c_double) :: DY
        integer  (c_int) :: I
        integer  (c_int) :: J
        
        integer  (c_int) :: N

        character(len=128), pointer :: f_string_output_filename
        integer :: string_length_output_filename

        ! Passing a character string from C to Fortran is only straighforward when the the string is of length 1 !!
        ! There are some fixes to this problem in the Fortran 2018 standard, but they have not yet been widely implemented
        ! Here we adopt a solution based on pointers which I dislike, but at least it works. 
      
        ! Do some work to get a filename we can use  ... 
        call c_f_pointer(c_loc( c_string_output_filename ), f_string_output_filename) ! convert c pointer to fortran pointer ... c_loc obtains the c address of the argument.
        string_length_output_filename = index(f_string_output_filename, c_null_char) - 1 ! The text that we need can be found in f_string_output_filename(1:string_length_ouput_filename). Hurrah.      
 
        open(ILP, file=f_string_output_filename(1:string_length_output_filename), &
             status="old", position="append", action="write") ! Open the relevant file, such that output will be appended 

        write(ILP,*)        
        write(ILP,*) "A plot of the Patterson origin peak along the axial directions."
        write(ILP,*) "---------------------------------------------------------------"
        write(ILP,*) 
        
          N=n_ax
                    
          !
          ! BLANK OUT THE PLOT ARRAY.
          !
          DO IX=0,NX
          DO IY=0,NY
            AA(IX,IY)=' '
          END DO
          END DO
          !
          ! FILL IN THE GRID MARKS.
          !
          DO IX=0,NX
            AA(IX, 0)='.'
            AA(IX,NY)='.'
          END DO
          DO IX=0,NX,10
            AA(IX, 0)='+'
            AA(IX,NY)='+'
          END DO
          DO IY=0,NY
            AA( 0,IY)='.'
            AA(NX,IY)='.'
          END DO
          DO IY=0,NY,10
            AA( 0,IY)='+'
            AA(NX,IY)='+'
          END DO
          !
          ! SCALE PLOT AXES.
          !
          XMIN=0
          XMAX=0.5
          YMIN=+1E9
          YMAX=-1E9
          DO I=1,3*N
            YMIN=MIN(YMIN,P(I))
            YMAX=MAX(YMAX,P(I))
          END DO
          !
          ! FILL IN DATA CURVE.
          !
          DX=(XMAX-XMIN)/N
          DO I=1,N
            XI=XMIN+(I-1)*DX   
            IX=NINT(((XI-XMIN)/(XMAX-XMIN))*NX)
            Y(3)=P(I)
            Y(2)=P(I+N)
            Y(1)=P(I+2*N)
            DO J=1,3
              YI=Y(J)
              IY=NINT(((YMAX-YI)/(YMAX-YMIN))*NY)
              IF (IX.GE.0.AND.IX.LE.NX.AND.IY.GE.0.AND.IY.LE.NY) THEN
                IF (J.EQ.1) AA(IX,IY)='W'
                IF (J.EQ.2) AA(IX,IY)='V'
                IF (J.EQ.3) AA(IX,IY)='U'
              END IF
            END DO
          END DO
          !
          ! PRINT THE PLOT ARRAY.
          !
          DY=(YMAX-YMIN)/NY
          DO IY=0,NY
            IF (MOD(IY,10).EQ.0) THEN
              WRITE (ILP,100) YMAX-IY*DY,(AA(IX,IY),IX=0,NX)
            ELSE                        
              WRITE (ILP,101)            (AA(IX,IY),IX=0,NX)
            END IF
          END DO
          DX=(XMAX-XMIN)/10
          WRITE (ILP,102) (XMIN+IX*DX,IX=0,10)
     100  FORMAT (1X,E10.3,101A1)
     101  FORMAT (1X,10X,  101A1)
     102  FORMAT (1X,1X,11F10.2)
     
          WRITE (ILP,'(/1X,''P(U,V,W) ALONG [100], [010], AND [001] FOR 0 .LE. U, V, W .LE. 0.5.'&
                     '//1X,''P(U,V,W) = (2/VCELL)*SUM(H,K,L) [FSQ(MEAS)/SUM(J) F(J)**2] COS [2*PI*(H*U + K*V + L*W)]'')')
          WRITE (ILP,'(/)')
     
          close(ILP)
          
          RETURN
    END SUBROUTINE PLOT1
!-----------------------------------------------------------------------
    SUBROUTINE PLOT2 (U,V,W,P,PMIN,P0,PP,RADIUS,G,c_string_output_filename) bind (c, name='PLOT2')
      
      use,intrinsic :: iso_c_binding
      implicit none

      real  (c_double), intent(in),  dimension(n_pc)   :: U
      real  (c_double), intent(in),  dimension(n_pc)   :: V
      real  (c_double), intent(in),  dimension(n_pc)   :: W
      real  (c_double), intent(in),  dimension(n_pc)   :: P
      real  (c_double), intent(in)                  :: PMIN
      real  (c_double), intent(in)                  :: P0
      real  (c_double), intent(in),  dimension(3,3) :: PP
      real  (c_double), intent(in)                  :: RADIUS
      real  (c_double), intent(in),  dimension(3,3) :: G
      !integer  (c_int), intent(in)                  :: ILP
      type     (c_ptr), intent(in), target          :: c_string_output_filename  ! the filename for the output. Note that we are defining a pointer to a C Char type variable here
      
      integer  (c_int) :: ILP = 15
      
      integer  (c_int), parameter :: NX=100
      integer  (c_int), parameter :: NY=50
      
      character*1, dimension(0:NX,0:NY) :: AA
      
      integer  (c_int) :: IX
      integer  (c_int) :: IY
      real  (c_double) :: XMIN
      real  (c_double) :: XMAX
      real  (c_double) :: YMIN
      real  (c_double) :: YMAX
      real  (c_double), dimension(3) :: R
      real  (c_double) :: XI
      real  (c_double) :: YI
      real  (c_double) :: DX
      real  (c_double) :: DY
      integer  (c_int) :: I
      integer  (c_int) :: J
      integer  (c_int) :: K
      
      integer  (c_int) :: N

      character(len=128), pointer :: f_string_output_filename
      integer :: string_length_output_filename
      
      ! Passing a character string from C to Fortran is only straighforward when the the string is of length 1 !!
      ! There are some fixes to this problem in the Fortran 2018 standard, but they have not yet been widely implemented
      ! Here we adopt a solution based on pointers which I dislike, but at least it works. 
    
      ! Do some work to get a filename we can use  ... 
      call c_f_pointer(c_loc( c_string_output_filename ), f_string_output_filename) ! convert c pointer to fortran pointer ... c_loc obtains the c address of the argument.
      string_length_output_filename = index(f_string_output_filename, c_null_char) - 1 ! The text that we need can be found in f_string_output_filename(1:string_length_ouput_filename). Hurrah.      

      open(ILP, file=f_string_output_filename(1:string_length_output_filename), &
           status="old", position="append", action="write") ! Open the relevant file, such that output will be appended 


       write(ILP,*)        
       write(ILP,*) "Fitting of a Gaussian function to the Patterson Origin Peak"
       write(ILP,*) "-----------------------------------------------------------"
       write(ILP,*) 
      
      N = n_pc
    
      !
      ! BLANK OUT THE PLOT ARRAY.
      !
      DO IX=0,NX
      DO IY=0,NY
        AA(IX,IY)=' '
      END DO
      END DO
      !
      ! FILL IN THE GRID MARKS.
      !
      DO IX=0,NX
        AA(IX, 0)='.'
        AA(IX,NY)='.'
      END DO
      DO IX=0,NX,10
        AA(IX, 0)='+'
        AA(IX,NY)='+'
      END DO
      DO IY=0,NY
        AA( 0,IY)='.'
        AA(NX,IY)='.'
      END DO
      DO IY=0,NY,10
        AA( 0,IY)='+'
        AA(NX,IY)='+'
      END DO
      !
      ! SCALE PLOT AXES.
      !
      XMIN=0
      XMAX=RADIUS
      YMIN=0
      YMAX=P(1)-PMIN
      !
      ! FILL IN THE ORIGIN PEAK PROFILE AND FITTED GAUSSIAN.
      !
      DO I=1,N
        R(1)=U(I)
        R(2)=V(I)
        R(3)=W(I)
        XI=0
        DO J=1,3
        DO K=1,3
          XI=XI+R(J)*R(K)*G(J,K)
        END DO
        END DO
        XI=SQRT(XI)
        IX=NINT(((XI-XMIN)/(XMAX-XMIN))*NX)
        !
        ! FITTED GAUSSIAN
        !
        YI=PUVW(U(I),V(I),W(I),P0,PP)
        IY=NINT(((YMAX-YI)/(YMAX-YMIN))*NY)
        IF (IX.GE.0.AND.IX.LE.NX.AND.IY.GE.0.AND.IY.LE.NY) AA(IX,IY)='*'
        !
        ! PEAK PROFILE
        !
        YI=P(I)-PMIN
        IY=NINT(((YMAX-YI)/(YMAX-YMIN))*NY)
        IF (IX.GE.0.AND.IX.LE.NX.AND.IY.GE.0.AND.IY.LE.NY) AA(IX,IY)='O'
      END DO
      !
      ! PRINT THE PLOT ARRAY.
      !
      DY=(YMAX-YMIN)/NY
      DO IY=0,NY
        IF (MOD(IY,10).EQ.0) THEN
          WRITE (ILP,100) YMAX-IY*DY,(AA(IX,IY),IX=0,NX)
        ELSE                        
          WRITE (ILP,101)            (AA(IX,IY),IX=0,NX)
        END IF
      END DO
      DX=(XMAX-XMIN)/10
      WRITE (ILP,102) (XMIN+IX*DX,IX=0,10)
 100  FORMAT (1X,E10.3,101A1)
 101  FORMAT (1X,10X,  101A1)
 102  FORMAT (1X,1X,11F10.3)
 
      WRITE (ILP,'(/1X,''RADIAL PROJECTION OF THE PATTERSON ORIGIN PEAK "O" SUPERIMPOSED ON THE FITTED GAUSSIAN "*"'')')
      WRITE (ILP,'(/)') 
      
      close(ILP)
      
      RETURN
    END SUBROUTINE PLOT2    
!-----------------------------------------------------------------------
    SUBROUTINE PRELIM (n_data, index_h, index_k, index_l, normalized_intensity, &
                       CELL,VCELL,GINV,P,PMIN,P0,PP,R) bind(c, name='PRELIM')


          !
          ! CALCULATE THE PATTERSON DENSITY ALONG THE A, B, AND C AXES OUT TO
          ! U, V, AND W = 0.5, AND FORM A FIRST APPROXIMATION TO A TRIVARIATE
          ! GAUSSIAN ORIGIN PEAK FUNCTION.
          !

          use,intrinsic :: iso_c_binding
          implicit none
  
          integer  (c_int),  intent(in)                       :: n_data
          integer  (c_int),  intent(in),    dimension(n_data) :: index_h
          integer  (c_int),  intent(in),    dimension(n_data) :: index_k
          integer  (c_int),  intent(in),    dimension(n_data) :: index_l
          real  (c_double),  intent(in),    dimension(n_data) :: normalized_intensity
          real  (c_double),  intent(in),  dimension(6)      :: CELL
          real  (c_double),  intent(in)                     :: VCELL
          real  (c_double),  intent(in),  dimension(3,3)    :: GINV
          real  (c_double), intent(out),  dimension(3*n_ax) :: P  ! Stores the axial Patterson density - along the A,B  and C axes - from 1 through n, n+1 through 2n, 2n+1 through 3n
          real  (c_double), intent(out)                     :: PMIN  ! The minimum Patterson density along the axial directions      
          real  (c_double), intent(out)                     :: P0 ! An estimate for the Amplitude of the trivariate Gaussian function - see eqn 6 of Blessing and Lang
          real  (c_double), intent(out),  dimension(3,3)    :: PP ! An estimate for the Matrix specifying the trivariate Gaussian function - see eqn 6 of Blessing and Lang
          real  (c_double), intent(out)                     :: R  ! The half-height radius (Angstroms) of the Patterson origin peak.
          
          real  (c_double)                 :: DELTA
          real  (c_double)                 :: CONST         
          integer  (c_int)                 :: I
          integer  (c_int)                 :: J
          integer  (c_int)                 :: IH
          integer  (c_int)                 :: IK
          integer  (c_int)                 :: IL
          real  (c_double)                 :: FSQ
          real  (c_double)                 :: FH
          real  (c_double)                 :: FK
          real  (c_double)                 :: FL        
          real  (c_double)                 :: U
          real  (c_double)                 :: V
          real  (c_double)                 :: W         
          integer  (c_int)                 :: IU
          integer  (c_int)                 :: IV
          integer  (c_int)                 :: IW
          real  (c_double)                 :: Q         
          real  (c_double), dimension(3)   :: SIGMA        
          real  (c_double), dimension(3,3) :: QQ
          real  (c_double)                 :: DETQQ

          integer  (c_int)                 :: N
          integer  (c_int)                 :: idata
                                    
          !
          ! SUM THE AXIAL PATTERSON DENSITIES FROM 0 TO 0.5.
          !
          
          N=n_ax
          
          DO I=1,3*N
            P(I)=0
          END DO
          DELTA=0.5/N
          ! rlkPI=3.141593
          CONST=2*PI*DELTA
!rlk          REWIND IFILE
!rlk     1    READ (IFILE,END=9) IH,IK,IL,FSQ
          ! loop over the data supplied by argument - this replaces the read from file
          do idata = 1,n_data
         
            IH = index_h(idata)
            IK = index_k(idata)
            IL = index_l(idata)
            FSQ = normalized_intensity(idata)
            
            FH=CONST*IH
            FK=CONST*IK
            FL=CONST*IL
            DO I=0,N-1
              U=FH*I
              V=FK*I
              W=FL*I
              IU=I+1
              IV=I+N+1
              IW=I+2*N+1
              P(IU)=P(IU)+FSQ*COS(U)
              P(IV)=P(IV)+FSQ*COS(V)
              P(IW)=P(IW)+FSQ*COS(W)
            END DO
            
          end do
!rlk          GO TO 1
     9    DO I=1,3*N
            P(I)=(2/VCELL)*P(I)
          END DO
          !
          ! FIND THE MINIMUM DENSITY.
          !
          PMIN=0
          DO I=1,3*N
            PMIN=MIN(PMIN,P(I))
          END DO
          !
          ! FIND THE HALF-HEIGHT RADII (ANGSTROMS) OF THE ORIGIN PEAK.
          !
          Q=0.5*(P(1)+PMIN)
          U=0
          V=0
          W=0
          DO I=2,N
            IU=I
            IV=I+N
            IW=I+2*N
            IF (U.EQ.0.AND.P(IU).LE.Q) U=IU
            IF (V.EQ.0.AND.P(IV).LE.Q) V=IV
            IF (W.EQ.0.AND.P(IW).LE.Q) W=IW
            IF (U*V*W.NE.0) GO TO 5
          END DO
     5    IU=U
          IV=V
          IW=W
          U=IU-(P(IU)-Q)/(P(IU)-P(IU-1))
          V=IV-(P(IV)-Q)/(P(IV)-P(IV-1))-N
          W=IW-(P(IW)-Q)/(P(IW)-P(IW-1))-2*N
          U=U*DELTA*CELL(1)
          V=V*DELTA*CELL(2)
          W=W*DELTA*CELL(3)
          !
          ! APPROXIMATE THE ORIGIN PEAK FUNCTION, P(U) - PMIN = P0*EXP(-UT*PP*U).
          !
          P0=P(1)-PMIN
          !
          ! THE TRIVARIATE GAUSSIAN EXPONENTIAL ARGUMENT -U(J)*PP(I,J)*U(I)
          ! CORRESPONDS TO THE EXPONENTIAL AGRUMENT OF A UNIVARIATE NORMAL
          ! DISTRIBUTION FUNCTION,
          !
          ! P(X) = 1/(SIGMA*SQRT(2*PI))*EXP(-(X - MU)**2/(2*SIGMA**2)),
          !
          ! AND FROM THE HALF-HEIGHT RADIUS R,
          !
          ! SIGMA = R/SQRT(-2*LN(0.5)).
          !
          SIGMA(1)=U/1.177410
          SIGMA(2)=V/1.177410
          SIGMA(3)=W/1.177410
          DO I=1,3
          DO J=I,3
            QQ(I,J)=2*SIGMA(I)*SIGMA(J)*GINV(I,J)
            QQ(J,I)=QQ(I,J)
          END DO
          END DO
          CALL MATINV3 (QQ,PP,DETQQ)
          !
          ! FOR A NORMAL GAUSSIAN, P(X) < 0.003 IF ABS(X - MU) > 3*SIGMA.
          !
          R=3*MAX(SIGMA(1),SIGMA(2),SIGMA(3))
          RETURN
    END SUBROUTINE PRELIM
!-----------------------------------------------------------------------
    SUBROUTINE PSUM (n_data, index_h, index_k, index_l, normalized_intensity, &
                       VCELL,U,V,W,P) bind(c, name='PSUM')

        !
        ! CALCULATE THE PATTERSON DENSITY ON A PRE-CALCULATED HEMISPHERICAL
        ! POLAR COORDINATE GRID AROUND THE ORIGIN.
        !

        use,intrinsic :: iso_c_binding
        implicit none
  
        integer  (c_int),  intent(in)                       :: n_data
        integer  (c_int),  intent(in),     dimension(n_data) :: index_h
        integer  (c_int),  intent(in),     dimension(n_data) :: index_k
        integer  (c_int),  intent(in),     dimension(n_data) :: index_l
        real  (c_double),  intent(in),     dimension(n_data) :: normalized_intensity
        real  (c_double),  intent(in)                        :: VCELL
        real  (c_double),  intent(inout),  dimension(n_pc)   :: U
        real  (c_double),  intent(inout),  dimension(n_pc)   :: V
        real  (c_double),  intent(inout),  dimension(n_pc)   :: W
        real  (c_double),  intent(out),    dimension(n_pc)   :: P

        integer  (c_int)                     :: I
        integer  (c_int)                     :: IH
        integer  (c_int)                     :: IK
        integer  (c_int)                     :: IL
        real  (c_double)                     :: FSQ

        integer  (c_int)                     :: idata

      integer  (c_int)                     :: N  
      N= n_pc
      
      DO I=1,N
        U(I)=TWOPI*U(I)
        V(I)=TWOPI*V(I)
        W(I)=TWOPI*W(I)
        P(I)=0
      END DO
      
!      REWIND IFILE
! 1    READ (IFILE,END=9) IH,IK,IL,FSQ
      ! loop over the data supplied by argument - this replaces the read from file
      do idata = 1,n_data

        IH = index_h(idata)
        IK = index_k(idata)
        IL = index_l(idata)
        FSQ = normalized_intensity(idata)
      
        DO I=1,N
          P(I)=P(I)+FSQ*COS(IH*U(I)+IK*V(I)+IL*W(I))
        END DO

      end do
!      GO TO 1

 9    DO I=1,N
        U(I)=U(I)/TWOPI
        V(I)=V(I)/TWOPI
        W(I)=W(I)/TWOPI
        P(I)=(2/VCELL)*P(I)
      END DO
      RETURN
    END SUBROUTINE PSUM

!-----------------------------------------------------------------------
    real (c_double) FUNCTION PUVW(U,V,W,P0,PP) bind(c, name='PUVW')
        
      use,intrinsic :: iso_c_binding
      implicit none
      
      real  (c_double), intent(in)                 :: U
      real  (c_double), intent(in)                 :: V
      real  (c_double), intent(in)                 :: W
      real  (c_double), intent(in)                 :: P0
      real  (c_double), intent(in), dimension(3,3) :: PP

      real  (c_double), dimension(3) :: UU
      real  (c_double)               :: Q
      integer  (c_int)                        :: I
      integer  (c_int)                        :: J
                      
      UU(1)=U
      UU(2)=V
      UU(3)=W
      Q=0
      DO I=1,3
      DO J=1,3
        Q=Q+UU(J)*PP(I,J)*UU(I)
      END DO
      END DO
      PUVW=P0*EXP(-Q)
      RETURN
    END FUNCTION PUVW
!-----------------------------------------------------------------------
    real (c_double) FUNCTION SINTHL (IH,IK,IL,GINV) bind(c, name='SINTHL')
    
      use,intrinsic :: iso_c_binding
      implicit none
 
      integer  (c_int)         , intent(in)                 :: IH
      integer  (c_int)         , intent(in)                 :: IK
      integer  (c_int)         , intent(in)                 :: IL
      real  (c_double)         , intent(in), dimension(3,3) :: GINV
      
      real  (c_double), dimension(3)          :: H
      real  (c_double)                        :: Q      
      integer  (c_int)                        :: I
      integer  (c_int)                        :: J
      
      H(1)=IH
      H(2)=IK
      H(3)=IL
      Q=0
      DO I=1,3
      DO J=1,3
        Q=Q+H(J)*GINV(I,J)*H(I)
      END DO
      END DO
      SINTHL=0.5*SQRT(Q)
      RETURN
    END FUNCTION SINTHL
    
!-----------------------------------------------------------------------
    SUBROUTINE SOLVE (P0,P,VCELL,G,GINV,crystal_family,SCALEK,B,U_CIF,U_STAR,BISO,UISO,&
                      c_string_output_filename) bind(c, name='SOLVE')
      
      !
      ! OBTAIN THE SCALE FACTOR, SCALEK, AND OVERALL ANISOTROPIC THERMAL
      ! VIBRATION TENSORS, B(I,J), U_CIF(I,J) and U_STAR(I,J) DEFINED BY
      !
      !   F(ABSOLUTE) = SCALEK*F(RELATIVE)
      !
      ! AND
      !
      !   F(T) = F(0)*EXP(-SUM(I) SUM(J) H(I)*H(J)*B(I,J))
      !
      ! WHERE
      !
      !   B(I,J) = 2*PI**2*ASTAR(I)*ASTAR(J)*U_CIF(I,J))
      !
      !   AND
      !
      !   B(I,J) = 2*PI**2*U_STAR(I,J))
      
      use,intrinsic :: iso_c_binding
      implicit none
      
      real    (c_double), intent(in)                  :: P0
      real    (c_double), intent(in),  dimension(3,3) :: P
      real    (c_double), intent(in)                  :: VCELL
      real    (c_double), intent(in),  dimension(3,3) :: G
      real    (c_double), intent(in),  dimension(3,3) :: GINV
      character (c_char), intent(in),  dimension(1)   :: crystal_family
      real    (c_double), intent(out)                 :: SCALEK
      real    (c_double), intent(out), dimension(3,3) :: B
      real    (c_double), intent(out), dimension(3,3) :: U_CIF
      real    (c_double), intent(out), dimension(3,3) :: U_STAR
      real    (c_double), intent(out)                 :: BISO
      real    (c_double), intent(out)                 :: UISO
      !integer    (c_int), intent(in)                  :: ILP
      type       (c_ptr), intent(in), target          :: c_string_output_filename  ! the filename for the output. Note that we are defining a pointer to a C Char type variable here

      integer  (c_int) :: ILP = 15
      
      real  (c_double), dimension(3,3) :: PINV
      real  (c_double)                 :: DETP
      integer  (c_int)                          :: I
      integer  (c_int)                          :: J
      
      
      character(len=128), pointer :: f_string_output_filename
      integer :: string_length_output_filename
      
      ! Passing a character string from C to Fortran is only straighforward when the the string is of length 1 !!
      ! There are some fixes to this problem in the Fortran 2018 standard, but they have not yet been widely implemented
      ! Here we adopt a solution based on pointers which I dislike, but at least it works. 
    
      ! Do some work to get a filename we can use  ... 
      call c_f_pointer(c_loc( c_string_output_filename ), f_string_output_filename) ! convert c pointer to fortran pointer ... c_loc obtains the c address of the argument.
      string_length_output_filename = index(f_string_output_filename, c_null_char) - 1 ! The text that we need can be found in f_string_output_filename(1:string_length_ouput_filename). Hurrah.      

      open(ILP, file=f_string_output_filename(1:string_length_output_filename), &
           status="old", position="append", action="write") ! Open the relevant file, such that output will be appended 
      
      write(ILP,*)
      write(ILP,*) "Determining scale factor and overall atomic displacement tensor from the fitted Gaussian"
      write(ILP,*)
            
      CALL MATINV3 (P,PINV,DETP)
      SCALEK=SQRT(SQRT(DETP)/(SQRT(PI**3)*VCELL*P0))

!  RLK ... Split the following into two loops, and perform symmetrization in the middle
!      DO I=1,3
!      DO J=I,3
!        B(I,J)=PI**2*PINV(I,J)/2
!        B(J,I)=B(I,J)
!        U(I,J)=B(I,J)/(2*PI**2*SQRT(GINV(I,I)*GINV(J,J)))
!        U(J,I)=U(I,J)
!      END DO
!      END DO

      ! First evaluate the unique elements of B
      DO I=1,3
        DO J=I,3
          B(I,J)=PI**2*PINV(I,J)/2
        END DO
      END DO
 
      write(ILP,fmt="(a,/,6(E16.9,1x),/)") "Unique elements of B before applying symmmetry constraints",&
      B(1,1), B(2,2), B(3,3), B(1,2), B(1,3), B(2,3)
      
      ! Then impose exact symmetry constraints on the elements of B
      ! These symmetry constraints are discussed in many places e.g.
      ! Methods in X-ray crystallography by J.W. Jefferey, Appendix IV
      ! Fundamentals of X-ray crystallography 2nd ed. Giacovazzo et al. Tables 3.8.1 and 3.8.2
      
      if      (crystal_family(1) .eq. "a") then
        write(ILP,fmt =10) "triclinic (a)   "
        
      else if (crystal_family(1) .eq. "m") then
        ! This assumes that the b axis is unique. Different settings will have different constraints
        write(ILP,fmt =20) "monoclinic (m)  ","A     B     C     0     E     0     "
        B(1,2) = 0.0
        B(2,3) = 0.0
        
      else if (crystal_family(1) .eq. "o") then
        write(ILP,fmt =20) "orthorhombic (o)","A     B     C     0     0     0     " 
        B(1,2) = 0.0
        B(1,3) = 0.0
        B(2,3) = 0.0
                      
      else if (crystal_family(1) .eq. "t") then
        write(ILP,fmt =20)  "tetragonal (t) ","A     A     C     0     0     0     "
        B(1,1) = ( B(1,1) + B(2,2) )/2.0
        B(2,2) = B(1,1)
        B(1,2) = 0.0
        B(1,3) = 0.0
        B(2,3) = 0.0
        
      else if (crystal_family(1) .eq. "h") then
        ! This assumes use of the hexagonal setting for the "rhombohedral" space groups. Constraints are different for the rhombohedral setting.
        write(ILP,fmt =20)  "hexagonal (h)  ","A     A     C     A/2   0     0     "
        B(1,1) = ( B(1,1) + B(2,2) + 2*B(1,2))/3.0
        B(2,2) = B(1,1)
        B(1,2) = B(1,1)/2.0
        B(1,3) = 0.0
        B(2,3) = 0.0
        
      else if (crystal_family(1) .eq. "c") then
        write(ILP,fmt =20)  "cubic (c)      ","A     A     A     0     0     0     "
        B(1,1) = ( B(1,1) + B(2,2) + B(3,3))/3.0
        B(2,2) = B(1,1)
        B(3,3) = B(1,1)
        B(1,2) = 0.0
        B(1,3) = 0.0
        B(2,3) = 0.0
        
      else
        write(ILP,*) "Illegal crystal family in function SOLVE. This should never happen"
        
      end if
      
 10   format(/,"The crystal family is ",a16,/," There are no symmetry constaints of the elements of U and B")
 20   format("The crystal family is ",a16,/,&
             "The following symmetry constraints exist on the elements of B (with equivalent constraints on the elements of U)",/,&
             "B11   B22   B33   B12   B13   B23   ",/,a36,/)

      write(ILP,fmt="(a,/,6(E16.9,1x),/)") "Unique elements of B after applying symmetry constraints",&
      B(1,1), B(2,2), B(3,3), B(1,2), B(1,3), B(2,3)

      ! Now complete the evaluation of B and U using the symmetrized values
      ! The representation of U originally returned here was been termed U-CIF by Grosse-Kunstleve RW, Adams PD.  J Appl Cryst. 2002;35(4):477–80. UCIF is defined in eqns 4 and 5
      ! See also Trueblood  et al Acta Cryst (1996) A52, 770-781 eqn 25
      ! The more conventional representation U_STAR is also now returned
       
      DO I=1,3
        DO J=I,3
          B(J,I)=B(I,J)
          U_STAR(I,J)=B(I,J)/(2*PI**2) 
          U_CIF(I,J)=U_STAR(I,J)/SQRT(GINV(I,I)*GINV(J,J)) 
          U_STAR(J,I)=U_STAR(I,J)
          U_CIF(J,I)=U_CIF(I,J)
        END DO
      END DO
      
      ! Compute the isotropic temperature factor equivalents

      BISO=0
      DO I=1,3
      DO J=1,3
        BISO=BISO+G(I,J)*B(I,J) ! This is the trace of B*G = Tr(B*G)
      END DO
      END DO
      BISO=BISO*4/3 ! 4/3 Tr(B*G) see Giacovazzo, Fundamentals of Crystallography eqn 3.B.10c
      UISO=BISO/(8*PI**2) ! Tr(B*G) / 6π^2 see  Giacovazzo, Fundamentals of Crystallography eqn 3.B.10d   

      close(ILP)
      
      RETURN
    END SUBROUTINE SOLVE
!-----------------------------------------------------------------------
    SUBROUTINE WEIGHT (P,P0,PP,U,V,W,Z,WT) bind(c, name='WEIGHT')
        
      !
      ! DOWN-WEIGHT OFF-ORIGIN PEAKS.
      !

      use,intrinsic :: iso_c_binding
      implicit none

      real  (c_double),  intent(in)                  :: P
      real  (c_double),  intent(in)                  :: P0
      real  (c_double),  intent(in),  dimension(3,3) :: PP
      real  (c_double),  intent(in)                  :: U
      real  (c_double),  intent(in)                  :: V
      real  (c_double),  intent(in)                  :: W
      real  (c_double),  intent(in)                  :: Z
      real  (c_double),  intent(out)                 :: WT
       
      IF (P.LT.0.5*P0.AND.P-PUVW(U,V,W,P0,PP).GT.Z) THEN
        WT=0.1
      ELSE
        WT=1
      END IF
      RETURN
    END SUBROUTINE WEIGHT
!-----------------------------------------------------------------------
    SUBROUTINE PRINT_SOLUTION (RADIUS, NCYCLE, R, SCALEK, B, U, BISO, UISO, &
                               c_string_output_filename) bind(c, name='PRINT_SOLUTION')
      ! A subroutine abstracted from the main program of the original Fortran 77 code.

      use,intrinsic :: iso_c_binding
      implicit none
      
      real  (c_double), intent(in)                 :: RADIUS
      integer  (c_int), intent(in)                 :: NCYCLE
      real  (c_double), intent(in)                 :: R
      real  (c_double), intent(in)                 :: SCALEK
      real  (c_double), intent(in), dimension(3,3) :: B
      real  (c_double), intent(in), dimension(3,3) :: U
      real  (c_double), intent(in)                 :: BISO
      real  (c_double), intent(in)                 :: UISO
!      integer  (c_int), intent(in)                 :: ILP
      type     (c_ptr), intent(in), target         :: c_string_output_filename  ! the filename for the output. Note that we are defining a pointer to a C Char type variable here

      integer  (c_int) :: ILP = 15

      integer  (c_int)                 :: N

      character(len=128), pointer :: f_string_output_filename
      integer :: string_length_output_filename

      ! Passing a character string from C to Fortran is only straighforward when the the string is of length 1 !!
      ! There are some fixes to this problem in the Fortran 2018 standard, but they have not yet been widely implemented
      ! Here we adopt a solution based on pointers which I dislike, but at least it works. 
    
      ! Do some work to get a filename we can use  ... 
      call c_f_pointer(c_loc( c_string_output_filename ), f_string_output_filename) ! convert c pointer to fortran pointer ... c_loc obtains the c address of the argument.
      string_length_output_filename = index(f_string_output_filename, c_null_char) - 1 ! The text that we need can be found in f_string_output_filename(1:string_length_ouput_filename). Hurrah.      

      open(ILP, file=f_string_output_filename(1:string_length_output_filename), &
           status="old", position="append", action="write") ! Open the relevant file, such that output will be appended 

      N = n_pc
      
      WRITE (ILP,692) RADIUS,N,NCYCLE,R
      WRITE (ILP,693) SCALEK,B(1,1),B(2,2),B(3,3),B(1,2),B(1,3),B(2,3), &
      U(1,1),U(2,2),U(3,3),U(1,2),U(1,3),U(2,3),UISO,BISO
 692  FORMAT ('OVERALL RESULTS:'/' --------------'/&
     'RADIUS OF FIT OF THE TRIVARIATE GAUSSIAN TO THE PATTERSON ORIGIN PEAK'/' RADIUS = ',F5.2,' ANGSTROM'/&
     'NUMBER OF GRID POINTS CALCULATED'/' NPOINTS = ',I6/'STATISTICS-OF-FIT:'/&
     'R = SQRT(SUM(WT*(POBS - PCALC)**2)/SUM(WT*POBS**2))'/&
     'NCYCLE = ',I10/' R      = ',E10.3)
 693  FORMAT ('SCALEK, B(I,J), AND U(I,J) ARE DEFINED ACCORDING TO'/&
     '  F(ABSOLUTE) = SCALEK*F(RELATIVE)'/'AND'/'  f(T) = f(0)*EXP(-SUM(I) SUM(J) H(I)*H(J)*B(I,J))'/'WHERE'/&
     '  B(I,J) = 2*PI**2*ASTAR(I)*ASTAR(J)*U(I,J))'/'      SCALEK'/' ',E12.4/&
     '      B(1,1)      B(2,2)      B(3,3)      B(1,2)      B(1,3)      B(2,3)'/' ',6E12.4/&
     '      U(1,1)      U(2,2)      U(3,3)      U(1,2)      U(1,3)      U(2,3)'/' ',6E12.4/&
     'EQUIVALENT ISOTROPIC TEMPERATURE FACTORS:'/&
     ' UISO = <U**2>         = ',E12.4,' SQUARED ANGSTROMS'/&
     ' BISO = 8*PI**2*<U**2> = ',E12.4,' SQUARED ANGSTROMS')
      
     close(ILP)
     
    END SUBROUTINE PRINT_SOLUTION 
    
end module rogers_analysis
