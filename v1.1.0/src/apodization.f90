      module apodization
      
      use,intrinsic :: iso_c_binding
      implicit none
  
      contains

        subroutine get_gaussian_width(smax, target_integral, sigma) bind(c)

!         This subroutine returns the standard deviation (sigma) of a Gaussian function with mean zero and unit height 
!         whose integral between 0 and smax is given by target_integral

!         In other words we are solving gaussian_integral(sigma,smax) - target_integral = 0, for sigma

!         Obviously 0 < target_integral < smax or we are in deep trouble 
!         As we do this only once and don't care greatly about speed, we'll use a simple bisection (half-interval) method to locate the solution
!         We can safely assume the solution for sigma is bracketed between 10^-4 and 10^2 ...
!   
         
          use,intrinsic :: iso_c_binding  
          implicit none
          
          integer, parameter :: single = selected_real_kind(p=6)
          integer, parameter :: double = selected_real_kind(p=15)

          real  (c_double), intent(in)     :: smax
          real  (c_double), intent(in)     :: target_integral
          real  (c_double), intent(out)    :: sigma          

          integer, parameter :: n_iterations = 75 ! This will take us down to machine precision for sure
          real  (c_double), parameter :: l_bracket = 1.0E-04
          real  (c_double), parameter :: r_bracket = 1.0E+02

          real  (c_double)            :: fl
          real  (c_double)            :: fm
          real  (c_double)            :: fr

          real  (c_double)            :: root_estimate
          real  (c_double)            :: ds
          real  (c_double)            :: midpoint
          integer                     :: i
           
          fl = gaussian_integral(l_bracket,smax) - target_integral ! should be negative
          fr = gaussian_integral(r_bracket,smax) - target_integral ! should be positive
!          write(*,*) fl,fr, fl*fr 
                  
          if ( fl*fr > 0 ) then                       
            write(*,*) "Something has gone horribly wrong in subroutine get_gaussian_width ... stopping"
            stop
          end if
          
          root_estimate = fl ! take the left bracket as our initial estimate of the root, function is negative here
          ds = r_bracket - l_bracket
          
          do i=1,n_iterations
            ds = ds*0.5 
            midpoint = root_estimate + ds
            fm = gaussian_integral(midpoint,smax) - target_integral
            if (fm <= 0) root_estimate = midpoint 
!            write(*,*) i, ds, midpoint, fm, root_estimate, root_estimate + 2*ds
          end do
          
          sigma = root_estimate
                                              
          return
             
        end subroutine get_gaussian_width
        
        real (c_double) function gaussian_integral(sigma, smax) bind(c)

!       computes the integral of a Gaussian function (with mean zero, unit height, and standard deviation sigma) between 0 and smax
!       Uses the Fortran2008 intrinsic function Erf

          use,intrinsic :: iso_c_binding  
          implicit none
          
          real  (c_double), intent(in)     :: sigma 
          real  (c_double), intent(in)     :: smax
         
          real  (c_double) :: pi
          
          pi = acos(-1.0)    
          gaussian_integral =  sigma * Sqrt(Pi/2.0) * Erf( smax/(Sqrt(2.0)*sigma) )

        end function gaussian_integral


      end module apodization
      