!-----------------------------------------------------------------
!     Main part. 
!     Reads input / Sets defaults
!     Call NEOART.
!     Output
!----------------------------------------------------------------- 
 
 program main

      use init, only : get_ns_and_ncm, set_defaults, read_input, l_test, l_scan
      use memory, only : alloc_all!, dealloc_all
      use init, only : bn,e,q,t, cff, ns ,nc, den, ic,vnlin, check
      use collision, only :tau
      use friction, only : fric_coeff
      use geometry, only : rbt, dpsidr
      IMPLICIT NONE  
      
      integer :: i,j
      real :: norm
      real :: bootstrap
      call get_ns_and_ncm()  !read number of species and charge states from input
      call alloc_all()       !allocate and initialize all dynamic arrays
      call set_defaults()    !set_defaults
      call read_input()      !read all input parameters
      call vnlin()           ! reads vlin_drive from file
      call check()           ! check species namelist and limited functionality
      call fric_coeff()      !calculate la and lab, also calls geom atm.
      call neoart(cff)       !main calculation
      
      if(l_test) then
      write(*,*) 'TESTCASE, (set l_test.false (==default) if not running test1 or test2)'       
      write(*,*) 'HEAT FLUX'
      do i= 1, ns
        write(*,*) 'by species', i
        do j= 1, nc(i)
          write(*,*) 'charge state', j
          if (ic==1) then
          NORM = 1.6E-22*BN**2*SQRT(E**3)/(2*(Q)**2 &
     &       *T(i)*TAU(i,2))  
          end if
          if (ic==2) then
          NORM = 1.6E-22*BN**2/(4*(Q)**2 &
     &       *T(i)*TAU(i,2))
          end if
          if (i.ne.1) norm=norm*sqrt(2.)
          write(*,*) cff(i,j,2)*NORM

        end do
      end do
      write(*,*) 'THE BOOTSTRAP CURRENT'
      bootstrap=0
      do i= 1, ns
        do j= 1, nc(i)
          bootstrap = bootstrap + cff(i,j,3)
        end do
      end do
      
      NORM = SQRT(E)*BN/(1600.*Q*DEN(1,1)*T(1))
      write(*,*) bootstrap*norm
      write(*,*)
      
      write(*,*) 'THE POLOIDAL FLOW'
      NORM = 1.E-3*BN**2*E/(Q*T(2))
       write(*,*) cff(2,1,4)*norm
       write(*,*)

      else if(l_scan) then   
      bootstrap=0
      do i= 1, ns
        do j= 1, nc(i)
          bootstrap = bootstrap + cff(i,j,3)
        end do
      end do
!         999   FORMAT(F9.3)       
!         write(*,999,advance='no') q
!         write(*,999,advance='no') e
!         write(*,999) bootstrap*1.E-6
      999   FORMAT(F9.4)       
!      write(*,999,advance='no') 2*1.2566E-6*DEN(2,1)*T(2)*1.6E3/BN/BN
      write(*,999) bootstrap*1.E-6
      
      else
      
      write(*,*) rbt, dpsidr
      write(*,*)
      write(*,*) 'THE BOOTSTRAP CURRENT in MA/m^2'
      bootstrap=0
      do i= 1, ns
        do j= 1, nc(i)
          bootstrap = bootstrap + cff(i,j,3)
        end do
      end do 
       write(*,999) bootstrap*1.E-6
       write(*,*)
       
       end if
       
       



end program      
      
      
