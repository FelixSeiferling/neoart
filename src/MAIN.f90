!-----------------------------------------------------------------
!     Main part. 
!     Reads input / Sets defaults
!     Call NEOART.
!     Output
!----------------------------------------------------------------- 
 
 program main

      use init, only : get_ns_and_ncm, set_defaults, read_input 
      use memory, only : alloc_all, dealloc_all
      use init, only : bn,e,q,t, cff
      use collision, only :tau
      use friction, only : fric_coeff

      IMPLICIT NONE  
      
      
      real :: norm
      call get_ns_and_ncm()  !read number of species and charge states from input
      call alloc_all()       !allocate and initialize all dynamic arrays
      call set_defaults()    !set_defaults
      call read_input()      !read all input parameters
      call fric_coeff()      !calculate la and lab, also calls geom atm.
      call neoart(cff)       !main calculation
 
 
      NORM = 1.6E-22*BN**2*SQRT(E**3)/(2*(Q)**2 &
     &         *T(2)*TAU(2,2))
 
      WRITE(*,*) CFF
   !   WRITE(*,*) E, SQRT(2.)*CFF(2,1,4)*NORM

  !    call dealloc_all()
      STOP 0
      
end program      
      
      