!-----------------------------------------------------------------
!     Main part. 
!     Reads input / Sets defaults
!     Call NEOART.
!     Output
!----------------------------------------------------------------- 
 
 program main

 
      use init, only : read_input, neo_abort
      use init, only :  ns,eps, nreg, sigma, nleg, nenergy, &
                       & ncof, neogeo, neofrc, ic, ncm,  &
                       & nc, zsp, m, T, den ,ds, ISEL,      &
                       & ISHOT, rho, e, q, Rn, Bn, eparr
      use collision, only : tau, xi
      use collision, only :colxi
      IMPLICIT NONE
      
      integer :: NMAXGR
      
      integer :: i,j,k, ierr

            

      PARAMETER(NMAXGR = 1000) 
        
      
      real, save, allocatable :: CFF1(:,:,:),CFF2(:,:,:),CFF3(:,:,:),CFF4(:,:,:)
      real, save, allocatable :: COEFF(:,:,:)
      real :: norm
       
      call read_input()
      
      allocate(COEFF(ns,ncm,4),stat=ierr)
      if (ierr /= 0) call neo_abort('could not allocate array')
      coeff=.0
      allocate(CFF1(ns,ncm,4),stat=ierr)
      if (ierr /= 0) call neo_abort('could not allocate array')
      cff1=.0
      allocate(CFF2(ns,ncm,4),stat=ierr)
      if (ierr /= 0) call neo_abort('could not allocate array')
      cff2=.0
      allocate(CFF3(ns,ncm,4),stat=ierr)
      if (ierr /= 0) call neo_abort('could not allocate array')
      cff3=.0
      allocate(CFF4(ns,ncm,4),stat=ierr)
      if (ierr /= 0) call neo_abort('could not allocate array')
      cff4=.0

      
      if(isel==2) then
!     COPY THEM INTO THE VALUES USED BY THE CODE
       CALL CIRCGEOM(1,RHO,RN,E,Q,BN)
!     USE THE CIRCULAR GEOMETRY APPROXIMATION
      else 
       call neo_abort('understand and implement other geometry')
      end if
     

      do k = 1, 13
 
        E = 0.04 + 0.24*(K-1.)/12. 

!       COPY GEOM VALUES INTO THE VALUES USED BY THE CODE
        CALL CIRCGEOM(1,RHO,RN,E,Q,BN)
  
        NEOFRC = .false. 
        NEOGEO = .TRUE.
        
          CALL NEOART(NEOGEO,NEOFRC,IC,CFF4)
  
 
          CALL COLXI(ns,nc,ncm,ZSP,DEN,T,M)
  
         NORM = 1.6E-22*BN**2*SQRT(E**3)/(2*(Q)**2 &
      &         *T(2)*TAU(2,2))
 
         WRITE(*,*) E, SQRT(2.)*CFF4(2,1,2)*NORM
 
       end do
      
!       deallocate(coeff)
!       deallocate(cff1)
!       deallocate(cff2)
!       deallocate(cff3)
!       deallocate(cff4)
!       deallocate(xi)
!       deallocate(tau)
!       deallocate(nc)
!       deallocate(zsp)
!       deallocate(m)
!       deallocate(T)
!       deallocate(den)
!       deallocate(ds)

      STOP 0
      
end program      
      
      
