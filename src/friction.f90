module friction

     implicit none
     real, save, allocatable :: la(:,:,:), lab(:,:,:,:)

contains

subroutine fric_coeff()
  
      use init, only : ns,nc,ncm,t,den,m,zsp
      use init, only : rho, sigma
      use init, only : nenergy,ncof,ic,nleg,nreg
      use geometry, only : geom
      use geometry, only : bav, b2av, bi2a, rbt, bgradp
      use geometry, only : dpsidr, rnq, fc, gclass, fm
      use geometry, only : mmx, r2i   
      use collision, only : tau, xi
      use collision, only : colxi
   implicit none
   
   real :: dum, ta,tb,ma,mb
   integer :: i,j,k,l,ierr
   real, dimension(3,3) :: mm,nn


       if (RHO.EQ.0.) CALL PERR(2)

       if ((IC.EQ.1).OR.(IC.EQ.3)) THEN
        if (NREG.NE.0) CALL PERR(11)
        if (NLEG.NE.3) CALL PERR(12)
        if (NENERGY.NE.1) CALL PERR(13)
       end if
       if ((IC.EQ.2).OR.(IC.EQ.3)) THEN
        DUM = 1.
        do i = 1, 4
          DUM = DUM*SIGMA(I)
        end do
        if (DUM.NE.1.) CALL PERR(15)
       end if
       if (NCOF.NE.1) CALL PERR(14)


      call geom()
! 
      call colxi()
! 
!     SWITCH OF ION - ELECTRON COLLISIONS? 
      if (NCOF.EQ.0) THEN
!       YES SWITCH OF
        j = 0
        do i = 1, ns
          if (ZSP(i,1).LT.0.) J = I
        end do
        if (j.NE.0) then
         do i = 1, ns
          if (J.NE.I) TAU(I,J) = 0.
         end do
        end if
      end if

!    CALCULATED THE FRICTION COEFFICIENTS
      la=0.
      do i = 1, ns
        do j = 1, ns
          TA = T(I)
          TB = T(J)
          MA = M(I)
          MB = M(J)
           CALL MENN(TA,TB,MA,MB,MM,NN)
           do k = 1, 3
             do l = 1, 3
                LAB(K,L,I,J) = TAU(I,J)*NN(K,L)
                LA(K,L,I) = LA(K,L,I) + TAU(I,J)*MM(K,L)
             end do
           end do
        end do
       end do

end subroutine



end module friction
