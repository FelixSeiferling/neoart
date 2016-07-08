
     
      
! C--------------------------------------------------------------------
! C     THIS SUBROUTINE CALCULATES THE ENERGY DEPENDENT COLLISION
! C     FREQUENCY AND THE ENERGY DEPENDT ENERGY SCATTERING.
! C
! C     INPUT  : NSM    THE MAXIMUM NUMBER OF SPECIES
! C              NCM    THE MAXIMUM NUMBER OF CHARGE STATES
! C              NS     THE ACTUAL NUMBER OF SPECIES
! C              js     THE SPECIES NUMBER FOR WHICH THE COLLISION 
! C                     FREQUENCY IS TO BE CALCULATED
! C              jc     THE CHARGE STATE OF THE SPECIES FOR WHICH 
! C                     THE COLLISION FREQUENCY IS TO BE CALCULATED.
! C              TAU    NORMALIZED COLLISION FREQUENCY (N_I M_I /
! C                     TAU_IJ)  ARRAY(NSM,NSM)
! C              M      ARRAY(NSM) THE MASS OF THE SPECIES IN KG
! C              T      ARRAY(NSM) THE TEMPERATURE OF THE SPECIES 
! C                     IN KEV
! C              XI     ARRAY(NSM,NCM) THE RELATIVE WEIGHT OF EVERY
! C                     CHARGE STATE
! C              NC     ARRAY(NSM) THE NUMBER OF CHARGE STATES FOR 
! C                     EVERY SPECIES
! C              X      NORMALIZED (TO THERMAL) VELOCITY FOR WHICH
! C                     THE ENERGY DEPENDENT COLLISION FREQUENCY IS
! C                     TO BE CALCULATED
! C              DEN    ARRAY(NSM,NCM) DENSITY OF EVERY COMPONENT 
! C                     IN 10^19 M^-3
! C     OUTPUT   NUD    THE PITCH ANGLE SCATTERING FREQUENCY
! C              NUE    THE ENERGY SCATTERING FREQUENCY
! C
! C     THE ROUTINE CALLS THE FOLLOWING FUNCTION
! C     ERF   : CALCULATES THE ERROR FUNCTION 
! C
! C--------------------------------------------------------------------
subroutine viscol(js, jc, X, NUD, NUE)
      use init, only: ns,nc,ncm,m,t,den
      use collision, only : tau, xi
      implicit none

      integer :: i, js,jc
      real :: x,nud, nue, erf
      real :: PH, G, PI, XAB


      PI = 4.*ATAN(1.)

      NUD = 0.
      NUE = 0.

!    THE LOOP OVER SPECIES
      do i = 1, ns

!       CALCULATE XAB = VTHB / VTHA
        XAB = SQRT(M(js)*T(i)/(M(i)*T(js)))

        PH = ERF(X/XAB)
        G  = (PH - 2*X*EXP(-(X/XAB)**2)/(XAB*SQRT(PI))) &
     &       /(2*(X/XAB)**2)

        NUD = NUD + TAU(js,I) * ( PH - G )/ X**3
        NUE = NUE + TAU(js,I) * (-2.*PH/X**3 + 4.*(T(js)/T(I) + &
     &        1./XAB**2) *G/X)

      end do
      
      NUD = NUD * 0.75 * SQRT(PI) * XI(js,jc) / (1.E19*DEN(js,jc)* &
     &      M(js))
      NUE = NUE * 0.75 * SQRT(PI) * XI(js,jc) / (1.E19*DEN(js,jc)* &
     &      M(js))

      return
end
      
