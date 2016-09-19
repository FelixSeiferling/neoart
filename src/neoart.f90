!--------------------------------------------------------------------
!     SUBROUTINE THIS ROUTINE CALCULATES THE TRANSPORT COEFFICIENTS
!     OF NEOCLASSICAL THEORY. 
!
!     INPUT   NS        : THE NUMBER OF SPECIES
!             NC        : ARRAY(NS) THAT GIVES THE NUMBER OF 
!                         CHARGE STATES PER SPECIES
!             NAR       : ACTUAL (FIRST) DIMENSION OF THE ZSP,DEN, 
!                         AND DS ARRAY, USED FOR CHECKING CONSISTENCY
!             NZM       : ACTUAL (SECOND) DIMENS. OF ARRAY (CONSIS.)
!             ZSP       : ARRAY(NAR,NZM) OF REAL THAT CONTAINS
!                         THE CHARGE NUMBER OF EVERY COMPONENT
!             M         : ARRAY(NS) THAT CONTAINS THE MASS OF 
!                         EVERY SPECIES IN KG
!             T         : ARRAY(NS) THAT CONTAINS THE TEMPERATURE
!                         OF THE SPECIES IN KEV
!             DEN       : ARRAY(NAR,NZM) THAT CONTAINS THE DENSITY
!                         OF EVERY COMPONENT IN UNITS OF 1 10^19
!             DS        : ARRAY(NAR,NZM,2) THAT CONTAINS THE FIRST
!                         AND SECOND THERMODYNAMIC FORCE (PRESSURE
!                         GRADIENT AND TEMPERATURE GRADIENT)
!
!                                          D LN P_IJ 
!                              DS(I,J,1) = ---------   
!                                            D rho
!
!                                          D LN T_I
!                              DS(I,J,2) = --------
!                                           D rho
!
!                         WHERE RHO IS THE FLUX LABEL ACTUALLY USED
!                         THE FLUXES ARE ALSO CALCULATED IN THIS 
!                         COORDINATE. FOR THIS TO BE POSSIBLE THE 
!                         GEOMETRY ROUTINE MUST SUPPLY DPSIDR = 
!                         D PSI / D RHO WHERE PSI IS THE POLOIDAL 
!                         FLUX. 
!             RHO       : THE FLUX SURFACE LABEL
!             EPS       : REQUIRED ACCURACY.
!             ISEL      : PARAMETER THAT SELECTS HOW THE GEOMETRY
!                         DEPENDENT PARAMETERS ARE CALCULATED 
!                         1 ASSUME HAMADA COORDINATES AND USE THE
!                           ROUTINE VISGEOM
!                         2 ASSUME CIRCULAR GEOMETRY AND USE CIRC-
!                           GEOM TO DETERMINE THE PARAMETERS
!                         3 READ FROM FILE
!                         4 READ FROM CHEASE OUTPUT, IN THAT CASE
!                           THE FLUX SURFACE LABEL IS 
!                           RHO=(RMAX-RMIN)/2/R0EXP as in GKW
!             ISHOT     : SHOT NUMBER, ONLY USED WHEN ISEL = 3
!             NREG      : PARAMETER THAT FORCES A CERTAIN REGIME
!                         IN THE CALCULATION OF THE VISCOSITY.
!                         0     THE NORMAL VALUE. THE BANANA AND
!                               PFIRSCH SCHLUETER CONTRIBUTIONS
!                               ARE WEIGHTED
!                         1     THE BANANA REGIME IS FORCED.
!                         OTHER THE PFIRSCH SCHLUETER REGIME IS
!                               FORCED.
!                         NOTE : FORCING THE BANANA REGIME IS NOT 
!                                THE SAME THING AS CALCULATING THE
!                                BANANA PLATEAU CONTRIBUTION. FOR
!                                HIGH COLLISION FREQ. THIS CONTRI-
!                                BUTION USUALLY DECREASES. WHEN
!                                NREG = 1 THERE IS NO SUCH DECR.
!             SIGMA(4)    THE SIGMA'S DETERMINE WHICH OF THE COUPLING 
!                         TERMS IN THE EQUATION FOR THE PFIRSCH 
!                         SCHLUETER REGIME ARE TAKEN INTO ACOUNT. 
!                         SIGMA(1)= 0,1  IS THE CROSS COUPLING BETWEEN 
!                                     THE HEAT EQUATION AND THE EQUATION
!                                     FOR THE THIRD LAGUERRE HARM.
!                         SIGMA(2)= 0.,1 IS THE CROSS COUPLING BETWEEN 
!                                     THE EQUATION FOR THE THIRD HARMONIC
!                                     AND THE TEMPERATURE PERTURB.
!                         SIGMA(3)= 0,1 IS THE ENERGY EXCHANGE BETWEEN 
!                                     THE DIFFERENT SPECIES
!                         SIGMA(4)= 0,1 IS THE COUPLING TO THE "DENSITY 
!                                     PERTURB." IN THE EQUATION FOR THE 
!                                     THRIRD LAGUERRE HARMONIC
!                         IN THE APPROXIMATIONS OF HIRSHMANN AND SIGMAR
!                         SIGMA(1)=0 SIGMA(2)=0 SIGMA(3)=1 SIGMA(4)=1
!                         THAT IS THE FIRST TWO TERMS ARE NEGLECTED.
!                         THE USUAL VALUES OF ALL THESE COEFFICIENTS IS 1.
!             NLEG      : NUMBER OF LEGENDRE POLYNOMALS IN THE 
!                         EXPANSION (MAXIMUM AND NORMAL VALUE IS 3)
!             NENERGY   : PARAMETER THAT DETERMINES WHETHER ENERGY
!                         SCATTERING IS TAKEN INTO ACCOUNT IN THE 
!                         CALCULATION OF THE VISCOSITY 
!                         0 NOT ACCOUNTED FOR 
!                         1 ACCOUNTED FOR (NORMAL VALUE)
!             NCOF        PARAMETERS THAT DETERMINES WHETHER ION-
!                         ELECTRON COLLISIONS ARE TAKEN INTO ACCOUNT
!                         0 NOT TAKEN INTO ACCOUNT
!                         1 TAKEN INTO ACCOUNT (NORMAL VALUE, 1) 
!                         0 ONLY TO OBTAIN SOME ANALYTIC RESULTS)
!             NEOGEO      LOGICAL IF TRUE THEN THE GEOMETRY DEP.
!                         PARAMETERS ARE RECALCULATED.
!             NEOFRC      LOGICAL IF TRUE THEN THE FRICTION AND 
!                         VISCOSITY MATRIX ARE NOT NEWLY CALCULATED
!             IC        : THE CONTRIBUTION FOR WHICH THE COEFFICIENTS
!                         ARE CALCULATED. 
!                         0 THE CLASSICAL PARTICLE FLUX
!                         1 BANANA PLATEAU CONTRIBUTION
!                         2 PFIRSCH SCHLUETER CONTRIBUTION
!                         3 BOTH BANANA-PLATEAU, CLASSICAL, AND P.S.
!             EPARR     : THE PARALLEL ELECTRIC FIELD TIMES THE 
!                         MAJOR RADIUS, I.E. THE LOOP VOLTAGE. 
!             COEFF     : ARRAY(NSM,ncm,4) THAT CONTAINS THE TRANSP.
!                         COEFFICIENTS OF SPECIES NSM,ncm THE LAST
!                         INDEX GIVES
!                         1 PARTICLE FLUX   \GAMMA_IJ
!                         2 ENERGY FLUX     Q_IJ / T_I
!                         3 PARALLEL FLOW   IN UNITS OF CURRENT
!                                           <J.B><B>/<B**2>
!                         4 POLOIDAL FLOW.  <U.VT>/<B.VT>
!  
!     THE ROUTINE CALLS THE FOLLOWING ROUTINES
!     PERR    :  ERROR HANDLING
!     GEOM    :  CALCULATE THE GEOMETRY DEPENDENT QUANTITIES
!     MENN    :  CALCULATE THE FRICTION COEFFICIENTS
!     PS      :  CALCULATE THE PFIRSCH SCHLUETER CONTRIBUTION
!     BP      :  CALCULATE THE BANANA PLATEAU CONTRIBUTION
!--------------------------------------------------------------------
subroutine neoart(coeff)
     
     
     use init, only :  ns,nc, ncm, zsp, m, t, den, ds, ic,    &
                       & rho ,eps, ISEL, ISHOT, nreg, & 
                       & sigma, nleg, nenergy, ncof, eparr, &
                       & vnlin_drive
      use collision, only : tau, xi
      use collision, only : colxi
      use friction, only : la, lab
      use friction, only : fric_coeff
      use error, only : neo_abort, neo_warn
      use geometry, only : geom
      use geometry, only : bav, b2av, bi2a, rbt, bgradp
      use geometry, only : dpsidr, rnq, fc, gclass, fm
      use geometry, only : mmx, r2i
      IMPLICIT NONE     
 
      integer :: i,j,k,l,ikeep, ierr
      real :: EPARN, DUM
      real, dimension(ns,ncm,4) :: coeff, coeffc
      real, dimension(ns,ncm,3) :: uai
   
      
      coeff=0.
      coeffc=0.
     !! confusing renormalization:
     !! 1.E-3 : Gradients are in keV
     !! Houlberg gradient terms with 2*pi*R B_t/dpsidr; we multiply whole eq. with dpsidr/rbt
     !! i can't find the 2pi
!     NORMALIZE THE EPARR   
      EPARN = 1.E-3*DPSIDR*R2I*EPARR 
!     renormalize vnlin terms to match other terms
      do i = 1, ns ; do j = 1, nc(i) ; do k = 1, 3
        vnlin_drive(i,j,k) = vnlin_drive(i,j,k)*1.E-3*dpsidr/rbt
      end do ; end do ; end do 
! ONLY BANANA PLATEAU AT THE MOMENT   !!!!


! !     CALCULATE THE CLASSICAL TRANSPORT CONTRIBUTION
!       if ((IC.EQ.0).OR.(IC.EQ.3)) THEN
!       if(isel==1) call neo_warn('classical contribution set to 0 for isel=1. &
!                     &            not yet implemented')
!       call class(coeffc)
! !       ADD THE CORRECT NORMALIZATION FACTOR, PUT RESULT 
! !       IN COEFF AND CLEAR COEFFC
!         do i = 1, ns
!           do j = 1, nc(i)
!             do k = 1, 2
!               COEFF(I,J,K) = COEFF(I,J,K)+ COEFFC(I,J,K)* &
!      &          GCLASS/(1.6022E-22*ZSP(I,J)*DPSIDR**2) 
!               COEFFC(I,J,K) = 0.
!             end do
!           end do
!         end do
!       end if     
!    
! 
! !     CALCULATE THE PFIRSCH SCHLUETER CONTRIBUTION
!       coeffc=0.
!       if ((IC.EQ.2).OR.(IC.EQ.3)) THEN
!         call PS(uai,coeffc)
! !       ADD THE CORRECT NORMALIZATION FACTOR, PUT RESULT IN 
! !       COEFF AND CLEAR COEFFC
!         do i = 1, ns
!           do j = 1, nc(i)
!             do k = 1, 2
!               coeff(i,j,k) = coeff(i,j,k) - coeffc(i,j,k)* &
!      &          (1.- B2AV*BI2A)*(RBT/DPSIDR)**2/ &
!      &          (1.60E-22*ZSP(I,J)*b2av)
!               coeffc(i,j,k) = 0.
!             end do
!           end do
!         end do
!       end if
! 
! 
      coeffc=0.
      if ((IC.EQ.1).OR.(IC.EQ.3)) then
       call bp(eparn,coeffc)
!       ADD THE CORRECT NORMALIZING FACTOR. PUT RESULT IN 
!       COEFF AND CLEAR COEFFC
        do i = 1, ns
          do j = 1, nc(i)
            do k = 1, 2
              coeff(i,j,k) = coeff(i,j,k) - coeffc(i,j,k)* &
     &          (RBT/DPSIDR)**2/(1.60E-22*ZSP(I,J)* &
     &           b2av)
              coeffc(i,j,k) = 0.
            end do
          end do
        end do
      end if 

        do i = 1, ns
          do j = 1, nc(i)
          coeff(i,j,3) = coeffc(i,j,3)*1600*den(i,j)*zsp(i,j)* &
     &                   rbt*bav/(dpsidr*b2av) 
          coeff(i,j,4) = 1e3*coeffc(i,j,4)*rbt/(b2av* &
     &                   dpsidr)
          end do
        end do


  return 
end subroutine

