
! C--------------------------------------------------------------------
! C     THIS ROUTINE CALCULATES THE VISCOSITY COEFFICIENTS USING THE 
! C     EXPRESSIONS OF K.C. SHAING, PHYS. PLASMAS 3 965 (1996)
! C     EXTENDED TO THE 3 LAGUERRE POLYNOMIAL EXPANSION
! C
! C     INPUT   ncm            : MAXIMUM NUMBER OF CHARGE STATES, USED
! C                              FOR CHECK ON INTERNAL CONSISTENCY
! C             NS             : THE ACTUAL NUMBER OF SPECIES
! C             js             : THE SPECIES NUMBER FOR WHICH THE 
! C                              VISCOSITY IS TO BE CALCULATED
! C             jc             : THE CHARGE STATE FOR WHICH THE VISC-
! C                              OSITY IS TO BE CALCULATED.
! C             NREG           : INTEGER THAT DETERMINES WHAT REGIME
! C                              IS TO BE FORCED IN THE CALCULATION 
! C                              OF THE VISCOSITY
! C                              0 = WEIGTED OVER ALL REGIMES
! C                              1 = THE BANANA REGIMES IS FORCED
! C                              OTHER = THE PFIRSCH SCHLUETER REGIME
! C                              IS FORCED. 
! c             BGRADP         : THE FLUX SURFACE QUANTITY N IN GRAD T
! C                              (THE INNER PRODUCT OF THE UNIT VECTOR
! C                              IN THE DIRECTION OF THE MAGNETIC FIELD
! C                              AND THE GRADIENT OF THE POLOIDAL 
! C                              ANGLE)
! C             FC             : THE NUMBER OF PASSING PARTICLES
! C             FM             : THE COEFFICIENTS NECESSARY TO CALCULATE
! C                              THE VISCOSITY IN THE PFIRSCH SCHLUETER
! C                              REGIME.
! C             MMX            : THE NUMBER OF COEFFICIENTS FM
! C             EPS            : THE ACCURACY WITH WHICH THE NUMERICAL
! C                            : COEFICIENTS SHOULD BE CALCULATED. 
! C                              NOTE THAT THE ACCURACY OF THE ANALYTIC
! C                              EQUATIONS IMPLEMENTED HERE IS NOT 
! C                              EXPECTED TO BE LARGER THAN 10% IN 
! C                              CERTAIN COLLISIONAL REGIMES. WHEN EPS = 0., 
! C                              THE VALUE 0.001 WILL BE ASSUMED.
! C             TAU            : ARRAY(ns,ns) THE COLLISION FREQUENCY
! C                              N_I M_I / TAU_IJ
! C             M              : ARRAY(ns) THE MASSES OF THE SPECIES 
! C                              IN KG
! C             T              : ARRAY(ns) TEMPERATURE OF THE SPECIES 
! C                              IN KEV
! C             XI             : ARRAY(ns,ncm) RELATIVE WEIGTH OF A 
! C                              CHARGE STATE.
! C             DEN            : ARRAY(ns,ncm) DENSITY IN UNIT 10^19 M^-3
! C             NCC            : ARRAY(ns) THAT GIVES THE NUMBER OF 
! C                              CHARGE STATES PER SPECIES.
! C     OUTPUT  MU             : ARRAY(3,3) THE VISCOSITY COEFFICIENTS
! C
! C     THE ROUTINE CALLS THE FOLLOWING SUBROUTINE
! C     VISCOL    : ROUTINE THAT CALCULATES THE ENERGY DEPENDENT 
! C                 COLLISION AND ENERGY SCATTERING FREQUENCY
! C     PERR      : ROUTINE THAT DOES THE ERROR HANDLING.
! C--------------------------------------------------------------------
subroutine viscos(js, jc, mu)       
      use init, only : ns,ncm,nc,m,t,den
      use init, only : ic, nreg, eps, nenergy
      use geometry, only : bgradp, fc , fm, mmx
      use collision, only : tau, xi
      implicit none
      
      integer :: js,jc, nv
      integer :: i,j,k,l
      logical :: nltest, nlerr

      real, dimension(3,3) :: mu, muo
      real :: TWOPI, FT, PI, X, DUM, VTH, NUD
      real :: NUE, NUT, KB, KPS, OMMN, NTOM, KTOT
      real :: VOORF, DX, WEIGHT, MEANMUO, A(6)

      NLTEST = .FALSE.

!     CALCULATE VTH
      VTH = SQRT(2.*1.6E-16*T(js)/M(js))

!    THE CONSTANT PI, 2 PI
      PI = 4.*ATAN(1.)
      TWOPI = 2.*PI

!    THE CONSTANTS TO APPROXIMATE THE NUT*IRM FUNCTION
      A(1) = 2. / 5.
      A(2) = - 22. / 105.
      A(3) = 6. / 35.
      A(4) = - 34. / 231.
      A(5) = 166. / 1287.
      A(6) = - 82. / 715.


      FT = 1.- FC

      NV = 100

      mu=0.
      muo=0.
      
!       THE LOOP OVER THE VELOCITY
        dx = 10./REAL(nv)

        do i = 0 , nv

          X = I*DX
          OMMN = VTH * X * BGRADP
        
          if (I.EQ.0) THEN

            VOORF = 0.
 
          else 

! c         NOW CALCULATE THE COLLISION AND ENERGY SCATTERING 
! C         FREQUENCY. 
	  call viscol(js, jc, x, NUD, NUE)
          NUT = 3*NUD + REAL(NENERGY)*NUE

          KB  = FT * NUD / FC
          KPS = 0.
          do j = 1, mmx
            NTOM = NUT / ( OMMN * REAL(J) )
            if (NTOM.LT.10) THEN
              DUM  = -1.5*NTOM**2 - 4.5*NTOM**4 +  &
     &             (0.25+(1.5+2.25*NTOM**2)*NTOM**2)*2*NTOM * &
     &             ATAN(1/NTOM)
            else
              NTOM = (1. / NTOM)**2 
              DUM = A(1)
              do k = 2, 6
                DUM = DUM + A(K)*NTOM**(K-1)
              end do
            end if
            KPS = KPS + FM(J)*DUM
          end do
          KPS = KPS * 1.5 * (VTH * X)**2/NUT
          if (NREG.EQ.0) THEN
            KTOT = KB*KPS/(KB+KPS)
          else
            if (NREG.EQ.1) THEN
              KTOT = KB
            else 
              KTOT = KPS
            end if
          end if
          if ((NLTEST).AND.(NV.EQ.100)) THEN
            WRITE(*,*)NUD/(X*VTH), KTOT/(X*VTH)
          end if

          VOORF = X**4 * EXP(-X**2)*KTOT*DX/3.

          end if

          if ((i.eq.0).or.(i.eq.nv)) then
             weight=1.
           else
             weight=2.*2.**mod(i,2)
          endif

          MU(1,1) = MU(1,1) + weight*VOORF
          MU(1,2) = MU(1,2) + weight*VOORF*(X**2-2.5)
          MU(2,2) = MU(2,2) + weight*VOORF*(X**2-2.5)**2
          MU(1,3) = MU(1,3) + weight*VOORF*(35./8.-3.5*X**2 + &
     &              0.5*X**4)
          MU(2,3) = MU(2,3) + weight*VOORF*(35./8.-3.5*X**2 + &
     &              0.5*X**4)*(X**2-2.5)
          MU(3,3) = MU(3,3) + weight*VOORF*(35./8.-3.5*X**2 + &
     &              0.5*X**4)**2
                      
    end do

        VOORF = 1.E19*DEN(js,jc)*M(js)* 8. / (3.*SQRT(PI))
        MU(2,1) = MU(1,2)
        MU(3,1) = MU(1,3)
        MU(3,2) = MU(2,3)
        MEANMUO = 0.
        do k =1,3; do l=1,3
          MU(K,L) = MU(K,L) * VOORF
          MEANMUO = MEANMUO + MU(K,L)
	end do; end do
        MEANMUO = MEANMUO / 9.
        NLERR = .FALSE.
        do k =1,3; do l=1,3
         if (ABS(MU(K,L)-MUO(K,L)).GT.ABS(EPS*MU(K,L))) THEN 
            if (ABS(MU(K,L)).GT.EPS*MEANMUO) NLERR = .TRUE.
          end if
	end do; end do
        if (NLERR) THEN
        do k =1,3; do l=1,3
            MUO(K,L) = MU(K,L)
	end do; end do
          NV = NV*2
          if (NV.GT.4e6) CALL PERR(4)
        end if
         

    return
end