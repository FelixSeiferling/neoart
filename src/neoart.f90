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
subroutine neoart(neogeo,neofrc,ic,coeff)
     
     
     use init, only :  ns,nc, ncm, zsp, m, t, den, ds,    &
                       & ds, rho ,eps, ISEL, ISHOT, nreg, & 
                       & sigma, nleg, nenergy, ncof, eparr
      use collision, only : tau, xi
      use collision, only : colxi
      use init, only : neo_abort
      IMPLICIT NONE
      include 'elem_config90.inc'

      integer, intent(in) :: ic
      logical, intent(in) :: neogeo, neofrc
      real, intent(out) :: coeff
      integer :: NMAXGR, NKEEP
      PARAMETER(NKEEP = nionmax)
      PARAMETER(NMAXGR = 1000)
      

  

      integer :: i,j,k,l,ikeep, ierr, mmx
      integer :: ISK 
      real :: EPARN,R2I,BAV,B2AV,BGRADP,FC,FM,RHOK, &
     &       BI2A,RBT,DPSIDR,RNQ,UAI, DUM, MM, NN, &
     &       TA,TB,MA,MB,GCLASS,COEFFC
      DIMENSION :: GCLASS(NKEEP), COEFFC(ns,ncm,4),  &
     &          COEFF(ns,ncm,4),BAV(NKEEP),B2AV(NKEEP), &
     &          BGRADP(NKEEP),FC(NKEEP),RHOK(NKEEP),FM(NKEEP,NMAXGR), &
     &          MMX(NKEEP),R2I(NKEEP), ISK(NKEEP),BI2A(NKEEP), &
     &          RBT(NKEEP),DPSIDR(NKEEP),RNQ(NKEEP),UAI(ns,ncm,3), &
     &          MM(3,3),NN(3,3)

     real, save, allocatable :: la(:,:,:), lab(:,:,:,:), deni(:,:)
     integer, save :: IT

!     STORE THE GEOMETRY DEPENDENT VARIABLES
      COMMON / KGEOM / BAV, B2AV, BGRADP, FC, FM, MMX, RHOK, &
     &         BI2A, RBT, DPSIDR, RNQ, ISK

      DATA IKEEP / 0  / 
      
   if(.not. allocated(la)) then
      allocate(la(3,3,ns),stat=ierr)
      if (ierr /= 0) call neo_abort('could not allocate array')
      la= .0
   end if 
   if(.not. allocated(lab)) then
      allocate(lab(3,3,ns,ns),stat=ierr)
      if (ierr /= 0) call neo_abort('could not allocate array')
      lab= .0
   end if 
   if(.not. allocated(deni)) then
      allocate(deni(ns,ncm),stat=ierr)
      if (ierr /= 0) call neo_abort('could not allocate array')
      deni= .0
   end if 

      
      if  ((NEOFRC).AND.(NEOGEO)) CALL PERR(16)
 
!     if THE FRICTION AND VISCOSITY COEFFICIENTS ARE NEWLY CAL-
!     CULATED THEN GO THROUGH THE FOLLOWING LOOP. 
       if (.NOT.NEOFRC) THEN
! 
! !     STORE THE DENSITY VALUES
        do i = 1, ns
         do j=1,nc(i)
          DENI(I,J) = DEN(I,J)
         end do
       end do
      
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


!      IS THE LOGICAL NEOGEO SET ????
       if (NEOGEO) THEN
        IKEEP = 1
        IT = IKEEP
        RHOK(IT) = RHO
        ISK(IT) = ISEL
        CALL GEOM(NMAXGR,BAV(IT),B2AV(IT),BI2A(IT), &
     &    RBT(IT),BGRADP(IT),DPSIDR(IT),RNQ(IT),FC(IT),GCLASS(IT), &
     &    FM(IT,1),MMX(IT),R2I(IT))
       else
!       DOES THE VALUE OF RHO EXIST
        IT = 0
        do i = 1, IKEEP
          if ((RHO.EQ.RHOK(I)).AND.(ISEL.EQ.ISK(I)))  IT = I
        end do 
        if (IT.EQ.0) THEN
          IKEEP = IKEEP + 1
          if (IKEEP.GT.NKEEP) THEN
            IKEEP = NKEEP
            CALL PERR(17)
          end if
          IT = IKEEP
          RHOK(IT) = RHO
          ISK(IT) = ISEL
          CALL GEOM(NMAXGR,BAV(IT),B2AV(IT), &
     &              BI2A(IT),RBT(IT),BGRADP(IT),DPSIDR(IT),RNQ(IT), &
     &              FC(IT),GCLASS(IT),FM(IT,1),MMX(IT),R2I(IT))
        end if
      end if

      CALL COLXI(ns,nc,ncm,ZSP,DEN,T,M)

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
      do i = 1, ns
        do l = 1, 3
          do k = 1, 3
            LA(k,l,i) = 0.
          end do
        end do
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

!     END OF THE IF STATEMENT ON NEOFRC
      else
!      RESTORE THE DENSITY VALUES
       do i = 1, ns
         do j=1,nc(i)
          DEN(I,J) = DENI(I,J)
         end do
       end do
        end if      
! 
! 
!     MAKE THE COEFFICIENTS ZERO
      do i = 1, ns
        do j = 1, nc(i)
          do k = 1, 4
            COEFF(I,J,K) = 0.
            COEFFC(I,J,K) = 0.
          end do
        end do
      end do


!     NORMALIZE THE EPARR
      EPARN = 1E-3*DPSIDR(IT)*R2I(IT)*EPARR 
      
      
!     write(*,*) lab
!    write(*,*) la
!     CALCULATE THE CLASSICAL TRANSPORT CONTRIBUTION
      if ((IC.EQ.0).OR.(IC.EQ.3)) THEN
!      call CLASS(ns,ncm,ns,nc,T,ZSP,XI,LA,LAB,DS,COEFFC)
!       ADD THE CORRECT NORMALIZATION FACTOR, PUT RESULT 
!       IN COEFF AND CLEAR COEFFC
        do i = 1, ns
          do j = 1, nc(i)
            do k = 1, 2
              COEFF(I,J,K) = COEFF(I,J,K)+ COEFFC(I,J,K)* &
     &          GCLASS(IT)/(1.6022E-22*ZSP(I,J)*DPSIDR(IT)**2) 
              COEFFC(I,J,K) = 0.
            end do
          end do
        end do
      end if

!     CALCULATE THE PFIRSCH SCHLUETER CONTRIBUTION
      if ((IC.EQ.2).OR.(IC.EQ.3)) THEN
!        call PS(ns,ncm,nc,NEOFRC,tau,XI,LA,LAB,m,t,DENI, &
!     &           ZSP,sigma,ds,RNQ(it),uai,coeffc)
!       ADD THE CORRECT NORMALIZATION FACTOR, PUT RESULT IN 
!       COEFF AND CLEAR COEFFC
        do i = 1, ns
          do j = 1, nc(i)
            do k = 1, 2
              coeff(i,j,k) = coeff(i,j,k) - coeffc(i,j,k)* &
     &          (1.- B2AV(it)*BI2A(it))*(RBT(it)/DPSIDR(it))**2/ &
     &          (1.60E-22*ZSP(I,J)*b2av(it))
              coeffc(i,j,k) = 0.
            end do
          end do
        end do
      end if


      if ((IC.EQ.1).OR.(IC.EQ.3)) then
!          call bp(ns,ncm,nc,neofrc,t,m,den,zsp,la, &
!       &   lab,eps,nenergy,nleg,nreg,bgradp,fc,fm,mmx,ds,eparn, &
!       &   coeffc)
!       ADD THE CORRECT NORMALIZING FACTOR. PUT RESULT IN 
!       COEFF AND CLEAR COEFFC
        do i = 1, ns
          do j = 1, nc(i)
            do k = 1, 2
              coeff(i,j,k) = coeff(i,j,k) - coeffc(i,j,k)* &
     &          (RBT(it)/DPSIDR(it))**2/(1.60E-22*ZSP(I,J)* &
     &           b2av(it))
              coeffc(i,j,k) = 0.
            end do
          end do
        end do
      end if 

        do i = 1, ns
          do j = 1, nc(i)
          coeff(i,j,3) = coeffc(i,j,3)*1600*den(i,j)*zsp(i,j)* &
     &                   rbt(it)*bav(it)/(dpsidr(it)*b2av(it)) 
          coeff(i,j,4) = 1e3*coeffc(i,j,4)*rbt(it)/(b2av(it)* &
     &                   dpsidr(it))
          end do
        end do
        
        
!   deallocate(la)
!   deallocate(lab)
!   deallocate(deni)
     
  return
end subroutine

