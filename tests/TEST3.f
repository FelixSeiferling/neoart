
      SUBROUTINE TEST3
C-----------------------------------------------------------------
C     TEST THE REDUCED CHARGE STATE METHOD
C-----------------------------------------------------------------
      
      IMPLICIT NONE

      INTEGER NS,NC,NAR,ISEL,NREG,NLEG,NENERGY,NCOF,
     +        IC,NZM,I,J,NMAXGR, K, ISHOT
      REAL M,T,DEN,DS,CFF1,CFF2,EPS, SIGMA, ZSP, EPARR
      REAL RHO, RN, E, Q, BN, MAXA
      LOGICAL NEOGEO, NEOFRC

      include 'elem_config.inc'
      
      parameter(NAR = NELMAX+2)
      parameter(NZM = NIONMAX)
      PARAMETER(NMAXGR = 1000)


      DIMENSION NC(NAR),ZSP(NAR,NZM),M(NAR),T(NAR),DEN(NAR,NZM),
     +          DS(NAR,NZM,2),CFF1(NAR,NZM,4),CFF2(NAR,NZM,4),
     +          SIGMA(4)

C     SET THE PARAMETERS FOR THE CIRCULAR GEOMETRY
      RHO = 0.05
      E = 5E-3
      Q = 2
      RN = 1.65
      BN = 2.5
C     COPY THEM INTO THE VALUES USED BY THE CODE
      CALL CIRCGEOM(1,RHO,RN,E,Q,BN)
C     USE THE CIRCULAR GEOMETRY APPROXIMATION
      ISEL = 2
C     SET THE ACCURACY
      EPS = 1E-5
C     DO NOT FORCE A REGIME IN THE CALCULATION OF THE 
C     VISCOSITY
      NREG = 0
c     THE NUMBER OF LEGENDRE HARMONICS
      NLEG = 3
C     USE ENERGY SCATTERING IN THE CALC. OF VISCOSITY
      NENERGY = 1
C     USE ION-ELECTRON COLLISIONS
      NCOF = 1
C     IN THE FIRST CALL THE MATRICES HAVE TO BE CALCULATED
      NEOFRC = .FALSE.
C     RECALCULATE THE GEOMETRY PARAMETERS 
      NEOGEO = .TRUE. 
C     SET THE ELECTRIC FIELD TO ZERO 
      EPARR = 0. 
C     NO SHOT NUMBER GIVEN 
      ISHOT = 0 

C     USE ALL COUPLINGS IN THE PFIRSCH SCHLUETER REGIME
      SIGMA(1) = 1
      SIGMA(2) = 1
      SIGMA(3) = 1
      SIGMA(4) = 1

C     THE NUMBER OF SPECIES IS FIRST 3, TWO SPECIES HOWEVER,
C     HAVE THE SAME MASS AND TEMPERATURE. THE CALCULATIONS 
C     ARE THEN REPEATED WITH 2 SPECIES WHERE THE SECOND HAS
C     TWO CHARGE STATES. IF ALL GOES WELL THE RADIAL FLUXES
C     ARE THE SAME IN THE TWO CASES. 
      NS = 3
C     IONS AND ELECTRONS HAVE ONLY ONE CHARGE 
      NC(1) = 1
      NC(2) = 1
      NC(3) = 1
C     THE MASS OF THE PROTON AND CARBON ATOMS
      M(1) = 1.6727E-27
      M(2) = 12*1.6727E-27
      M(3) = 12*1.6727E-27
C     THE CHARGE OF THE IONS
      ZSP(1,1) = 1
      ZSP(2,1) = 2
      ZSP(3,1) = 5
      ZSP(2,2) = 5
C     THE DENSITY OF THE SPECIES IN 10^19 M^-3
      DEN(1,1) = 1.
      DEN(2,1) = 0.3
      DEN(3,1) = 0.75
      DEN(2,2) = 0.75 
C     THE TEMPERATURE IN KEV
      T(1) =  0.01
      T(2) =  0.01
      T(3) =  0.01
C     THE BANANA AND PS REGIME
      DO 10000 IC = 1, 2

C     SET THE MAXIMUM OF THE RELATIVE CHANGE TO ZERO
      MAXA = 0.

C     THE DIFFERENT THERMODYNAMIC FORCES
      DO 10001 K = 1, 6

      NS = 3
      NC(2) = 1

C     THE THERMODYNAMIC FORCES
      DO 101 I = 1, NS
        DO 101 J = 1, NC(I)
          DS(I,J,1) = 0.
          DS(I,J,2) = 0.
 101  CONTINUE
      IF (K.EQ.1) DS(1,1,1) = 1.
      IF (K.EQ.2) DS(1,1,2) = 1.
      IF (K.EQ.3) DS(2,1,1) = 1.
      IF (K.EQ.4) DS(2,1,2) = 1.
      IF (K.EQ.5) DS(3,1,1) = 1.
      IF (K.EQ.6) DS(3,1,2) = 1.

      CALL NEOART(NS,NC,NAR,NZM,ZSP,M,T,DEN,DS,RHO,EPS,
     +            ISEL,ISHOT,NREG,SIGMA,NLEG,NENERGY,NCOF,
     +            NEOGEO,NEOFRC,IC,EPARR,CFF1)

C     NOW REPEAT THE CALCULTION
      NS = 2
      NC(2) = 2

      DO 102 I = 1, NS
        DO 102 J = 1, NC(I)
          DS(I,J,1) = 0.
          DS(I,J,2) = 0.
 102  CONTINUE
      IF (K.EQ.1) DS(1,1,1) = 1.
      IF (K.EQ.2) DS(1,1,2) = 1.
      IF (K.EQ.3) DS(2,1,1) = 1.
      IF (K.EQ.4) DS(2,1,2) = 1.
      IF (K.EQ.5) DS(2,2,1) = 1.
      IF (K.EQ.6) DS(2,2,2) = 1.

      CALL NEOART(NS,NC,NAR,NZM,ZSP,M,T,DEN,DS,RHO,EPS,
     +            ISEL,ISHOT,NREG,SIGMA,NLEG,NENERGY,NCOF,
     +            NEOGEO,NEOFRC,IC,EPARR,CFF2)

      DO 103 I = 1, 6-2*IC 
        MAXA = MAX(MAXA,(CFF1(1,1,I)-CFF2(1,1,I))/CFF1(1,1,I))
        MAXA = MAX(MAXA,(CFF1(2,1,I)-CFF2(2,1,I))/CFF1(2,1,I))
        MAXA = MAX(MAXA,(CFF1(3,1,I)-CFF2(2,2,I))/CFF1(3,1,I))
 103  CONTINUE

10001 CONTINUE


      IF (IC.EQ.1) WRITE(*,1000)'BANANA REGIME'
      IF (IC.EQ.2) WRITE(*,1000)'PFIRSCH SCHLUETER REGIME'

      WRITE(*,1000)'THE MAXIMUM RELATIVE CHANGE',MAXA
 1000 FORMAT(A27,1X,1PE13.5)

10000 CONTINUE 

C THE RESULT CAN BE FOR T = 1 KEV
C              BANANA REGIME
C THE MAXIMUM RELATIVE CHANGE   8.00845E-12
C    PFIRSCH SCHLUETER REGIME
C THE MAXIMUM RELATIVE CHANGE   1.02652E-14
C GOOD AGREEMENT. HOWEVER FOR VERY LOW TEMPERATURES
C T = 1EV
C              BANANA REGIME
C THE MAXIMUM RELATIVE CHANGE   3.02278E-01
C    PFIRSCH SCHLUETER REGIME
C THE MAXIMUM RELATIVE CHANGE   4.48823E-03
C THE AGREEMENT IS STILL MORE OR LESS ACCEPTABLE FOR THE PS 
C FLUXES BUT NOT FOR THE BANANA PLATEAU FLUXES. FOR BOTH 
C REGIMES THIS LOSS OF ACCURACY CAN BE TRACED TO A LARGE 
C CONDITION NUMBER OF THE MATRICES, WHICH IS MOSTLY CONNECTED
C TO THE INCREASE IN THE FRICTION COEFFICIENTS. IT MUST BE 
C NOTED THAT THE BANANA FLUXES ARE EXTREMELY SMALL IN THIS 
C TEMPERATURE RANGE AND THE PROBLEM IS RATHER ACADEMIC. IT 
C SEEMS POSSIBLE TO CHANGE THE SOLUTION METHOD TO OVERCOME 
C THESE PROBLEMS BUT THIS IS NOT IMPLMENTED IN THIS VERSION 
C OF THE CODE
C
      RETURN 
      END
