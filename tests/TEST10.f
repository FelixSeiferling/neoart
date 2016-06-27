      SUBROUTINE TEST10
C-----------------------------------------------------------------
C     CALCULATE THE HEAT TRANSPORT OF A PLASMA CONTAINING
C     OF HYDROGEN AND OXYGEN, AND COMPARE THEM WITH THE ANALYTIC 
C     EXPRESSIONS OF HIRSHMANN AND SIGMAR.
C     HERE THE HYDROGEN HEAT FLUX IS CALCULATED AS A FUNCTION 
C     OF THE IMPURITY CONTENT (GIVEN BY ALPHA) FOR DIFFERENT 
C     CHARGE STATES OF THE OXYGEN.
C-----------------------------------------------------------------
      IMPLICIT NONE

      INTEGER NS,NC,NAR,ISEL,NREG,NLEG,NENERGY,NCOF,
     +        IC,NZM,I,J,NMAXGR,K, L,ISHOT
      REAL M,T,DEN,DS,CFF4,XI,TAU,C3,ZSP,EPARR,
     +       EPS, SIGMA, NORM, ALPHA, OMTH, OMTI, RESUL
      REAL RHO, RN, E, Q, BN
      LOGICAL NEOGEO, NEOFRC

      include 'elem_config.inc'
      
      parameter(NAR = NELMAX+2)
      parameter(NZM = NIONMAX)
      PARAMETER(NMAXGR = 1000)

      DIMENSION NC(NAR),ZSP(NAR,NZM),M(NAR),T(NAR),DEN(NAR,NZM),
     +          DS(NAR,NZM,2),CFF4(NAR,NZM,4),XI(NAR,NZM),
     +          TAU(NAR,NAR),SIGMA(4),RESUL(40,10)


C     SET THE PARAMETERS FOR THE CIRCULAR GEOMETRY
      RHO = 0.05
      E = 1E-4
      Q = 2
      RN = 1.65
      BN = 2.5
C     COPY THEM INTO THE VALUES USED BY THE CODE
      CALL CIRCGEOM(1,RHO,RN,E,Q,BN)
C     USE THE CIRCULAR GEOMETRY APPROXIMATION
      ISEL = 2
C     SET THE ACCURACY
      EPS = 1E-5
C     ION-ELECTRON COLLISIONS
      NCOF = 1
C     IN THE FIRST CALL THE MATRICES HAVE TO BE CALCULATED
      NEOFRC = .FALSE.
C     RECALCULATE THE GEOMETRY COEFFICIENTS 
      NEOGEO = .TRUE. 
C     NO SHOT NUMBER 
      ISHOT = 0
C     ZERO PARALLEL ELECTRIC FIELD 
      EPARR = 0. 
C     CALCULATE THE PS CONTRIBUTION
      IC = 2
C     SET THE COUPLING 
      SIGMA(1) = 0 
      SIGMA(2) = 0
      SIGMA(3) = 1
      SIGMA(4) = 1 
C     MAKE SURE THAT ONE IS IN THE COLLISIONAL LIMIT
      T(1) = 0.02
      T(2) = 0.02

C     THE NUMBER OF SPECIES IS 2
      NS = 2
C     IONS HAVE ONLY ONE CHARGE STATE
      NC(1) = 1
      NC(2) = 1
C     THE MASS OF THE HYDROGEN AND OXYGEN
      M(1) = 1.6727E-27
      M(2) = 16*1.6727E-27
C     THE CHARGE OF THE HYDROGEN AND OXYGEN
      ZSP(1,1) = 1
C     THE DENSITY OF THE SPECIES IN 10^19 M^-3
      DEN(1,1) = 1.

C     DO THE DIFFERENT CHARGE STATES
      DO 10000 K = 1, 4
      
      ZSP(2,1) = INT(2**(K-1))
      
C     LOOP OVER ALPHA
      DO 20000 L = 1, 20

      ALPHA = 0.1*EXP(3.*(L-1)/19*LOG(10.))
      DEN(2,1) = ALPHA*DEN(1,1)/ZSP(2,1)**2

      NEOFRC = .FALSE.
C     THE THERMODYNAMIC FORCES
      DO 204 I = 1, NS
        DO 204 J = 1, NC(I)
          DS(I,J,1) = 0.
          DS(I,J,2) = 0.
 204  CONTINUE
      DS(2,1,2) = 1.  
      DS(1,1,2) = 1.


      CALL NEOART(NS,NC,NAR,NZM,ZSP,M,T,DEN,DS,RHO,EPS,
     +            ISEL,ISHOT,NREG,SIGMA,NLEG,NENERGY,NCOF,
     +            NEOGEO,NEOFRC,IC,EPARR,CFF4)


      CALL COLXI(NAR,NZM,NS,NC,ZSP,DEN,T,M,TAU,XI)


      NORM = 1.6E-22*BN**2/(2.*Q**2*T(1))
      ALPHA  = DEN(2,1)*ZSP(2,1)**2/(DEN(1,1)*ZSP(1,1)**2)


      OMTH = SQRT(2*1.6E-16*T(1)*M(1))*DEN(1,1)*1E19/TAU(1,1)
     +       /(Q*RN)
      OMTI = SQRT(2*1.6E-16*T(2)*M(2))*DEN(2,1)*1E19/TAU(2,2)
     +       /(Q*RN)

C     USE ALPHA AS X COORDINATE
      RESUL(L,1) = ALPHA

C     CODE RESULT
      RESUL(L,1+K) = -CFF4(1,1,2)*NORM/TAU(1,2)

C     RESULT OF HIRSHMANN SIGMAR
      RESUL(L,5+K) = C3(0.D0,0.D0)*C3(ALPHA,0.D0)*(ZSP(2,1)**4 + 
     +               ALPHA*ZSP(2,1))/(ALPHA*C3(0.D0,0.D0)*ZSP(2,1)**4
     +               + ALPHA*C3(ALPHA,0.D0)*SQRT(M(1)/M(2)))

20000 CONTINUE 
10000 CONTINUE

C     OPEN(11,FILE = 'TEST10.DAT')
C     THE COLLUMS HAVE THE FOLLOWING MEANING
C     COLLUM  1. : ALPHA = N_I E_I^2 / N_H E_H^2
C     COLLUM  2. : CODE RESULT FOR E_I = E
C     COLLUM  3. : CODE RESULT FOR E_I = 2 E
C     COLLUM  4. : CODE RESULT FOR E_I = 4 E
C     COLLUM  5. : CODE RESULT FOR E_I = 8 E
C     COLLUM  6. : ANALYTIC RESULT FOR E_I = E
C     COLLUM  7. : ANALYTIC RESULT FOR E_I = 2 E
C     COLLUM  8. : ANALYTIC RESULT FOR E_I = 4 E
C     COLLUM  9. : ANALYTIC REUSLT FOR E_I = 8 E
      DO 30000 I = 1, 20
        WRITE(*,1000)(RESUL(I,L),L=1,9)
30000 CONTINUE

 1000 FORMAT(14(1X,1PE13.5))

      RETURN 
      END
