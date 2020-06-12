C dseed: inital seed (dseed is automatically updated)
C nr: dimension of the double precision array r
C r: real*8 r(nr)
C GGUBS generates nr independent uniform random variables
C GGNML generates nr independent normal variables


      SUBROUTINE GGUBS_STEP(DSEED,NR,R)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   GGUBS GENERATES NR SINGLE PRECISION RANDOM VARIATES UNIFORM C
C ON (0,1) BY A LINEAR CONGRUENTIAL SCHEME.  THIS ROUTINE IS    C
C DEPENDENT ON MACHINE WORD SIZE.                               C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   LAST MODIFIED                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C   ON ENTRY                                                    C
C       DSEED   DOUBLE PRECISION                                C
C               SEED FOR GENERATOR                              C
C       NR      INTEGER                                         C
C               NUMBER OF VARIATES TO GENERATE                  C
C   ON RETURN                                                   C
C       R       REAL (NR)                                       C
C               SINGLE PRECISION ARRAY CONTAINING THE VARIATES  C
C   GGUBS CALLS                                                 C
C               DMOD                                            C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
      INTEGER            NR
      REAL*8             R(NR)
      DOUBLE PRECISION   DSEED
C
C                              LOCAL
C
      INTEGER            I
      DOUBLE PRECISION   D2P31M,D2P31,DMULTX
C
C                              MACHINE CONSTANTS
C                              D2P31M=(2**31) - 1
C                              D2P31 =(2**31)(OR AN ADJUSTED VALUE)
C
      DATA               D2P31M/2147483647.D0/
      DATA               D2P31/2147483711.D0/
      DATA               DMULTX/16807.0D+00/
C
      DO 5 I=1,NR
         DSEED=DMOD(DMULTX*DSEED,D2P31M)
         R(I) =DSEED / D2P31
  5   CONTINUE
C
C                               END OF GGUBS
C
      RETURN
      END
C===============================================================C
C===============================================================C
      SUBROUTINE GGNML(DSEED,NR,R)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   GGNML GENERATES NR SINGLE PRECISION N(0,1) RANDOM VARIATES  C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   LAST MODIFIED                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C   ON ENTRY                                                    C
C       DSEED   DOUBLE PRECISION                                C
C               SEED FOR GENERATOR                              C
C       NR      INTEGER                                         C
C               NUMBER OF VARIATES TO GENERATE                  C
C   ON RETURN                                                   C
C       R       REAL (NR)                                       C
C               SINGLE PRECISION ARRAY CONTAINING THE VARIATES  C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
      INTEGER            NR
      REAL*8             R(NR)
      DOUBLE PRECISION   DSEED
C                              LOCAL
      INTEGER             IER
C
C                              GET NR RANDOM NUMBERS
C                              UNIFORM (0,1)
C
      CALL GGUBS(DSEED,NR,R)
C
C                              TRANSFORM EACH UNIFORM DEVIATE
C
      DO 5 I=1,NR
         CALL MDNRIS(R(I),R(I),IER)
    5 CONTINUE
C
C                               END OF GGNML
C
      RETURN
      END


      SUBROUTINE MDNRIS (P,Y,IER)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                               C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   LAST MODIFIED                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL*8             P,Y
      INTEGER            IER
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL*8             EPS,G0,G1,G2,G3,H0,H1,H2,A,W,WI,SN,SD
      REAL*8             SIGMA,SQRT2,X,XINF
      DATA               XINF/1.7014E+38/
      DATA               SQRT2/1.414214/
      DATA               EPS/1.1921E-07/
      DATA               G0/.1851159E-3/,G1/-.2028152E-2/
      DATA               G2/-.1498384/,G3/.1078639E-1/
      DATA               H0/.9952975E-1/,H1/.5211733/
      DATA               H2/-.6888301E-1/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (P .GT. 0.0 .AND. P .LT. 1.0) GO TO 5
      IER = 129

      if (p .lt. 0.0D+00) then
         sigma = -1.0D+00
      else
         sigma = 1.0D+00
      endif

C      SIGMA = SIGN(1.0,P)
C      write(6,666) p,sigma
C666   format(' Sign #1: ',2f15.8)
C      pause

      Y = SIGMA * XINF
      GO TO 20
    5 IF(P.LE.EPS) GO TO 10
      X = 1.0 -(P + P)
      CALL MERFI (X,Y,IER)
      Y = -SQRT2 * Y
      GO TO 20
C                                  P TOO SMALL, COMPUTE Y DIRECTLY
   10 A = P+P
      W = SQRT(-dLOG(A+(A-A*A)))
C                                  USE A RATIONAL FUNCTION IN 1./W
      WI = 1./W
      SN = ((G3*WI+G2)*WI+G1)*WI
      SD = ((WI+H2)*WI+H1)*WI+H0
      Y = W + W*(G0+SN/SD)
      Y = -Y*SQRT2
C                               END OF MDNRIS
  20  RETURN
      END
C===============================================================C
      SUBROUTINE MERFI (P,Y,IER)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                               C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   LAST MODIFIED                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL*8             P,Y
      INTEGER            IER
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL*8             A,B,X,Z,W,WI,SN,SD,F,Z2,RINFM,A1,A2,A3,B0,B1,
     *                   B2,B3,C0,C1,C2,C3,D0,D1,D2,E0,E1,E2,E3,F0,F1,
     *                   F2,G0,G1,G2,G3,H0,H1,H2,SIGMA
      DATA               A1/-.5751703/,A2/-1.896513/,A3/-.5496261E-1/
      DATA               B0/-.1137730/,B1/-3.293474/,B2/-2.374996/
      DATA               B3/-1.187515/
      DATA               C0/-.1146666/,C1/-.1314774/,C2/-.2368201/
      DATA               C3/.5073975E-1/
      DATA               D0/-44.27977/,D1/21.98546/,D2/-7.586103/
      DATA               E0/-.5668422E-1/,E1/.3937021/,E2/-.3166501/
      DATA               E3/.6208963E-1/
      DATA               F0/-6.266786/,F1/4.666263/,F2/-2.962883/
      DATA               G0/.1851159E-3/,G1/-.2028152E-2/
      DATA               G2/-.1498384/,G3/.1078639E-1/
      DATA               H0/.9952975E-1/,H1/.5211733/
      DATA               H2/-.6888301E-1/
      DATA               RINFM/1.7014E+38/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      X = P

      if (x .lt. 0.0D+00) then
         sigma = -1.0D+00
      else
         sigma = 1.0D+00
      endif

C      SIGMA = SIGN(1.0,X)
C      write(6,666) x,sigma
C666   format(' Sign #2: ',2f15.8)
C      pause

C                                  TEST FOR INVALID ARGUMENT
      IF (.NOT.(X.GT.-1. .AND. X.LT.1.)) GO TO 30
      Z = ABS(X)
      IF (Z.LE. .85) GO TO 20
      A = 1.-Z
      B = Z
C                                  REDUCED ARGUMENT IS IN (.85,1.),
C                                     OBTAIN THE TRANSFORMED VARIABLE
    5 W = SQRT(-dLOG(A+A*B))
      IF (W.LT.2.5) GO TO 15
      IF (W.LT.4.) GO TO 10
C                                  W GREATER THAN 4., APPROX. F BY A
C                                     RATIONAL FUNCTION IN 1./W
      WI = 1./W
      SN = ((G3*WI+G2)*WI+G1)*WI
      SD = ((WI+H2)*WI+H1)*WI+H0
      F = W + W*(G0+SN/SD)
      GO TO 25
C                                  W BETWEEN 2.5 AND 4., APPROX. F
C                                     BY A RATIONAL FUNCTION IN W
   10 SN = ((E3*W+E2)*W+E1)*W
      SD = ((W+F2)*W+F1)*W+F0
      F = W + W*(E0+SN/SD)
      GO TO 25
C                                  W BETWEEN 1.13222 AND 2.5, APPROX.
C                                     F BY A RATIONAL FUNCTION IN W
   15 SN = ((C3*W+C2)*W+C1)*W
      SD = ((W+D2)*W+D1)*W+D0
      F = W + W*(C0+SN/SD)
      GO TO 25
C                                  Z BETWEEN 0. AND .85, APPROX. F
C                                     BY A RATIONAL FUNCTION IN Z
   20 Z2 = Z*Z
      F = Z+Z*(B0+A1*Z2/(B1+Z2+A2/(B2+Z2+A3/(B3+Z2))))
C                                  FORM THE SOLUTION BY MULT. F BY
C                                     THE PROPER SIGN
   25 Y = SIGMA*F
      IER = 0
      GO TO 40
C                                  ERROR EXIT. SET SOLUTION TO PLUS
C                                     (OR MINUS) INFINITY
   30 IER = 129
      Y = SIGMA * RINFM
C                               END OF MERFI
   40 RETURN
      END
C===============================================================C

      SUBROUTINE DDNOR(Y,GAUSS)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 P(6), Q(5), A(9), B(8), C(5), D(4)
      DATA P(1)/-6.58749161529837803157D-04/,
     1     P(2)/-1.60837851487422766278D-02/,
     2     P(3)/-1.25781726111229246204D-01/,
     3     P(4)/-3.60344899949804439429D-01/,
     4     P(5)/-3.05326634961232344035D-01/,
     5     P(6)/-1.63153871373020978498D-02/
      DATA Q(1)/2.33520497626869185443D-03/,
     1     Q(2)/6.05183413124413191178D-02/,
     2     Q(3)/5.27905102951428412248D-01/,
     3     Q(4)/1.87295284992346047209D00/,
     4     Q(5)/2.56852019228982242072D00/
      DATA A(1)/1.23033935479799725272D03/,
     1     A(2)/2.05107837782607146532D03/,
     2     A(3)/1.71204761263407058314D03/,
     3     A(4)/8.81952221241769090411D02/,
     4     A(5)/2.98635138197400131132D02/,
     5     A(6)/6.61191906371416294775D01/,
     6     A(7)/8.88314979438837594118D00/,
     7     A(8)/5.64188496988670089180D-01/,
     8     A(9)/2.15311535474403846343D-08/
      DATA B(1)/1.23033935480374942043D03/,
     1     B(2)/3.43936767414372163696D03/,
     2     B(3)/4.36261909014324715820D03/,
     3     B(4)/3.29079923573345962678D03/,
     4     B(5)/1.62138957456669018874D03/,
     5     B(6)/5.37181101862009857509D02/,
     6     B(7)/1.17693950891312499305D02/,
     7     B(8)/1.57449261107098347253D01/
      DATA C(1)/3.209377589138469472562D03/,
     1     C(2)/3.774852376853020208137D02/,
     2     C(3)/1.138641541510501556495D02/,
     3     C(4)/3.161123743870565596947D00/,
     4     C(5)/1.857777061846031526730D-01/
      DATA D(1)/2.844236833439170622273D03/,
     1     D(2)/1.282616526077372275645D03/,
     2     D(3)/2.440246379344441733056D02/,
     3     D(4)/2.360129095234412093499D01/
      DATA ORPI/.5641895835477562869483D0/,
     1   ROOT2/.70710678118654752440083D0/
C  THIS SUBROUTINE USES CODY'S METHOD TO EVALUATE THE CUMULATIVE
C NORMAL DISTRIBUTION. IT IS PROBABLY ACCURATE TO 19 OR 20
C SIGNIFICANT DIGITS. IT WAS WRITTEN BY JAMES MACKINNON LATE IN
C 1977, BASED ON THE CODY ARTICLE REFERRED TO IN THE DOCUMENTATION
C FOR IMSL SUBROUTINE MDNOR.
      ISW = 1
      IF (Y.LT.-16.D0) Y = -16.D0
      IF (Y.GT.16.D0) Y = 16.D0
      X = -Y*ROOT2
      IF(X.GT.0.D0) GO TO 1
      IF(X.LT.0.D0) GO TO 2
      GAUSS = .5D0
      RETURN
    2 CONTINUE
      X = - X
      ISW = -1
    1 CONTINUE
      IF(X.LT..477D0) GO TO 10
      IF(X.LE.4.D0) GO TO 20
C  EVALUATE ERFC FOR X.GT.4.0
      X2 = X*X
      XM2 = 1.D0/X2
      XM4 = XM2*XM2
      XM6 = XM4*XM2
      XM8 = XM4*XM4
      XM10 = XM6*XM4
      TOP = P(1) + P(2)*XM2 + P(3)*XM4 + P(4)*XM6 + P(5)*XM8 + P(6)*XM10
      BOT = Q(1) + Q(2)*XM2 + Q(3)*XM4 + Q(4)*XM6 + Q(5)*XM8 + XM10
      CRAP = ORPI + TOP/(BOT*X2)
      ERFC = DEXP(-X2)*CRAP/X
C
      IF(ISW.EQ.-1) ERFC = 2.D0 - ERFC
      GAUSS = ERFC*.5D0
      RETURN
   20 CONTINUE
C  EVALUATE ERFC FOR .477.LT.X.LE.4.0
      X2 = X*X
      X3 = X2*X
      X4 = X2*X2
      X5 = X3*X2
      X6 = X3*X3
      X7 = X3*X4
      X8 = X4*X4
      TOP = A(1) + A(2)*X + A(3)*X2 + A(4)*X3 + A(5)*X4 + A(6)*X5 +
     1 A(7)*X6 + A(8)*X7 + A(9)*X8
      BOT = B(1) + B(2)*X + B(3)*X2 + B(4)*X3 + B(5)*X4 + B(6)*X5 +
     1 B(7)*X6 + B(8)*X7 + X8
      ERFC = DEXP(-X2)*TOP/BOT
C
      IF(ISW.EQ.-1) ERFC = 2.D0 - ERFC
      GAUSS = ERFC*.5D0
      RETURN
   10 CONTINUE
C  EVALUATE ERF FOR X.LT..477
      X2 = X*X
      X4 = X2*X2
      X6 = X4*X2
      X8 = X4*X4
      TOP = C(1) + C(2)*X2 + C(3)*X4 + C(4)*X6 + C(5)*X8
      BOT = D(1) + D(2)*X2 + D(3)*X4 + D(4)*X6 + X8
      ERF = X*TOP/BOT
C
      ERF = ERF*ISW
      ERFC = 1.D0 - ERF
      GAUSS = ERFC*.5D0
      RETURN
      END
