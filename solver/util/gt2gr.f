      PROGRAM GT2GR
*
*     derived from fvec.f on 1998/11/29
*     modified on 2001/05/18 for Q
*     This routine transforms Gt3 output from a linear run
*     to GrADS file same as in the response.grd.
*
*     variables: PSI, XI, U, V, T, Z, Ps, Q, dT, dQ
*
      include 'dim.f'
*
*
      REAL*4       XMISS
      PARAMETER ( XMISS = -999. )
*
      REAL*4       DAT1( IMAX*JMAX, KMAX)
      REAL*4       DAT2( IMAX, JMAX     )
      REAL*4       DAT3( IMAX*JMAX, KMAX)
      REAL*4       G   ( IMAX, JMAX     )  
      REAL*4       GPS ( IMAX*JMAX      )
      CHARACTER    HEAD( 64 ) * 16
*
*     [work]
      REAL*4      P( KMAX )          !! pressure level
      REAL*4      S( KMAX )          !! sigma level
      REAL*8      SIG( KMAX )        !! sigma level
      REAL*8      SIGM( KMAX+1 )     !! sigma level increment
      REAL*8      DSIG( KMAX )       !! unused
      REAL*8      DSIGM( KMAX+1 )    !! unused
      REAL*8      FACT( 11 )         !! factor
      CHARACTER   HSIG *(16)         !! unused
      CHARACTER   HSIGM *(16)        !! unused
      INTEGER     IFS, IFC, IFU, IFV, IFW, IFT, IFZ
      INTEGER     IFP, IFQ, IFX, IFY, IFO, IFB
      INTEGER     IJ, I, J, K
      INTEGER     NOUT
      DATA        NOUT / 6 /
*
*     [intrinsic]
*
      NAMELIST    /NMFGT/   CFS, CFC, CFU, CFV, CFW,
     &                      CFT, CFZ, CFP, CFQ, CFX, 
     &                      CFY, CFO, FACT, OPL
      NAMELIST    /NMBS/    CBS0, CBS
      NAMELIST    /NMCLS/   OCLASSIC
*
      CHARACTER*900 CFS         !! input file name (gtool3)
      CHARACTER*900 CFC         !! input file name (gtool3)
      CHARACTER*900 CFU         !! input file name (gtool3)
      CHARACTER*900 CFV         !! input file name (gtool3)
      CHARACTER*900 CFW         !! input file name (gtool3)
      CHARACTER*900 CFT         !! input file name (gtool3)
      CHARACTER*900 CFZ         !! input file name (gtool3)
      CHARACTER*900 CFP         !! input file name (gtool3)
      CHARACTER*900 CFQ         !! input file name (gtool3)
      CHARACTER*900 CFX         !! input file name (gtool3)
      CHARACTER*900 CFY         !! input file name (gtool3)
      CHARACTER*900 CFO         !! output file name (grads)
      CHARACTER*900 CBS0        !! basic state (dummy)
      CHARACTER*900 CBS         !! basic state
      LOGICAL      OPL          !! convert to p-level?
      LOGICAL      OCLASSIC     !! conventional?

      DATA         FACT      / 11*1.D0 /
      DATA         OPL       / .TRUE. /
      DATA         OCLASSIC  / .TRUE. /
*
*     open the NAMELIST file
*
      REWIND( 1 )
      OPEN ( 0, FILE='IEEE_ERROR')
      OPEN ( 1, FILE = 'SETPAR', STATUS = 'OLD')
*
      READ( 1, NMFGT) 
      REWIND( 1 )
      READ( 1, NMBS) 
      REWIND( 1 )
      READ( 1, NMCLS) 
*
*
      WRITE( NOUT, * ) '### start Gtool3 ---> GrADS ###'
      WRITE( NOUT, * )
*------------------------------------------
*     hujun, 20130701
      write(NOUT,*) 'read basic state start!'
*------------------------------------------
*
*     read basic state Ps
*
      IFB = 55
      OPEN( IFB, FILE = CBS, FORM = 'UNFORMATTED', STATUS = 'UNKNOWN' )
*------------------------------------------
*     hujun, 20130701
      write(NOUT,*) 'open CBS done!'
*------------------------------------------
      DO 500 I = 1, KMAX*3
         READ( IFB ) G
  500 CONTINUE
*------------------------------------------
*     hujun, 20130701
      write(NOUT,*) 'read G first time!'
*------------------------------------------
      READ( IFB ) G
*------------------------------------------
*     hujun, 20130701
      write(NOUT,*) 'read G second time!'
*------------------------------------------
      DO 510 J = 1, JMAX
         DO 520 I = 1, IMAX
            IJ = IJ + 1
            GPS( IJ ) = DBLE( G( I, J) )
  520    CONTINUE
  510 CONTINUE
      CLOSE( IFB )
*------------------------------------------
*     hujun, 20130701
      write(NOUT,*) 'read basic state done!'
*------------------------------------------
*
*     sigma coodinate
*
      CALL SETCONS
      CALL SETSIG
     O     ( SIG   , DSIG , HSIG  , 
     O       SIGM  , DSIGM, HSIGM  )
      DO 530 K = 1, KMAX
         S( K ) = SNGL( SIG( K ) )
         P( K ) = SNGL( PLEV( K ) )
  530 CONTINUE
*------------------------------------------
*     hujun, 20130701
      write(NOUT,*) 'sigma coodinate done!'
*------------------------------------------
*
*     open files
*
      IFS = 10
      IFC = 11
      IFU = 12
      IFV = 13
      IFW = 14
      IFT = 15
      IFZ = 16
      IFP = 17
      IFQ = 18
      IFX = 19
      IFY = 20
      IFO = 98
      OPEN ( IFS, FILE = CFS, FORM = 'UNFORMATTED',
     $     STATUS = 'OLD' )
      OPEN ( IFC, FILE = CFC, FORM = 'UNFORMATTED',
     $     STATUS = 'OLD' )
      OPEN ( IFU, FILE = CFU, FORM = 'UNFORMATTED',
     $     STATUS = 'OLD' )
      OPEN ( IFV, FILE = CFV, FORM = 'UNFORMATTED',
     $     STATUS = 'OLD' )
      OPEN ( IFW, FILE = CFW, FORM = 'UNFORMATTED',
     $     STATUS = 'OLD' )
      OPEN ( IFT, FILE = CFT, FORM = 'UNFORMATTED',
     $     STATUS = 'OLD' )
      OPEN ( IFZ, FILE = CFZ, FORM = 'UNFORMATTED',
     $     STATUS = 'OLD' )
      OPEN ( IFP, FILE = CFP, FORM = 'UNFORMATTED',
     $     STATUS = 'OLD' )
      IF( .NOT. OCLASSIC ) THEN
         OPEN ( IFQ, FILE = CFQ, FORM = 'UNFORMATTED',
     $        STATUS = 'OLD' )
         OPEN ( IFX, FILE = CFX, FORM = 'UNFORMATTED',
     $        STATUS = 'OLD' )
         OPEN ( IFY, FILE = CFY, FORM = 'UNFORMATTED',
     $        STATUS = 'OLD' )
      ENDIF
      OPEN ( IFO, FILE = CFO, FORM = 'UNFORMATTED', STATUS = 'UNKNOWN' )
      
      WRITE( NOUT, * ) 'Input data :', CFS
      WRITE( NOUT, * ) 'Input data :', CFC
      WRITE( NOUT, * ) 'Input data :', CFU
      WRITE( NOUT, * ) 'Input data :', CFV
      WRITE( NOUT, * ) 'Input data :', CFW
      WRITE( NOUT, * ) 'Input data :', CFT
      WRITE( NOUT, * ) 'Input data :', CFZ
      WRITE( NOUT, * ) 'Input data :', CFP
      IF( .NOT. OCLASSIC ) THEN
         WRITE( NOUT, * ) 'Input data :', CFQ
         WRITE( NOUT, * ) 'Input data :', CFX
         WRITE( NOUT, * ) 'Input data :', CFY
      ENDIF
      WRITE( NOUT, * ) '................'
      WRITE( NOUT, * ) 'Output data:', CFO
      WRITE( NOUT, * ) '................'
      WRITE( NOUT, * ) 
*------------------------------------------
*     hujun, 20130701
      write(NOUT,*) 'open data done!'
*------------------------------------------
*
*     read GrADS data
*
      WRITE( NOUT, * ) 'Read data'
      WRITE( NOUT, * ) '................'
      WRITE( NOUT, * ) 
*
 10   READ( IFS, END = 99 ) HEAD
      WRITE( NOUT, *) 'DATE:',HEAD( 27 )
      READ( IFS ) DAT1          !! psi
*------------------------------------------
*     hujun, 20130701
      write(NOUT,*) 'read data done!'
*------------------------------------------
*
*     sigma to pressure
* 
      IF( OPL ) THEN
         CALL S2P
     O        ( DAT3,
     I        DAT1 , GPS, P, S, IMAX*JMAX, KMAX, XMISS )
      ELSE
         CALL COPY4 
     O        ( DAT3,
     I        DAT1, IMAX*JMAX*KMAX )
      ENDIF
*
      DO 20 K = 1, KMAX
         IJ = 0
         DO 30 J = 1, JMAX
            DO 40 I = 1, IMAX
               IJ = IJ + 1
               G( I, J) = FACT( 1 ) * DAT3( IJ, K)
               IF( DAT3( IJ, K) .EQ. XMISS ) G( I, J) = XMISS
   40       CONTINUE
   30    CONTINUE
         WRITE( IFO ) G
   20 CONTINUE
*
      READ( IFC ) HEAD
      READ( IFC ) DAT1          !! xi
*
*     sigma to pressure
* 
      IF( OPL ) THEN
         CALL S2P
     O        ( DAT3,
     I        DAT1 , GPS, P, S, IMAX*JMAX, KMAX, XMISS )
      ELSE
         CALL COPY4 
     O        ( DAT3,
     I        DAT1, IMAX*JMAX*KMAX )
      ENDIF
*
      DO 50 K = 1, KMAX
         IJ = 0
         DO 60 J = 1, JMAX
            DO 70 I = 1, IMAX
               IJ = IJ + 1
               G( I, J) = FACT( 2 ) * DAT3( IJ, K)
               IF( DAT3( IJ, K) .EQ. XMISS ) G( I, J) = XMISS
   70       CONTINUE
   60    CONTINUE
         WRITE( IFO ) G
   50 CONTINUE
*
      READ( IFU ) HEAD
      READ( IFU ) DAT1          !! u
*
*     sigma to pressure
* 
      IF( OPL ) THEN
         CALL S2P
     O        ( DAT3,
     I        DAT1 , GPS, P, S, IMAX*JMAX, KMAX, XMISS )
      ELSE
         CALL COPY4 
     O        ( DAT3,
     I        DAT1, IMAX*JMAX*KMAX )
      ENDIF
*
      DO 80 K = 1, KMAX
         IJ = 0
         DO 90 J = 1, JMAX
            DO 100 I = 1, IMAX
               IJ = IJ + 1
               G( I, J) = FACT( 3 ) * DAT3( IJ, K)
               IF( DAT3( IJ, K) .EQ. XMISS ) G( I, J) = XMISS
  100       CONTINUE
   90    CONTINUE
         WRITE( IFO ) G
   80 CONTINUE
*
      READ( IFV ) HEAD
      READ( IFV ) DAT1          !! v
*
*     sigma to pressure
* 
      IF( OPL ) THEN
         CALL S2P
     O        ( DAT3,
     I        DAT1 , GPS, P, S, IMAX*JMAX, KMAX, XMISS )
      ELSE
         CALL COPY4 
     O        ( DAT3,
     I        DAT1, IMAX*JMAX*KMAX )
      ENDIF
*
      DO 110 K = 1, KMAX
         IJ = 0
         DO 120 J = 1, JMAX
            DO 130 I = 1, IMAX
               IJ = IJ + 1
               G( I, J) = FACT( 4 ) * DAT3( IJ, K)
               IF( DAT3( IJ, K) .EQ. XMISS ) G( I, J) = XMISS
  130       CONTINUE
  120    CONTINUE
         WRITE( IFO ) G
  110 CONTINUE
*
      READ( IFW ) HEAD
      READ( IFW ) DAT1          !! omega
*
*     sigma to pressure
* 
      IF( OPL ) THEN
         CALL S2P
     O        ( DAT3,
     I        DAT1 , GPS, P, S, IMAX*JMAX, KMAX, XMISS )
      ELSE
         CALL COPY4 
     O        ( DAT3,
     I        DAT1, IMAX*JMAX*KMAX )
      ENDIF
*
      DO 140 K = 1, KMAX
         IJ = 0
         DO 150 J = 1, JMAX
            DO 160 I = 1, IMAX
               IJ = IJ + 1
               G( I, J) = FACT( 5 ) * DAT3( IJ, K)
               IF( DAT3( IJ, K) .EQ. XMISS ) G( I, J) = XMISS
  160       CONTINUE
  150     CONTINUE
         WRITE( IFO ) G
  140 CONTINUE
*
      READ( IFT ) HEAD
      READ( IFT ) DAT1          !! T
*
*     sigma to pressure
* 
      IF( OPL ) THEN
         CALL S2P
     O        ( DAT3,
     I        DAT1 , GPS, P, S, IMAX*JMAX, KMAX, XMISS )
      ELSE
         CALL COPY4 
     O        ( DAT3,
     I        DAT1, IMAX*JMAX*KMAX )
      ENDIF
*
      DO 170 K = 1, KMAX
         IJ = 0
         DO 180 J = 1, JMAX
            DO 190 I = 1, IMAX
               IJ = IJ + 1
               G( I, J) = FACT( 6 ) * DAT3( IJ, K)
               IF( DAT3( IJ, K) .EQ. XMISS ) G( I, J) = XMISS
  190       CONTINUE
  180    CONTINUE
         WRITE( IFO ) G
  170 CONTINUE
*
      READ( IFZ ) HEAD
      READ( IFZ ) DAT1          !! z 
*
*     sigma to pressure
* 
      IF( OPL ) THEN
         CALL S2P
     O        ( DAT3,
     I        DAT1 , GPS, P, S, IMAX*JMAX, KMAX, XMISS )
      ELSE
         CALL COPY4 
     O        ( DAT3,
     I        DAT1, IMAX*JMAX*KMAX )
      ENDIF
*
      DO 200 K = 1, KMAX
         IJ = 0
         DO 210 J = 1, JMAX
            DO 220 I = 1, IMAX
               IJ = IJ + 1
               G( I, J) = FACT( 7 ) * DAT3( IJ, K)
               IF( DAT3( IJ, K) .EQ. XMISS ) G( I, J) = XMISS
  220       CONTINUE
  210    CONTINUE
         WRITE( IFO ) G
  200 CONTINUE
*
      READ( IFP ) HEAD
      READ( IFP ) DAT2          !! Ps
      DO 230 J = 1, JMAX
         DO 240 I = 1, IMAX
            G( I, J) = FACT( 8 ) * DAT2( I, J)
            IF( DAT2( I, J) .EQ. XMISS ) G( I, J) = XMISS
  240    CONTINUE
  230 CONTINUE
      WRITE( IFO ) G
*     
      IF( OCLASSIC ) GOTO 10
*
      READ( IFQ ) HEAD
      READ( IFQ ) DAT1          !! q
*
*     sigma to pressure
* 
      IF( OPL ) THEN
         CALL S2P
     O        ( DAT3,
     I        DAT1 , GPS, P, S, IMAX*JMAX, KMAX, XMISS )
      ELSE
         CALL COPY4 
     O        ( DAT3,
     I        DAT1, IMAX*JMAX*KMAX )
      ENDIF
*
      DO 300 K = 1, KMAX
         IJ = 0
         DO 310 J = 1, JMAX
            DO 320 I = 1, IMAX
               IJ = IJ + 1
               G( I, J) = FACT( 9 ) * DAT3( IJ, K)
               IF( DAT3( IJ, K) .EQ. XMISS ) G( I, J) = XMISS
  320       CONTINUE
  310    CONTINUE
         WRITE( IFO ) G
  300 CONTINUE
*
      READ( IFX ) HEAD
      READ( IFX ) DAT1          !! dT
*
*     sigma to pressure
* 
      IF( OPL ) THEN
         CALL S2P
     O        ( DAT3,
     I        DAT1 , GPS, P, S, IMAX*JMAX, KMAX, XMISS )
      ELSE
         CALL COPY4 
     O        ( DAT3,
     I        DAT1, IMAX*JMAX*KMAX )
      ENDIF
*
      DO 330 K = 1, KMAX
         IJ = 0
         DO 340 J = 1, JMAX
            DO 350 I = 1, IMAX
               IJ = IJ + 1
               G( I, J) = FACT( 10 ) * DAT3( IJ, K)
               IF( DAT3( IJ, K) .EQ. XMISS ) G( I, J) = XMISS
  350       CONTINUE
  340    CONTINUE
         WRITE( IFO ) G
  330 CONTINUE
*
      READ( IFY ) HEAD
      READ( IFY ) DAT1          !! dQ
*
*     sigma to pressure
* 
      IF( OPL ) THEN
         CALL S2P
     O        ( DAT3,
     I        DAT1 , GPS, P, S, IMAX*JMAX, KMAX, XMISS )
      ELSE
         CALL COPY4 
     O        ( DAT3,
     I        DAT1, IMAX*JMAX*KMAX )
      ENDIF
*
      DO 360 K = 1, KMAX
         IJ = 0
         DO 370 J = 1, JMAX
            DO 380 I = 1, IMAX
               IJ = IJ + 1
               G( I, J) = FACT( 11 ) * DAT3( IJ, K)
               IF( DAT3( IJ, K) .EQ. XMISS ) G( I, J) = XMISS
  380       CONTINUE
  370    CONTINUE
         WRITE( IFO ) G
  360 CONTINUE
*
      GOTO 10
*     
   99 CLOSE( IFS )
      CLOSE( IFC )
      CLOSE( IFU )
      CLOSE( IFV )
      CLOSE( IFW )
      CLOSE( IFT )
      CLOSE( IFZ )
      CLOSE( IFP )
      IF( .NOT. OCLASSIC ) THEN
         CLOSE( IFQ )
         CLOSE( IFX )
         CLOSE( IFY )
      ENDIF
      CLOSE( IFO )
      WRITE( NOUT, * )
      WRITE( NOUT, * ) '### finish Gtool3 ---> GrADS ###'
*
*      
      STOP
      END
*##############################################
      SUBROUTINE S2P
     O         ( GDO  ,
     I           GDI, GDPS, PLI, SLI, IJDIM, KMAX, XMISS )
*
*   [PARAM] 
      INTEGER    IJDIM
      INTEGER    KMAX
*   [OUTPUT] 
      REAL       GDO ( IJDIM, KMAX )
*   [INPUT] 
      REAL       GDPS( IJDIM )
      REAL       GDI ( IJDIM, KMAX )
      REAL       PLI ( KMAX )
      REAL       SLI ( KMAX )
      LOGICAL    OMISB
      LOGICAL    OMIST
      REAL       XMISS
      DATA       OMISB / .FALSE. /
      DATA       OMIST / .FALSE. /
*
*   [INTERNAL WORK] 
      INTEGER    KMAXD
      PARAMETER  ( KMAXD = 100 )
      REAL       PILN ( KMAXD )
      REAL       POLN ( KMAXD )
      REAL       GDIZ ( KMAXD )
      REAL       GDOZ ( KMAXD )
      INTEGER    IJ, KI, KO
      REAL       TMP
*
      DO 2000 IJ = 1, IJDIM
         DO 2010 KI = 1, KMAX
            GDIZ(KI) = GDI(IJ,KI)
            PILN(KI) = LOG( GDPS(IJ)*SLI(KI) )
 2010    CONTINUE 
         IF(PILN(1).LT.PILN(2))THEN
            DO 2020 KI = 1, KMAX/2
               TMP              = GDIZ(KMAX+1-KI)
               GDIZ(KMAX+1-KI) = GDIZ(KI)
               GDIZ(KI)         = TMP
               TMP              = PILN(KMAX+1-KI)
               PILN(KMAX+1-KI) = PILN(KI)
               PILN(KI)         = TMP
 2020       CONTINUE 
         ENDIF
         DO 2030 KO = 1, KMAX
            POLN(KO) = LOG( PLI(KO) )
 2030    CONTINUE
         CALL SPLINE
     O        ( GDOZ,
     I                 POLN,  KMAX,
     I          GDIZ,  PILN,  KMAX,
     I          OMISB, OMIST, XMISS  )
         DO 2100 KO = 1, KMAX
            GDO(IJ,KO) = GDOZ(KO)
            IF( GDPS(IJ) .LT. PLI(KO) ) GDO(IJ,KO) = XMISS
 2100    CONTINUE 
 2000 CONTINUE
*
      RETURN
      END
*********************************************************************
      SUBROUTINE SPLINE
     O         ( ZI   ,
     I                  XI   , LMAX ,
     I           Z    , X    , KMAX , 
     I           OMISB, OMIST, XMISS  )
*
      INTEGER    LMAX, KMAX
      REAL       ZI  ( LMAX )
      REAL       XI  ( LMAX )
      REAL       Z   ( KMAX )
      REAL       X   ( KMAX )
      LOGICAL    OMISB
      LOGICAL    OMIST
      REAL       XMISS

*   [INTERNAL WORK] 
      INTEGER    KMAXD
      PARAMETER  ( KMAXD = 1024 )
      REAL       Y2  ( KMAXD )
      REAL       U   ( KMAXD )
      REAL       SIG, P, QN, UN, H, A, B 
      INTEGER    K, L, KU
*
      IF ( KMAX .GT. KMAXD ) THEN
         WRITE (6,*) ' ### KMAXD IS TOO SMALL ', KMAXD, KMAX
         STOP
      ENDIF
*
      Y2(1)=0.
      U (1)=0.
      DO 120 K=2, KMAX-1
         SIG   = (X(K)-X(K-1))/(X(K+1)-X(K-1))
         P     = SIG*Y2(K-1)+2.
         Y2(K) = (SIG-1.)/P
         U (K) = (6.*( (Z(K+1)-Z(K))/(X(K+1)-X(K))
     &                -(Z(K)-Z(K-1))/(X(K)-X(K-1)))
     &              /(X(K+1)-X(K-1))
     &             - SIG*U(K-1)                     )/P
  120 CONTINUE
      QN = 0.
      UN = 0.
      Y2(KMAX) = (UN-QN*U(KMAX-1))/(QN*Y2(KMAX-1)+1.)
      DO 130 K= KMAX-1, 1, -1
         Y2(K) = Y2(K)*Y2(K+1)+U(K)
  130 CONTINUE
*
      DO 500 L = 1, LMAX
         KU = 1
         DO 300 K = 1, KMAX
            IF( X(K) .LT. XI(L) ) THEN
               KU = K
               GOTO 310
            ENDIF
  300    CONTINUE
         KU = KMAX+1
  310    CONTINUE
*
         IF      ( KU .EQ. 1 ) THEN
            IF ( OMISB ) THEN
               ZI(L) = XMISS
            ELSE
               ZI(L) = Z(1)
            ENDIF
         ELSE IF ( KU .EQ. KMAX+1 ) THEN
            IF ( OMIST ) THEN
               ZI(L) = XMISS
            ELSE
               ZI(L) = Z(KMAX)
            ENDIF            
         ELSE
            KU   = MAX(KU,2)
            H    = X(KU)-X(KU-1)
            A    = (X(KU)-XI(L))/H
            B    = (XI(L)-X(KU-1))/H
            ZI(L)= A*Z(KU-1)+B*Z(KU)
     &           + (A*(A*A-1)*Y2(KU-1)+B*(B*B-1)*Y2(KU))*(H*H)/6.
         ENDIF
  500 CONTINUE
*
      RETURN
      END
***********************************************************************
      SUBROUTINE COPY4      !! copy matrix
     O         ( DATAO  ,
     I           DATAI  ,
     D           IDIM     )
*
*   [PARAM]
      INTEGER    IDIM
*
*   [OUTPUT]       
      REAL*4     DATAO ( IDIM ) !! output data
*
*   [INPUT] 
      REAL*4     DATAI ( IDIM ) !! input data
*
*   [INTERNAL WORK] 
      INTEGER    I
*
      DO 1100 I = 1, IDIM
        DATAO ( I ) = DATAI ( I )
 1100 CONTINUE
*
      RETURN
      END
