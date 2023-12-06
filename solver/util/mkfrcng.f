      PROGRAM MKFRCNG
*
*     derived from out_agcm/donly/mktfrc.f on 1998/11/08
*     This routine makes simplified forcing function
*     for the 3-D linear response.
*
*     Alternative selection of the forcing shape is:
*      Horizontal: Elliptic function/zonally uniform function 
*      Vertical  : Sinusoidal function/Gamma function/uniform
*
      include 'dim.f'
C     2022.11.11: Ziming Chen 
      include '/global/homes/c/chenzm/intel/netcdf_intel/include/netcdf.
     *inc'
C     2022.11.11: Ziming Chen 
*
*
      INTEGER      MAXN
      PARAMETER (  MAXN = 2*NMAX*(NVAR*KMAX+1) )
*
      REAL*8       GFRCT( IDIM*JMAX, KMAX)  
      REAL*8       WFRCT( NMDIM, KMAX)  
      REAL*8       WFRCF( NMDIM, KMAX)  
      REAL*8       WXVOR( MAXN, KMAX, 0:NTR)  
      REAL*8       WXDIV( MAXN, KMAX, 0:NTR)  
      REAL*8       WXTMP( MAXN, KMAX, 0:NTR)  
      REAL*8       WXPS ( MAXN, 0:NTR )
      REAL*8       WXSPH( MAXN, KMAX, 0:NTR)  
*
      REAL*8       ALON( IDIM, JMAX )
      REAL*8       ALAT( IDIM, JMAX )
      REAL*8       SIG( KMAX )
      REAL*8       DLON( IDIM, JMAX ) !! unused
      REAL*8       SIGM( KMAX+1 )     !! unused
      REAL*8       DSIG( KMAX )       !! unused
      REAL*8       DSIGM( KMAX+1 )    !! unused
      CHARACTER    HALON *(16)        !! unused
      CHARACTER    HSIG *(16)         !! unused
      CHARACTER    HSIGM *(16)        !! unused
*
*     [work]
      REAL*8       AZ( KMAX )
      REAL*8       AXY( IDIM, JMAX)
      REAL*8       DUM( IMAX, JMAX)
      REAL*8       DUMX( IMAX, JMAX)
      REAL*8       PI
      REAL*8       V1, V2, V
      INTEGER      NMO   ( 2, 0:MMAX, 0:LMAX ) !! order of spect. suffix
      INTEGER      IFM, IFG
      INTEGER      I, J, K, L, M, IJ, KV, IW, JW( 0:NTR ), LEND
      INTEGER      NOUT
      REAL*8       S2D
*
*     [intrinsic]
      INTRINSIC    DSIN, DSQRT, DBLE, SNGL
*
C     2022.11.10: Ziming Chen 
      REAL         r_Vor(IMAX, JMAX, KMAX), r_Div(IMAX, JMAX, KMAX), 
     &             r_Temp(IMAX, JMAX, KMAX), r_Hum(IMAX, JMAX, KMAX), 
     &             r_PS(IMAX, JMAX)
      REAL         lonin(IMAX), latin(JMAX), pre(KMAX)
      integer      status, ncid, retval
      INTEGER      ilon_dimid, nlonin, ilon_varid
      INTEGER      ilat_dimid, nlatin, ilat_varid
      INTEGER      ipre_dimid, nprein, ipre_varid
      INTEGER      nlon, nlat, np
      INTEGER      ivor_varid, idiv_varid, itemp_varid, ips_varid
C     2022.11.10: Ziming Chen 

      CHARACTER    CHPR( MAXH )*15
      CHARACTER    CVPR( MAXV )*15
*
      LOGICAL      OVOR         !! set vorticity forcing ? 
      LOGICAL      ODIV         !! set divergence forcing ? 
      LOGICAL      OTMP         !! set tmperature forcing ? 
      LOGICAL      OPS          !! set sfc. pressure forcing ? 
      LOGICAL      OSPH         !! set humidity forcing ? 
      LOGICAL      OWALL        !! write all the wavenumber
      LOGICAL      OCLASSIC     !! classic dry model
C       2022.11.04: Ziming Chen
      CHARACTER*500 CFM         !! output file name (matrix)
      CHARACTER*500 CFG         !! output file name (grads)
      CHARACTER*500 CFG_mine    !! my forcing pattern
C       2022.11.04: Ziming Chen
      INTEGER      KHPR         !! type of horizontal profile
      INTEGER      KVPR         !! type of vertical profile
      REAL*8       HAMP         !! horizontal amplitude
      REAL*8       VAMP         !! vertical amplitude
      REAL*8       XDIL         !! x-dilation for Elliptic-shape
      REAL*8       YDIL         !! y-dilation for Elliptic-shape
      REAL*8       VDIL         !! dilation for Gamma-shape
      REAL*8       XCNT         !! x-center for Elliptic-shape
      REAL*8       YCNT         !! y-center for Elliptic-shape
      REAL*8       VCNT         !! vertical center for Gamma-shape
      REAL*8       FACT( 5 )    !! factor (unused)

      NAMELIST    /NMFIN/   CFM, CFG, CFG_mine, FACT
      NAMELIST    /NMVAR/   OVOR, ODIV, OTMP, OPS, OSPH
      NAMELIST    /NMHPR/   KHPR, HAMP, XDIL, YDIL, XCNT, YCNT
      NAMELIST    /NMVPR/   KVPR, VAMP, VDIL, VCNT
      NAMELIST    /NMALL/   OWALL
      NAMELIST    /NMCLS/   OCLASSIC

*
      DATA         NOUT / 6 /
      DATA         S2D / 86400.D0 /
      DATA         CHPR
     $           /'Elliptic       ','Zonally uniform'/
      DATA         CVPR
     $           /'Sinusoidal     ','Gamma          ',
     $            'Uniform        '/
*                  123456789012345   123456789012345   
      DATA    OVOR, ODIV, OTMP, OPS, OSPH
     &       / .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /
*
*     open the NAMELIST file
*
      REWIND( 1 )
      OPEN ( 0, FILE='IEEE_ERROR')
      OPEN ( 1, FILE = 'SETPAR', STATUS = 'OLD')
*
      READ( 1, NMFIN) 
      REWIND( 1 )
      READ( 1, NMVAR) 
      REWIND( 1 )
      READ( 1, NMHPR) 
      REWIND( 1 )
      READ( 1, NMVPR) 
      REWIND( 1 )
      READ( 1, NMALL) 
      REWIND( 1 )
      READ( 1, NMCLS) 
*
*
      WRITE( NOUT, * ) '### MAKE FORCING MATRIX ###'
*
      CALL SETCONS
      CALL SPSTUP               !! spherical harmonic functions
      CALL SETNMO2
     O     ( NMO   ,
     D     MMAX  , LMAX  , NMAX  , MINT    )
      

      WRITE( NOUT, * )
      WRITE( NOUT, * ) '................'
      WRITE( NOUT, * ) 'Selected shape:'
      WRITE( NOUT, * ) '  Horizontal:',CHPR( KHPR )
      WRITE( NOUT, * ) '  Vertical  :',CVPR( KVPR )
      WRITE( NOUT, * ) '................'
      WRITE( NOUT, * )
*
*     open files
*
C       
C     2022.10.10: Ziming Chen
C     Read by using NetCDF module
      WRITE(*, *) 'Forcing Pattern: ', CFG_mine
      retval = nf_open(CFG_mine, nf_nowrite, ncid)
      retval = nf_inq_dimid(ncid, 'lon', ilon_dimid)
      retval = nf_inq_dimlen(ncid, ilon_dimid, nlonin)
      retval = nf_inq_dimid(ncid, 'lat', ilat_dimid)
      retval = nf_inq_dimlen(ncid, ilat_dimid, nlatin)
      retval = nf_inq_dimid(ncid, 'lev', ipre_dimid)
      retval = nf_inq_dimlen(ncid, ipre_dimid, nprein)

      nlon   = IMAX
      nlat   = JMAX
      if (nlonin .ne. nlon .or. nlatin .ne. nlat) then
        print*, 'grid error: nlonin,lon,nlatin,lat', nlonin, nlon, 
     &                                               nlatin, nlat
        stop
      endif

      np     = KMAX
      if (nprein .ne. np) then 
        print*, 'pressure levels do not match', nprein, np
        stop
      endif

      retval = nf_inq_varid(ncid, 'lon', ilon_varid)
      retval = nf_get_var(ncid, ilon_varid, lonin)
      retval = nf_inq_varid(ncid, 'lat', ilat_varid)
      retval = nf_get_var(ncid, ilat_varid, latin)
      retval = nf_inq_varid(ncid, 'lev', ipre_varid)
      retval = nf_get_var(ncid, ipre_varid, pre)
C 
      retval = nf_inq_varid(ncid, 'v', ivor_varid)
      retval = nf_get_var_real(ncid, ivor_varid, r_Vor)
C 
	   retval = nf_inq_varid(ncid, 'd', idiv_varid)
      retval = nf_get_var_real(ncid, idiv_varid, r_Div)
C 
      retval = nf_inq_varid(ncid, 't', itemp_varid)
      retval = nf_get_var_real(ncid, itemp_varid, r_Temp)
C       
	   retval = nf_inq_varid(ncid, 'p', ipre_varid)
      retval = nf_get_var_real(ncid, ipre_varid, r_PS)

      retval = nf_close(ncid)
C     Read by using NetCDF module

C C     Read the grd
C       WRITE(*, *) 'Forcing Pattern: ', CFG_mine
C       OPEN ( 20, FILE = CFG_mine, status='OLD', FORM = 'UNFORMATTED', 
C      &  ACCESS = 'direct', RECL = 128*64*20*4)      
C C       READ( 12) r_Temp
C       DO 111 K = 1, KMAX
C           READ( 20, rec = 1 ) ((r_Vor(I, J, K), I = 1, IMAX), 
C      &                            J = 1, JMAX)
C  111  CONTINUE
C       DO 112 K = 1, KMAX
C           READ( 20, rec = 1 ) ((r_Div(I, J, K), I = 1, IMAX), 
C      &                            J = 1, JMAX)
C  112  CONTINUE
C       DO 113 K = 1, KMAX
C           READ( 20, rec = 1 ) ((r_Temp(I, J, K), I = 1, IMAX), 
C      &                            J = 1, JMAX)
C  113  CONTINUE
C       READ( 20, rec = 1 ) ((r_PS(I, J), I = 1, IMAX), J = 1, JMAX)
C       DO 114 K = 1, KMAX
C           READ( 20, rec = 1 ) ((r_Hum(I, J, K), I = 1, IMAX), 
C      &                            J = 1, JMAX)
C  114  CONTINUE

C       CLOSE(20)
C C C     Read the grd
C       WRITE(*, *)'Temp: ',((r_Temp(I, J, 0), I = 1, IMAX), J = 1, JMAX)
C       WRITE(*, *)"Temp: ", MAXVAL(r_Temp), MINVAL(r_Temp)
C       STOP
C C     2022.10.10: Ziming Chen

      IFM = 10
      IFG = 20
      OPEN ( IFM, FILE = CFM, FORM = 'UNFORMATTED', STATUS = 'UNKNOWN' ) 
      OPEN ( IFG, FILE = CFG, FORM = 'UNFORMATTED', STATUS = 'UNKNOWN' )
      WRITE( NOUT, * ) 'Matrix file         :', CFM
      WRITE( NOUT, * ) 'Matrix file (GrADS) :', CFG
      WRITE( NOUT, * ) '................'
      WRITE( NOUT, * ) 
*
      PI = ATAN( 1.0D0 ) * 4.0D0
      CALL SETLON
     O         ( ALON  , DLON , HALON  )
      CALL SETLAT
     O         ( ALAT  , DLON , HALON  )
      CALL SETSIG
     O         ( SIG   , DSIG , HSIG  , 
     O           SIGM  , DSIGM, HSIGM  )
*
      DO 1 J = 1, JMAX
         DO 1 I = 1, IDIM
            ALON( I, J) = ALON( I, J) * 360.D0 / ( 2.D0 * PI )
            ALAT( I, J) = ALAT( I, J) * 360.D0 / ( 2.D0 * PI )
 1    CONTINUE
*
      CALL SETZ( GFRCT, IDIM*JMAX*KMAX )
      CALL SETZ( DUM, IMAX*JMAX )
      CALL SETZ( WFRCT, NMDIM*KMAX )
      CALL SETZ( WFRCF, NMDIM*KMAX )
*
*     vertical profile
*
      IF( KMAX .GT. 1 ) THEN
         DO 10 K = 1, KMAX
            IF( KVPR .EQ. 1 ) AZ( K ) = VAMP * DSIN( PI * SIG( K ) )
            IF( KVPR .EQ. 2 ) AZ( K ) = VAMP * 
     $           DEXP( - VDIL * ( SIG( K ) - VCNT )**2 )
            IF( KVPR .EQ. 3 ) AZ( K ) = VAMP 
 10      CONTINUE
         WRITE( NOUT, * ) 'Set vertical shape'
         WRITE( NOUT, * ) '................'
         WRITE( NOUT, * ) 
      ELSE
         AZ( K ) = VAMP
      ENDIF
*      AZ(1) = 0.3535404
*      AZ(2) = 0.356526
*      AZ(3) = 0.3643257
*      AZ(4) = 0.3846925
*      AZ(5) = 0.4370963
*      AZ(6) = 0.5610808
*      AZ(7) = 0.8154266
*      AZ(8) = 1.232312
*      AZ(9) = 1.712831
*      AZ(10) = 2.098864
*      AZ(11) = 2.304671
*      AZ(12) = 2.332472
*      AZ(13) = 2.240164
*      AZ(14) = 2.07528
*      AZ(15) = 1.903895
*      AZ(16) = 1.784213
*      AZ(17) = 1.709818
*      AZ(18) = 1.658892
*      AZ(19) = 1.607487
*      AZ(20) = 1.522195
*
*     horizontal profile
*
      IJ = 0
      DO 20 J = 1, JMAX
         DO 30 I = 1, IDIM
            IJ = IJ + 1
*
            IF( KHPR .EQ. 1 ) THEN !! Elliptic
               IF( XDIL .LE. 0.D0 .OR. YDIL .LE. 0.D0 ) THEN
                  WRITE( NOUT, *) 
     $                 '### Parameter XDIL/YDIL not correct ###' 
                  STOP
               ENDIF
*
               V1 = ( ALON( I, J ) - XCNT )**2 / XDIL**2
               IF( ALON(I,J).LT.180. .AND. XCNT+XDIL.GE.360. )
     $              V1 = ( 360. + ALON( I, J ) - XCNT )**2 / XDIL**2
               IF( ALON(I,J).GT.180. .AND. XCNT-XDIL.LE.0. )
     $              V1 = ( ALON( I, J ) - 360. - XCNT )**2 / XDIL**2
               V2 = ( ALAT( I, J ) - YCNT )**2 / YDIL**2
               V = V1 + V2
               IF( V .GT. 1. ) THEN
                  AXY( I, J) = 0.D0
               ELSE
                  V = DSQRT( V ) 
                  AXY( I, J) = HAMP * ( 1.D0 - V )
               ENDIF
            ENDIF
               
            IF( KHPR .EQ. 2 ) THEN !! zonal uniform
*
               IF( ALAT( I, J) .GE. YCNT - YDIL .AND.
     $             ALAT( I, J) .LE. YCNT + YDIL ) THEN
                  V = HAMP
               ELSE
                  V = 0.D0
               ENDIF

               AXY( I, J) = V
            ENDIF

            KV = KMAX
            IF( OPS  ) KV = 1 
            DO 40 K = 1, KV
               GFRCT( IJ, K) = AXY( I, J) * AZ( K ) / S2D
 40         CONTINUE
*
 30      CONTINUE
 20   CONTINUE
      WRITE( NOUT, * ) 'Set horizontal shape'
      WRITE( NOUT, * ) '................'
      WRITE( NOUT, * ) 
*
*     write GrADS file
*
      DO 100 K = 1, KMAX
         IF( OVOR ) THEN
            DO 110 J = 1, JMAX
               DO 110 I = 1, IMAX
                  IJ = (J-1)*IDIM + I
C                 2022.11.11: Ziming Chen
C                   DUMX( I, J) = GFRCT( IJ, K)
                  DUMX(I, J) = r_Vor(I, J, K)
C                 2022.11.11: Ziming Chen
 110        CONTINUE
            WRITE( IFG ) ((SNGL(DUMX(I,J)),I=1,IMAX),J=1,JMAX)
         ELSE
            WRITE( IFG ) ((SNGL(DUM(I,J)),I=1,IMAX),J=1,JMAX)
         ENDIF
 100  CONTINUE
*
      DO 120 K = 1, KMAX
         IF( ODIV ) THEN
            DO 130 J = 1, JMAX
               DO 130 I = 1, IMAX
                  IJ = (J-1)*IDIM + I
C                 2022.11.11: Ziming Chen
C                   DUMX( I, J) = GFRCT( IJ, K)
                  DUMX(I, J) = r_Div(I, J, K)
C                 2022.11.11: Ziming Chen
 130        CONTINUE
            WRITE( IFG ) ((SNGL(DUMX(I,J)),I=1,IMAX),J=1,JMAX)
         ELSE
            WRITE( IFG ) ((SNGL(DUM(I,J)),I=1,IMAX),J=1,JMAX)
         ENDIF
 120  CONTINUE
*
      DO 140 K = 1, KMAX
         IF( OTMP ) THEN
            DO 150 J = 1, JMAX
               DO 150 I = 1, IMAX
                  IJ = (J-1)*IDIM + I
C                 2022.11.11: Ziming Chen
C                   DUMX( I, J) = GFRCT( IJ, K)
                  DUMX(I, J) = r_Temp(I, J, K)
C                 2022.11.11: Ziming Chen
 150        CONTINUE
            WRITE( IFG ) ((SNGL(DUMX(I,J)),I=1,IMAX),J=1,JMAX)
         ELSE
            WRITE( IFG ) ((SNGL(DUM(I,J)),I=1,IMAX),J=1,JMAX)
         ENDIF
 140  CONTINUE
*
      IF( OPS ) THEN
         DO 160 J = 1, JMAX
            DO 160 I = 1, IMAX
               IJ = (J-1)*IDIM + I
C                 2022.11.11: Ziming Chen
C                DUMX( I, J) = GFRCT( IJ, 1)
               DUMX(I, J) = r_PS(I, J)
C                 2022.11.11: Ziming Chen
 160     CONTINUE
         WRITE( IFG ) ((SNGL(DUMX(I,J)),I=1,IMAX),J=1,JMAX)
      ELSE
         WRITE( IFG ) ((SNGL(DUM(I,J)),I=1,IMAX),J=1,JMAX)
      ENDIF
*
      IF( .NOT. OCLASSIC ) THEN
         DO 170 K = 1, KMAX
            IF( OSPH ) THEN
               DO 180 J = 1, JMAX
                  DO 180 I = 1, IMAX
                     IJ = (J-1)*IDIM + I
C                 2022.11.11: Ziming Chen
C                      DUMX( I, J) = GFRCT( IJ, K)
                     DUMX(I, J) = r_Hum(I, J, K)
C                 2022.11.11: Ziming Chen
 180           CONTINUE
               WRITE( IFG ) ((SNGL(DUMX(I,J)),I=1,IMAX),J=1,JMAX)
            ELSE
               WRITE( IFG ) ((SNGL(DUM(I,J)),I=1,IMAX),J=1,JMAX)
            ENDIF
 170     CONTINUE
      ENDIF
*
      CLOSE( IFG )
      WRITE( NOUT, * ) 'Written to GrADS file'
      WRITE( NOUT, * ) '................'
      WRITE( NOUT, * ) 
*
*     grid to wave (forcing matrix)
*
      IF( OPS ) THEN
         CALL G2W
     O        ( WFRCT   ,
     I          GFRCT   , '    ', 'POSO', 1    )
      ELSE
         CALL G2W
     O        ( WFRCT   ,
     I          GFRCT   , '    ', 'POSO', KMAX )
      ENDIF
*
*     write down
*
      IW = 0
      DO 200 M = 0, NTR
         LEND = MIN( LMAX, NMAX-M)
         DO 220 K = 1, KMAX
            IW = 0
            DO 210 L = 0, LEND
               IF( M .EQ. 0 .AND. L .EQ. 0 ) GOTO 210
               I = NMO( 1, M, L)
               J = NMO( 2, M, L)
               IW = IW + 1
               WXVOR(IW,K,M)            = WFRCF(I,K)
               WXDIV(IW,K,M)            = WFRCF(I,K)
               WXTMP(IW,K,M)            = WFRCF(I,K)
               WXPS (IW,M  )            = WFRCF(I,1)
               WXSPH(IW,K,M)            = WFRCF(I,K)
               IF( OVOR ) WXVOR(IW,K,M) = WFRCT(I,K)
               IF( ODIV ) WXDIV(IW,K,M) = WFRCT(I,K)
               IF( OTMP ) WXTMP(IW,K,M) = WFRCT(I,K)
               IF(  OPS ) WXPS (IW,M  ) = WFRCT(I,1)
               IF( OSPH ) WXSPH(IW,K,M) = WFRCT(I,K)
               IF( M .EQ. 0 ) GOTO 210
               IW = IW + 1
               WXVOR(IW,K,M)            = WFRCF(J,K)
               WXDIV(IW,K,M)            = WFRCF(J,K)
               WXTMP(IW,K,M)            = WFRCF(J,K)
               WXPS (IW,M  )            = WFRCF(J,1)
               WXSPH(IW,K,M)            = WFRCF(J,K)
               IF( OVOR ) WXVOR(IW,K,M) = WFRCT(J,K)
               IF( ODIV ) WXDIV(IW,K,M) = WFRCT(J,K)
               IF( OTMP ) WXTMP(IW,K,M) = WFRCT(J,K)
               IF(  OPS ) WXPS (IW,M  ) = WFRCT(J,1)
               IF( OSPH ) WXSPH(IW,K,M) = WFRCT(J,K)
  210       CONTINUE
  220    CONTINUE
         IF( .NOT. OWALL ) THEN
            IF( OCLASSIC ) THEN
               WRITE( IFM ) 
     $              ((WXVOR(I,K,M),I=1,IW),K=1,KMAX),
     $              ((WXDIV(I,K,M),I=1,IW),K=1,KMAX),
     $              ((WXTMP(I,K,M),I=1,IW),K=1,KMAX),
     $              ( WXPS (I,M)  ,I=1,IW          )
            ELSE
               WRITE( IFM ) 
     $              ((WXVOR(I,K,M),I=1,IW),K=1,KMAX),
     $              ((WXDIV(I,K,M),I=1,IW),K=1,KMAX),
     $              ((WXTMP(I,K,M),I=1,IW),K=1,KMAX),
     $              ( WXPS (I,M)  ,I=1,IW          ),
     $              ((WXSPH(I,K,M),I=1,IW),K=1,KMAX)
            ENDIF
         ELSE
            JW( M ) = IW
         ENDIF
  200 CONTINUE
      IF( OWALL ) THEN
         IF( OCLASSIC ) THEN
            WRITE( IFM ) 
     $           (((WXVOR(I,K,M),I=1,JW(M)),K=1,KMAX),
     $            ((WXDIV(I,K,M),I=1,JW(M)),K=1,KMAX),
     $            ((WXTMP(I,K,M),I=1,JW(M)),K=1,KMAX),
     $            ( WXPS (I,M)  ,I=1,JW(M)          ),M=0,NTR)
         ELSE
            WRITE( IFM ) 
     $           (((WXVOR(I,K,M),I=1,JW(M)),K=1,KMAX),
     $            ((WXDIV(I,K,M),I=1,JW(M)),K=1,KMAX),
     $            ((WXTMP(I,K,M),I=1,JW(M)),K=1,KMAX),
     $            ( WXPS (I,M)  ,I=1,JW(M)          ),
     $            ((WXSPH(I,K,M),I=1,JW(M)),K=1,KMAX),M=0,NTR)
         ENDIF
      ENDIF
      CLOSE( IFM )
      IF( OWALL ) THEN
         WRITE( NOUT, * ) 'Written to matrix file (all)'
      ELSE
         WRITE( NOUT, * ) 'Written to matrix file (pwm)'
      ENDIF
      WRITE( NOUT, * ) '................'
      WRITE( NOUT, * ) 
*
      WRITE( NOUT, * )
      WRITE( NOUT, * ) '### END OF EXECUTION ###'
      WRITE( NOUT, * )
*
*
  999 STOP
      END
*######################
      SUBROUTINE SETZ( A, IA )
*
*
      INTEGER IA
      REAL*8  A( IA )

      REAL*8  ZERO
      DATA    ZERO / 0.D0 /
      INTEGER I
*
      DO 10 I = 1, IA
         A( I ) = ZERO
 10   CONTINUE

      RETURN
      END
*##########################
      SUBROUTINE SETNMO2    !! order of matrix
     O         ( NMO   ,
     D           MMAX  , LMAX  , NMAX  , MINT    )
*
*   [PARAM] 
      INTEGER    MMAX
      INTEGER    LMAX
      INTEGER    NMAX
      INTEGER    MINT
*
*   [OUTPUT]
      INTEGER    NMO   ( 2, 0:MMAX, 0:LMAX ) !! order of spect. suffix
*
*   [INTERNAL WORK] 
      INTEGER    L, M, MEND, NMH
*
      NMH  = 0
      DO 2200 L = 0, LMAX
         MEND = MIN( MMAX, NMAX-L )
         DO 2100 M = 0, MEND, MINT
            NMH = NMH + 1
            IF ( MMAX .EQ. 0 ) THEN
               NMO ( 1, M, L ) = NMH
               NMO ( 2, M, L ) = NMH
            ELSE
               NMO ( 1, M, L ) = 2* NMH - 1
               NMO ( 2, M, L ) = 2* NMH
            ENDIF
 2100    CONTINUE
 2200 CONTINUE
*
      RETURN
      END
