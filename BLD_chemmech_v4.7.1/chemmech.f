
C***************************************************************************
C  Significant portions of Models-3/CMAQ software were developed by        *
C  Government employees and under a United States Government contract.     *
C  Portions of the software were also based on information from non-       *
C  Federal sources, including software developed by research institutions  *
C  through jointly funded cooperative agreements. These research institu-  *
C  tions have given the Government permission to use, prepare derivative   *
C  works, and distribute copies of their work to the public within the     *
C  Models-3/CMAQ software release and to permit others to do so. EPA       *
C  therefore grants similar permissions for use of Models-3/CMAQ software, *
C  but users are requested to provide copies of derivative works to the    *
C  Government without re-strictions as to use by others.  Users are        *
C  responsible for acquiring their own copies of commercial software       *
C  associated with the Models-3/CMAQ release and are also responsible      *
C  to those vendors for complying with any of the vendors' copyright and   *
C  license restrictions. In particular users must obtain a Runtime license *
C  for Orbix from IONA Technologies for each CPU used in Models-3/CMAQ     *
C  applications.                                                           *
C                                                                          *
C  Portions of I/O API, PAVE, and the model builder are Copyrighted        *
C  1993-1997 by MCNC--North Carolina Supercomputing Center and are         *
C  used with their permissions subject to the above restrictions.          *
C***************************************************************************

C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header$

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)CHEMMECH.F 1.1 /project/mod3/MECH/src/driver/mech/SCCS/s.CHEMMECH.F 02 Jan 1997 15:26:41

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      PROGRAM CHEMMECH

      IMPLICIT NONE
      INCLUDE 'PARMS.e'
      INCLUDE 'CHMECH.e'
      CHARACTER(  1 ) :: CHR
      CHARACTER( 16 ) :: WORD
      CHARACTER( 81 ) :: INBUF
      CHARACTER( 32 ) :: MECHNAME = 'MECHDEF'
      CHARACTER(  3 ) :: END
      CHARACTER( 16 ) :: SPCLIS( MAXSPEC )
      INTEGER, EXTERNAL :: INDEX1
      INTEGER IMECH, LPOINT, IEOL
      INTEGER ICOL, ISPC, IRX

      INTEGER NR, IP, NXX, NS, IPR, IPHOTAB, NC
      INTEGER MXPRD                            ! max no. products
      INTEGER NPRDCT( MAXRXNUM )               ! no. of products for rx j
      INTEGER NREACT( MAXRXNUM )               ! no. of reactants for rx j
      INTEGER IPH( MAXPHOTRXNS,3 )
      INTEGER NPHOTAB                          ! no. of photolysis tables
      CHARACTER( 16 ) :: PHOTAB( MAXPHOTRXNS ) ! photolysis table label
      INTEGER IRR( MAXRXNUM,MAXPRODS+3 )
      REAL    SC ( MAXRXNUM,MAXPRODS )
      CHARACTER( 16 ) :: LABEL( MAXRXNUM,2 ) ! LABEL(NXX,1) 1st label found in rx NXX
                                             ! LABEL(NXX,2) 2nd label found in rx NXX
      INTEGER SPC1RX( MAXSPEC )              ! rx index of 1st occurence of species
                                             ! in mechanism table
      CHARACTER( 120 ) :: EQNAME_MECH
      CHARACTER(  32 ) :: DESCRP_MECH
      CHARACTER(  16 ) :: NAMCONSTS( MAXCONSTS ) = (/
     &                    'ATM_AIR',
     &                    'ATM_O2',
     &                    'ATM_N2',
     &                    'ATM_H2',
     &                    'ATM_CH4' /)
      CHARACTER(  16 ) :: CLABEL                  ! mechanism constants label
      REAL( 8 )        :: CONSTVAL                ! retrieved constant
      REAL( 8 )        :: CVAL( MAXCONSTS )       ! mechanism constants value
      INTEGER, PARAMETER :: LUNOUT = 6

      INTEGER, EXTERNAL :: JUNIT
      EXTERNAL NAMEVAL

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Initialize mechanism array variables
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO 101 IRX = 1, MAXRXNUM
         DO ISPC = 1, MAXPRODS+3
            IRR( IRX,ISPC ) = 0
         END DO
         DO ISPC = 1, MAXPRODS
            SC( IRX,ISPC ) = 0.0
         END DO
         DO ISPC = 1, 3
            RTDAT( ISPC,IRX ) = 0.0
         END DO
         KTYPE( IRX ) = 0
         IORDER( IRX )  = 0
         IRXBITS( IRX ) = 0
         LABEL( IRX,1) = '<<<<<<<<<<<<<<<<'
         LABEL( IRX,2) = '>>>>>>>>>>>>>>>>'
         KRX1( IRX ) = 0
         KRX2( IRX ) = 0
         KRX3( IRX ) = 0
         KRX4( IRX ) = 0
         KRX5( IRX ) = 0
         KRX6( IRX ) = 0
         KRX7( IRX ) = 0
101   CONTINUE
      NFALLOFF = 0
      DO 103 IRX = 1, MAXFALLOFF
         IRRFALL( IRX ) = 0   
         DO ISPC = 1, 5
            RFDAT( ISPC,IRX ) = 0.0
         END DO
103   CONTINUE
      DO 105 IRX = 1, MAX3BODIES
         NRXWM( IRX )   = 0
         NRXWW( IRX )   = 0
         NRXWO2( IRX )  = 0
         NRXWN2( IRX )  = 0
         NRXWCH4( IRX ) = 0
         NRXWH2( IRX )  = 0
105   CONTINUE
      KTN1 = 0
      KTN2 = 0
      KTN3 = 0
      KTN4 = 0
      KTN5 = 0
      KTN6 = 0
      KTN7 = 0
!     KCNV = 0
      NWM  = 0
      NWW  = 0
      NWO2 = 0
      NWN2 = 0
      NWCH4 = 0
      NWH2 = 0
      DO ISPC = 1, MAXSPEC
         SPCLIS( ISPC ) = ' '
         SPC1RX( ISPC ) = 0
      END DO
      IP = 0
      NPHOTAB = 0
      DO ISPC = 1, MAXPHOTRXNS
         IPH( ISPC,1 ) = 0
         IPH( ISPC,2 ) = 0
         IPH( ISPC,3 ) = 0
         PHOTAB( ISPC ) = ' '
      END DO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Open mechanism input file and get the first non-comment line
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      NS = 0
      NXX = 0
      MXPRD = 0
      IMECH = JUNIT()
C symbolic link locates "MECHNAME"; setenv requires INQUIRE:
!     OPEN ( UNIT = IMECH, FILE = MECHNAME, STATUS = 'UNKNOWN' )
      CALL NAMEVAL ( MECHNAME, EQNAME_MECH )
      OPEN ( UNIT = IMECH, FILE = EQNAME_MECH, STATUS = 'UNKNOWN' )
      CALL RDLINE ( IMECH, INBUF, LPOINT, IEOL )
      CALL GETCHAR ( IMECH, INBUF, LPOINT, IEOL, CHR )
      WORD( 1:4 ) = '    '
      IF ( CHR .EQ. 'R' .OR. CHR .EQ. 'r' )
     &   WORD( 1:4 ) = INBUF( LPOINT:LPOINT+3 )
      IF ( WORD( 1:4 ) .NE. 'REAC' ) THEN
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C The first string is the mechanism descriptive name 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         DESCRP_MECH = TRIM( ADJUSTL( INBUF ) )
         CALL RDLINE ( IMECH, INBUF, LPOINT, IEOL )
         CALL GETCHAR ( IMECH, INBUF, LPOINT, IEOL, CHR )
         CALL GETWORD ( IMECH, INBUF, LPOINT, IEOL, CHR, WORD )
      ELSE
         DESCRP_MECH = '00000000'
         CALL GETWORD ( IMECH, INBUF, LPOINT, IEOL, CHR, WORD )
!        CALL GETCHAR ( IMECH, INBUF, LPOINT, IEOL, CHR )
      END IF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Make sure this word is REAC and then check for ppm or cm units
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( WORD( 1:4 ) .EQ. 'REAC' ) THEN
         IF ( CHR .EQ. '[' ) THEN
            CALL GETCHAR ( IMECH, INBUF, LPOINT, IEOL, CHR )
            CALL GETWORD ( IMECH, INBUF, LPOINT, IEOL, CHR, WORD )
            IF ( WORD( 1:2 ) .EQ. 'PP' .OR. WORD( 1:2 ) .EQ. 'pp' ) THEN
               KUNITS = 1
            ELSE IF ( WORD( 1:2 ) .EQ. 'CM' .OR. WORD( 1:2 ) .EQ. 'cm' ) THEN
               KUNITS = 2
            ELSE
               WRITE( LUNOUT, 2003 ) INBUF
               STOP ' *** CHEMMECH ERROR ***'
            END IF
            IF ( CHR .NE. ']' ) THEN
               WRITE( LUNOUT, 2005 ) INBUF
               STOP ' *** CHEMMECH ERROR ***'
            END IF
            CALL GETCHAR ( IMECH, INBUF, LPOINT, IEOL, CHR )
         END IF   ! CHR .EQ. '['
      ELSE
         WRITE( LUNOUT, 2007 ) INBUF
         STOP ' *** CHEMMECH ERROR ***'
      END IF      ! word .EQ. 'REAC'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Make sure an equal sign is present processing any reactions 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( CHR .EQ. '=' ) THEN
         CALL GETCHAR ( IMECH, INBUF, LPOINT, IEOL, CHR )
      ELSE
         WRITE( LUNOUT, 2009 ) INBUF
         STOP ' *** CHEMMECH ERROR ***'
      END IF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Start of reaction processing
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
201   CONTINUE
      NXX = NXX + 1
      IF ( NXX .GT. MAXRXNUM ) THEN
         WRITE( LUNOUT, 2011 ) INBUF
         STOP ' *** CHEMMECH ERROR ***'
      END IF 
      IF ( CHR .EQ. '<' ) THEN   ! label for this reaction
         CALL GETLABEL ( IMECH, INBUF, LPOINT, IEOL, CHR, LABEL( NXX,1 ) )
      END IF
      ICOL = 0
      IORDER( NXX ) = 0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Get the reactants
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO 301 ISPC = 1, 3
         CALL GETRCTNT ( IMECH, INBUF, IEOL, LPOINT, CHR, WORD,
     &                   NXX, NS, SPCLIS, SPC1RX,
     &                   NWM, NRXWM, NWW, NRXWW,
     &                   NWO2, NRXWO2, NWN2, NRXWN2,
     &                   NWCH4, NRXWCH4, NWH2, NRXWH2,
     &                   ICOL, IORDER, IRXBITS, IRR )
         IF ( CHR .EQ. '+' ) THEN
            CALL GETCHAR ( IMECH, INBUF, LPOINT, IEOL, CHR )
         ELSE
            GO TO 303
         END IF
301   CONTINUE
303   CONTINUE
      NREACT( NXX ) = ICOL
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Check for equal sign after all reactants read
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( CHR .EQ. '=' ) THEN
         CALL GETCHAR ( IMECH, INBUF, LPOINT, IEOL, CHR )
      ELSE
         WRITE( LUNOUT, 2013 ) INBUF( 1:IEOL )
         STOP ' *** CHEMMECH ERROR ***'
      END IF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Get the products
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ICOL = 3
      IF ( CHR .NE. '#' .AND. CHR .NE. '%' ) THEN
401      CONTINUE
         CALL GETPRDCT ( IMECH, INBUF, LPOINT, IEOL, CHR, WORD,
     &                   NXX, NS, SPCLIS, SPC1RX,
     &                   ICOL, IRXBITS, IRR, SC )
!        IF ( CHR .EQ. '+' .OR. CHR .EQ. '-' ) THEN
!           CALL GETCHAR ( IMECH, INBUF, LPOINT, IEOL, CHR )
         IF ( CHR .NE. '#' .AND. CHR .NE. '%' ) THEN
            GO TO 401
         END IF
      END IF
      ICOL = ICOL - 3
      NPRDCT( NXX ) = ICOL
      MXPRD = MAX ( MXPRD, ICOL )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Check for start of rate constant after all products read
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
!     IF ( CHR .NE. '#' .AND. CHR .NE. '%' ) THEN
!        WRITE( LUNOUT, 2013 ) INBUF( 1:IEOL )
!        STOP ' *** CHEMMECH ERROR ***'
!     END IF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Get rate constants and check for end of reaction symbol
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CALL GETRATE ( IMECH, INBUF, LPOINT, IEOL, CHR,
     &               NXX, LABEL,
     &               KTYPE, IRXBITS, NFALLOFF, IRRFALL,
     &               IPH, IP, NPHOTAB, PHOTAB,
     &               KTN1, KTN2, KTN3, KTN4, KTN5, KTN6, KTN7,
     &               KRX1, KRX2, KRX3, KRX4, KRX5, KRX6, KRX7,
     &               RTDAT, RFDAT )
      IF ( CHR .NE. ';' ) THEN
         WRITE( LUNOUT, 2015 ) INBUF
         STOP ' *** CHEMMECH ERROR ***'
      END IF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Check for end of mechanism; if not, go back to 201 and get the
C next reaction
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CALL GETCHAR ( IMECH, INBUF, LPOINT, IEOL, CHR )
      END = '   '
      IF ( CHR .EQ. 'E' .OR. CHR .EQ. 'e' ) END = INBUF( LPOINT:LPOINT+2 )
      IF ( END .NE. 'END' .AND. END .NE. 'end' ) GO TO 201
      NR = NXX
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Get mechanism constant values for NRXWM, NRXWO2, NRXWN2, NRXWCH4, and NRXWH2 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CVAL = 0.0D0   ! array
      NC = 0
      CALL RDLINE ( IMECH, INBUF, LPOINT, IEOL )
      IF ( IEOL .GE. 0 ) THEN   ! not end of mechanism
         WORD( 1:9 ) = INBUF( IEOL-8:IEOL )
         IF ( WORD( 1:9 ) .EQ. 'CONSTANTS' .OR.
     &        WORD( 1:9 ) .EQ. 'constants' .OR.
     &        WORD( 1:9 ) .EQ. 'Constants' ) THEN
         CALL RDLINE ( IMECH, INBUF, LPOINT, IEOL )
         CALL GETCHAR ( IMECH, INBUF, LPOINT, IEOL, CHR )
421         CONTINUE
            NC = NC + 1
            IF ( CHR .EQ. '<' ) THEN   ! label for this constant
               CALL GETLABEL ( IMECH, INBUF, LPOINT, IEOL, CHR, CLABEL )
            END IF
            CALL GETWORD ( IMECH, INBUF, LPOINT, IEOL, CHR, CLABEL )
            IF ( CHR .EQ. '=' ) THEN
               CALL GETCHAR ( IMECH, INBUF, LPOINT, IEOL, CHR )
               CALL GETREAL ( IMECH, INBUF, LPOINT, IEOL, CHR, CONSTVAL )
            END IF
            IPR = INDEX1 ( CLABEL, MAXCONSTS, NAMCONSTS )
            CVAL( IPR ) = CONSTVAL
            END = '   '
            IF ( CHR .EQ. 'E' .OR. CHR .EQ. 'e' ) END = INBUF( LPOINT:LPOINT+2 )
            IF ( END .NE. 'END' .AND. END .NE. 'end' ) GO TO 421
         END IF
      END IF   ! not end of mechanism

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Resolve all reactions label references
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      WRITE( LUNOUT, * ) ' '
      DO 501 IRX = 1, NR
         IF ( LABEL( IRX,2 ) .NE. '>>>>>>>>>>>>>>>>' ) THEN
C search all rx's for LABEL(,1) match ...
            DO NXX = 1, NR
               IF ( LABEL( NXX,1 ) .EQ. LABEL( IRX,2 ) ) THEN
Could be a linear dependency photolysis reaction ...
                  IF ( KTYPE( IRX ) .EQ. 0 ) THEN
                     DO IPR = 1, IP
                        IF ( IPH( IPR,3 ) .EQ. 0 ) THEN
                           IPH( IPR,2 ) = NXX
                           WRITE( LUNOUT, 1003 ) IRX, IPH( IPR,2 ), RTDAT( 1,IRX ) 
1003                       FORMAT(  3X, 'Reaction', I4,
     &                              1X, 'is proportional to photolysis',
     &                              1X, 'reaction', I4, ','
     &                            /T51, 'scaled by:', 1PG13.5 )
                           GO TO 501
                        END IF
                     END DO
                     WRITE( LUNOUT, 2017 ) IRX
                     STOP ' *** CHEMMECH ERROR ***'
                  END IF
Could be a reverse equilibrium reaction ...
                  IF ( KTYPE( IRX ) .EQ. 5 ) THEN
                     RTDAT( 3,IRX ) = FLOAT ( NXX )
                     WRITE( LUNOUT, 1005 ) IRX, NXX
1005                 FORMAT( 3X, 'Reaction', I4,
     &                       1X, 'is a reverse equilibrium reaction:',
     &                       1X, 'The forward reaction is', I4 )
                     GO TO 501
Could a linear dependency reaction ...
                  ELSE IF ( KTYPE( IRX ) .EQ. 6 ) THEN
                     RTDAT( 2,IRX ) = FLOAT ( NXX )
                     WRITE( LUNOUT, 1007 ) IRX, NXX, RTDAT( 1,IRX )
1007                 FORMAT( 3X, 'Reaction', I4,
     &                       1X, 'is proportional to reaction', I4, ',',
     &                       2X, 'scaled by:', 1PG13.5 )
                     GO TO 501
                  ELSE
                     WRITE( LUNOUT, 2019 ) IRX
                     STOP ' *** CHEMMECH ERROR ***'
                  END IF
               END IF  ! LABEL(NXX,1) .EQ. LABEL( IRX,2)
            END DO
C if we get here, LABEL(,1) match not found
            WRITE( LUNOUT, 2021 ) IRX
            STOP ' *** CHEMMECH ERROR ***'
         END IF  ! LABEL .NE.  ...
501   CONTINUE

C Error-check phot tables and report to log
      WRITE( LUNOUT, * ) ' '
      IPHOTAB = 0
      DO IPR = 1, IP
         IF ( IPH( IPR,3 ) .NE. 0 ) THEN ! table
            IPHOTAB = IPHOTAB + 1
            IRX = IPH( IPR,1 )
            NXX = IPH( IPR,2 )
            WRITE( LUNOUT, 1009 ) IRX, PHOTAB( NXX ), RTDAT( 1,IRX ) 
1009        FORMAT(  3X, 'Reaction', I4,
     &               1X, 'uses photolysis table: ', A16,
     &               1X, 'scaled by:', 1PG13.5 )
         END IF
      END DO
      WRITE( LUNOUT, 1011 ) IPHOTAB, NPHOTAB
1011  FORMAT(/ 5X, 'There are', I3,
     &         1X, 'photolysis table references out of', I3,
     &         1X, 'tables' / )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Create the Fortran INCLUDE files for the reactions data and the
C interim species include file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CALL WREXTS ( EQNAME_MECH,
     &              DESCRP_MECH,
     &              NS, SPCLIS, SPC1RX,
     &              NR,
     &              MXPRD,
     &              IRR,
     &              RTDAT,
     &              SC,
     &              NPRDCT,
     &              NREACT,
     &              KUNITS,
     &              KTYPE,
     &              IRXBITS,
     &              IORDER,
     &              KTN1, KRX1,
     &              KTN2, KRX2,
     &              KTN3, KRX3,
     &              KTN4, KRX4,
     &              KTN5, KRX5,
     &              KTN6, KRX6,
     &              KTN7, KRX7,
     &              NFALLOFF,
     &              IRRFALL,
     &              RFDAT,
     &              NWM,   NRXWM,
     &              NWW,   NRXWW,
     &              NWO2,  NRXWO2,
     &              NWN2,  NRXWN2,
     &              NWCH4, NRXWCH4,
     &              NWH2,  NRXWH2,
     &              LABEL( 1,1 ),
     &              IP,
     &              IPH,
     &              NPHOTAB,
     &              PHOTAB,
     &              NAMCONSTS,
     &              CVAL )

2003  FORMAT( / 5X, '*** ERROR: Units must be either cm, CM, PPM, or ppm'
     &        / 5X, 'Last line read was:' / A81 )
2005  FORMAT( / 5X, '*** ERROR: End bracket, ], missing for units code'
     &        / 5X, 'Last line read was:' / A81 )
2007  FORMAT( / 5X, '*** ERROR: First word of reaction input must be REAC'
     &        / 5X, 'Last line read was:' / A81 )
2009  FORMAT( / 5X, '*** ERROR: Equal sign expected after REACTIONS'
     &        / 5X, 'Last line read was:' / A81 )
2011  FORMAT( / 5X, '*** ERROR: Maximum number of reactions exceeded'
     &        / 5X, 'Last line read was:' / A81 )
2013  FORMAT( / 5X, '*** ERROR: Equal sign expected after reactants'
     &        / 5X, 'Last line read was:' / A81 )
!013  FORMAT( / 5X, '*** ERROR: Rate constant data must begin with a # or %'
!    &        / 5X, 'Last line read was:' / A81 )
2015  FORMAT( / 5X, '*** ERROR: Reactions line must end with a ;'
     &        / 5X, 'Last line read was:' / A81 )
2017  FORMAT( / 5X, '*** ERROR: Linear dependency photolysis reaction label',
     &          1X, 'points to undefined reaction'
     &        / 5X, 'Processing for reaction number:', I6 )
2019  FORMAT( / 5X, '*** ERROR: Reaction label refers to undefined reaction type'
     &        / 5X, 'Processing for reaction number:', I6 )
2021  FORMAT( / 5X, '*** ERROR: Label points to currently undefined reaction'
     &        / 5X, 'Processing for reaction number:', I6 )

      WRITE( LUNOUT, * ) '   Normal Completion of CHEMMECH'
      STOP
      END
