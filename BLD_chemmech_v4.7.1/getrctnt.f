
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
C @(#)GETRCTNT.F	1.1 /project/mod3/MECH/src/driver/mech/SCCS/s.GETRCTNT.F 02 Jan 1997 15:26:45

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE GETRCTNT ( IMECH, INBUF, IEOL, LPOINT, CHR, WORD,
     &                      NXX, NS, SPCLIS, SPC1RX,
     &                      NWM, NRXWM, NWW, NRXWW, 
     &                      NWO2, NRXWO2, NWN2, NRXWN2,
     &                      NWCH4, NRXWCH4, NWH2, NRXWH2,
     &                      ICOL, IORDER, IRXBITS, IRR )

C=======================================================================
C reads the reactants of each reaction
C
C input: NXX, NS, ICOL, IORDER
C needs: IMECH, INBUF, IEOL, LPOINT, CHR, WORD for GETWORD
C        SPCLIS for LKUPSPEC
C output: ICOL, IORDER, IRXBITS, IRR,
C         NWM
C         NRXWM
C         NWW
C         NRXWW
C         NWO2
C         NRXWO2
C         NWN2
C         NRXWN2
C         NWCH4
C         NRXWCH4
C         NWH2
C         NRXWH2
C
C (possibly updates: NS, INBUF, LPOINT, CHR, WORD)
C=======================================================================
      IMPLICIT NONE
      INCLUDE 'PARMS.e'
      CHARACTER(  1 ) :: CHR
!     CHARACTER(  4 ) :: WORD
      CHARACTER( 16 ) :: WORD
      CHARACTER( 81 ) :: INBUF
      INTEGER IMECH, IEOL, LPOINT
!     CHARACTER(  4 ) :: SPECIES, SPCLIS( MAXSPEC )
      CHARACTER( 16 ) :: SPECIES, SPCLIS( MAXSPEC )
      INTEGER NXX, NS, NSPEC
      INTEGER NWM,     NRXWM( MAX3BODIES )
      INTEGER NWW,     NRXWW( MAX3BODIES )
      INTEGER NWO2,   NRXWO2( MAX3BODIES )
      INTEGER NWN2,   NRXWN2( MAX3BODIES )
      INTEGER NWCH4, NRXWCH4( MAX3BODIES )
      INTEGER NWH2,   NRXWH2( MAX3BODIES )
      INTEGER ICOL
      INTEGER IORDER ( MAXRXNUM )
      INTEGER IRXBITS( MAXRXNUM )
      INTEGER IRR( MAXRXNUM,MAXPRODS+3 )
      INTEGER SPC1RX( MAXSPEC )

      IF ( CHR .NE. '=' ) THEN
         CALL GETWORD ( IMECH, INBUF, LPOINT, IEOL, CHR, WORD )
         SPECIES = WORD
         IORDER( NXX ) = IORDER( NXX ) + 1
         IF ( SPECIES( 1:4 ) .EQ. 'M   ' .OR.
     &        SPECIES( 1:4 ) .EQ. 'm   ' ) THEN
            NWM = NWM + 1
            NRXWM( NWM) = NXX
            IRXBITS( NXX ) = IBSET ( IRXBITS( NXX ), 2 )
         ELSE IF ( SPECIES( 1:4 ) .EQ. 'H2O ' .OR.
     &             SPECIES( 1:4 ) .EQ. 'h2o ') THEN
            NWW = NWW + 1 
            NRXWW( NWW ) = NXX
            IRXBITS( NXX ) = IBSET ( IRXBITS( NXX ), 3 )
         ELSE IF ( SPECIES( 1:4 ) .EQ. 'O2  ' .OR.
     &             SPECIES( 1:4 ) .EQ. 'o2  ') THEN
            NWO2 = NWO2 + 1
            NRXWO2( NWO2 ) = NXX
            IRXBITS( NXX ) = IBSET ( IRXBITS( NXX ), 4 )
         ELSE IF ( SPECIES( 1:4 ) .EQ. 'N2  ' .OR.
     &             SPECIES( 1:4 ) .EQ. 'n2  ') THEN
            NWN2 = NWN2 + 1
            NRXWN2( NWN2 ) = NXX
            IRXBITS( NXX ) = IBSET ( IRXBITS( NXX ), 5 )
         ELSE IF ( SPECIES( 1:4 ) .EQ. 'CH4 ' .OR.
     &             SPECIES( 1:4 ) .EQ. 'ch4 ') THEN
            NWCH4 = NWCH4 + 1
            NRXWCH4( NWCH4 ) = NXX
            IRXBITS( NXX ) = IBSET ( IRXBITS( NXX ), 6 )
         ELSE IF ( SPECIES( 1:4 ) .EQ. 'H2  ' .OR.
     &             SPECIES( 1:4 ) .EQ. 'h2  ') THEN
            NWH2 = NWH2 + 1
            NRXWH2( NWH2 ) = NXX
            IRXBITS( NXX ) = IBSET ( IRXBITS( NXX ), 7 )
         ELSE
            ICOL = ICOL + 1
            CALL LKUPSPEC ( NS, SPECIES, SPCLIS, NXX, SPC1RX, NSPEC )
            IRR( NXX,ICOL ) = NSPEC
         END IF
         IF ( IORDER( NXX ) .GT. 3 ) THEN
            WRITE( *,2001 ) INBUF
            STOP
         END IF
      ELSE
         IORDER( NXX ) = 0
      END IF      ! CHR .NE. '='

      RETURN
2001  FORMAT( / 5X, '*** ERROR: Too many reactants read in -- max=3'
     &        / 5X, 'Last line read was:' / A81 )
      END
