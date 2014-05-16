
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
      SUBROUTINE WRT_RATE_CONSTANT( NR, IP, LABEL, NS, SPCLIS  )


      USE MECHANISM_DATA
      USE KPP_DATA
      
      IMPLICIT NONE

      INTEGER,         INTENT( IN ) :: NR ! number of reactions
      INTEGER,         INTENT( IN ) :: IP ! number of photolysis reaction
      CHARACTER( 16 ), INTENT( IN ) :: LABEL( MAXRXNUM,2 ) ! LABEL(NXX,1) 1st label found in rx NXX
                                                            ! LABEL(NXX,2) 2nd label found in rx NXX
      INTEGER,         INTENT( IN ) :: NS ! number of species
      CHARACTER( 16 ), INTENT( IN ) :: SPCLIS( MAXSPEC )

c..local Variables for steady-state species

       
      CHARACTER(  1 ) :: CHR
      CHARACTER( 16 ) :: WORD
      CHARACTER( 37 ) :: PHRASE
      CHARACTER( 81 ) :: INBUF
      CHARACTER( 16 ) :: RXNS_MODULE_DATA = 'RXNS_MODULE_DATA'
      CHARACTER(  3 ) :: ENDD

      INTEGER, EXTERNAL :: INDEX1
      INTEGER, EXTERNAL :: INDEXES
      INTEGER            :: LPOINT, IEOL
      INTEGER            :: ICOL, ISPC, IRX, IDX
      INTEGER            :: NXX, IPR, IPHOTAB, NC
      INTEGER            :: MXRCT                         ! max no. of reactants
      INTEGER            :: DUMMY_COEF( MAXRXNUM )        ! Yields for the DUMMY variable in each reaction
      INTEGER            :: SS1RX( MAXNLIST )             ! First reaction occurrence for each SS species
      
c..Variables for species to be dropped from mechanism
      INTEGER         :: N_DROP_SPC = 0
      CHARACTER( 16 ) :: DROP_SPC( MAXNLIST )
      LOGICAL         :: LERROR
      LOGICAL         :: KPP_DUMMY   = .FALSE.
      LOGICAL         :: FIRST_TERM  = .TRUE.
      REAL( 8 )       :: WREXT_COEFFS( MAXSPECTERMS)
      INTEGER         :: WREXT_INDEX(  MAXSPECTERMS)

      INTEGER SPC1RX( MAXSPEC )              ! rx index of 1st occurence of species
                                             ! in mechanism table
      CHARACTER( 120 ) :: EQN_MECH_KPP
      CHARACTER( 120 ) :: SPC_MECH_KPP
      CHARACTER( 891 ) :: REACTION_STR(  MAXRXNUM )
      CHARACTER(  16 ) :: COEFF_STR
      CHARACTER(  16 ) :: NAMCONSTS( MAXCONSTS ) = (/
     &                    'ATM_AIR         ',
     &                    'ATM_O2          ',
     &                    'ATM_N2          ',
     &                    'ATM_H2          ',
     &                    'ATM_CH4         ' /)

      CHARACTER(  16 )    :: CLABEL                  ! mechanism constants label
      REAL( 8 )           :: CONSTVAL                ! retrieved constant
      REAL( 8 )            :: CVAL( MAXCONSTS )       ! mechanism constants value
      INTEGER, PARAMETER  :: LUNOUT = 6


      CHARACTER(  12 ) :: EXFLNM_SPCS = 'SPCSDATX'
      CHARACTER(  12 ) :: EXFLNM_RXDT = 'RXNSDATX'
      CHARACTER(  12 ) :: EXFLNM_RXCM = 'RXNSCOMX'

      INTEGER, EXTERNAL :: JUNIT
      INTEGER            :: ICOUNT, IPRODUCT, ISP
      EXTERNAL           :: NAMEVAL
      
      CHARACTER( 120 )   :: MSG, XMSG
      INTEGER            :: STATUS

      INTEGER, SAVE :: IFNEVER = 0     ! Flag for counter initialization
      INTEGER, SAVE :: NDLMAX  = 0     ! Max # of PD loss terms in any reaction
      INTEGER, SAVE :: NDPMAX  = 0     ! Max # of PD prod terms in any reaction

   
      INTEGER :: ICLO( NCS2 )        ! Pointer to # of ops in decomp loop 1
      INTEGER :: JCLO( NCS2 )        ! Pointer to # of ops in decomp loop 2
      INTEGER :: NSPECT( NCS )       ! Number of species in mechanism ncs
      INTEGER :: DBUFF( MXARRAY )
      INTEGER, ALLOCATABLE, SAVE :: ISAPORL( : )  ! Count of PD terms for each species

      INTEGER, ALLOCATABLE, SAVE :: ISPARDER( :,: )  ! Indicator of a PD term in the 
                                                     ! Jacobian matrix
      INTEGER, ALLOCATABLE, SAVE :: IZILCH  ( :,: )  ! # of nonzero calcs in decomp
                                                     ! loop 1
      INTEGER, ALLOCATABLE, SAVE :: JZILCH  ( :,: )  ! # of nonzero calcs in decomp
                                                     ! loop 2
      INTEGER, ALLOCATABLE, SAVE :: LZERO   ( :,: )  ! Symbolic Jacobian matrix

      INTEGER, ALLOCATABLE, SAVE :: IZEROI  ( : )  ! Pointer to decomp loop 1 i index
      INTEGER, ALLOCATABLE, SAVE :: IZEROK  ( : )  ! Pointer to decomp loop 1 k index
      INTEGER, ALLOCATABLE, SAVE :: JZERO   ( : )  ! Pointer to decomp loop 2 i index
      INTEGER IOS                  ! status
   

      INTEGER I,J,K,I1,J1,I2       ! Matrix loop indices
      INTEGER IA, IB               ! I,J index holders for decomp loop 2
      INTEGER INEW, JNEW           ! Index for sorted species number
      INTEGER IOLD, JOLD           ! Index for old species number
      INTEGER IPA, KPA             ! I,K index holders for decomp loop 1
      INTEGER IPB, KPB             ! I,K index holders for decomp loop 1
      INTEGER IPROD, JP            ! Species number of a product
      INTEGER IREACT, IR, JR       ! Species number of a reactant
      INTEGER ISP2                 ! Species loop indices
      INTEGER JRE, JPR, IRE        ! Indices for nonzero Jacobian entries 
      INTEGER JZ3, JZ4             ! Counter for calcs in backsub groupings
      INTEGER NP, IAP              ! Product loop indices
      INTEGER IAL, JAL             ! Reactant loop indices
      INTEGER IAR                  ! Pointer to location of PD term
      INTEGER IARRAY2              ! Final # of matrix entries w/ Sp. Mat
      INTEGER ICB                  ! Counter for # of terms in decomp loop 1
      INTEGER ICBSUM               ! Running count of calcs for j index 
                                   ! in decomp loop 1
      INTEGER ICCOUNT              ! Two term op count for decomp loop 1
      INTEGER ICNT                 ! Total op counter for decomp loop 1
      INTEGER ICNTA                ! op. counter for decomp loop 1 w/ Sp Mat 
      INTEGER ICNTB                ! op. counter for decomp loop 1 w/ Sp Mat
      INTEGER IFSUN                ! Day/night loop index
      INTEGER IJSTEP               ! Number of terms to calc in decomp loops
      INTEGER IMINNEW              ! Index holder for sort routine
      INTEGER IMINOLD              ! Index holder for sort routine
      INTEGER IPORR                ! Species number of a product or reactant
      INTEGER JCB                  ! Counter for # of terms in decomp loop 2
      INTEGER JCCOUNT              ! Two term op count for decomp loop 2
      INTEGER JCNT                 ! Total op counter for decomp loop 2 
      INTEGER JCNTA                ! op. counter for decomp loop 2 w/o Sp Mat
      INTEGER JCNTB                ! op. counter for decomp loop 2 w/ Sp Mat
      INTEGER JZ                   ! Loop index for backsub loops
      INTEGER KA                   ! Loop index for decomposition loops
      INTEGER KCNT                 ! op. counter for bksub loop 1 w/ Sp. Mat.
      INTEGER KCNTA                ! op. counter for bksub loop 1 w/o Sp Mat
      INTEGER KNTARRAY             ! Final # of matrix entries w/o Sp. Mat
      INTEGER KOUNT0               ! Initial # of matrix entries w/ Sp. Mat
      INTEGER KOUNT0A              ! Initial # of matrix entries w/o Sp. Mat
      INTEGER KZ                   ! # of nonzero calcs in backsub loop 1
      INTEGER NCSP                 ! Mechanism number NCS+1=day NCS+2=night
      INTEGER NK, NRT              ! Rate and Reactant number 
      INTEGER NLS                  ! Number of loss PD terms
      INTEGER NOCHANG              ! Count of number of species not reacting
      INTEGER NPR                  ! Number of prod PD terms
      INTEGER NQQ                  ! Loop index for Gear order      
      INTEGER NRPP                 ! Reactant plus product loop index
      INTEGER NRX                  ! Reaction loop index
      INTEGER NU                   ! Active reaction count holder
      INTEGER MCNT                 ! op. counter for bksub loop 2 w/ Sp. Mat.
      INTEGER MCNTA                ! op. counter for bksub loop 2 w/o Sp. Mat.
      INTEGER MINVALU              ! Current number of PD terms in sort
      INTEGER MXIARRAY              ! maximum # of components is sparse Jacobain vector
      INTEGER MZ                   ! # of nonzero calcs in backsub loop 2
      INTEGER SPECIAL_TERMS         ! Total # of terms in special rate
      INTEGER COUNT_TERMS           ! Active count of terms in a special rate


      LOGICAL LITE               ! option to omitted specific write statements
  
      INTERFACE 
       SUBROUTINE GETRCTNT ( IMECH, INBUF, IEOL, LPOINT, CHR, WORD,
     &                      NXX, NS, SPCLIS, SPC1RX,
     &                      ICOL, LABEL, N_DROP_SPC, DROP_SPC )
         USE KPP_DATA
         USE MECHANISM_DATA
         IMPLICIT NONE
         INTEGER,         INTENT(   IN  ) :: IMECH
         CHARACTER( 81 ), INTENT( INOUT ) :: INBUF
         INTEGER,         INTENT( INOUT ) :: LPOINT
         INTEGER,         INTENT( INOUT ) :: IEOL
         CHARACTER(  1 ), INTENT( INOUT ) :: CHR
         CHARACTER( 16 ), INTENT( INOUT ) :: WORD
         INTEGER,         INTENT(   IN  ) :: NXX
         INTEGER,         INTENT( INOUT ) :: NS
         CHARACTER( 16 ), INTENT( INOUT ) :: SPCLIS( MAXSPEC )
         INTEGER,         INTENT( INOUT ) :: SPC1RX( MAXSPEC )
         INTEGER,         INTENT( INOUT ) :: ICOL
         CHARACTER( 16 ), INTENT(   IN  ) :: LABEL( MAXRXNUM, 2 )
         INTEGER,         INTENT(   IN  ) :: N_DROP_SPC
         CHARACTER( 16 ), INTENT(   IN  ) :: DROP_SPC( MAXNLIST )
        END SUBROUTINE GETRCTNT
        SUBROUTINE GETPRDCT ( IMECH, INBUF, LPOINT, IEOL, CHR, WORD,
     &                      NXX, NS, SPCLIS, SPC1RX,
     &                      ICOL, N_DROP_SPC, DROP_SPC )
          USE MECHANISM_DATA
          IMPLICIT NONE
          INTEGER,         INTENT(   IN  ) :: IMECH
          CHARACTER( 81 ), INTENT( INOUT ) :: INBUF
          INTEGER,         INTENT( INOUT ) :: LPOINT
          INTEGER,         INTENT( INOUT ) :: IEOL
          CHARACTER(  1 ), INTENT( INOUT ) :: CHR
          CHARACTER( 16 ), INTENT( INOUT ) :: WORD
          INTEGER,         INTENT(   IN  ) :: NXX
          INTEGER,         INTENT( INOUT ) :: NS
          CHARACTER( 16 ), INTENT( INOUT ) :: SPCLIS( MAXSPEC )
          INTEGER,         INTENT( INOUT ) :: SPC1RX( MAXSPEC )
          INTEGER,         INTENT( INOUT ) :: ICOL
          INTEGER,         INTENT(   IN  ) :: N_DROP_SPC
          CHARACTER( 16 ), INTENT(   IN  ) :: DROP_SPC( MAXNLIST )
         END SUBROUTINE GETPRDCT
         SUBROUTINE GETRATE ( IMECH, INBUF, LPOINT, IEOL, CHR,
     &                         NXX, LABEL, IP )
           USE MECHANISM_DATA
           IMPLICIT NONE
           CHARACTER(  1 ), INTENT( INOUT ) :: CHR
           CHARACTER( 81 ), INTENT( INOUT ) :: INBUF
           INTEGER,         INTENT( IN )    :: IMECH
           INTEGER,         INTENT( INOUT ) :: LPOINT
           INTEGER,         INTENT( INOUT ) :: IEOL
           INTEGER,         INTENT( INOUT ) :: IP
           INTEGER,         INTENT( IN )    :: NXX
           CHARACTER( 16 ), INTENT( INOUT ) :: LABEL( MAXRXNUM,2 )
        END SUBROUTINE GETRATE
        SUBROUTINE WREXTS_FORTRAN90 (WRUNIT, EQNAME_MECH, DESCRP_MECH, NS, 
     &                      SPCLIS, SPC1RX, NR, IP,  NAMCONSTS, CVAL, SS1RX, LITE  ) 
          USE MECHANISM_DATA
          IMPLICIT NONE
          INTEGER,          INTENT( IN )  ::  WRUNIT     ! logical write unit no.
          CHARACTER( 120 ), INTENT ( IN ) :: EQNAME_MECH
          CHARACTER(  32 ), INTENT ( IN ) :: DESCRP_MECH
          INTEGER,          INTENT ( IN ) :: NS                ! no. of species found in mechanism table
          CHARACTER(  16 ), INTENT ( IN ) :: SPCLIS( MAXSPEC ) ! species list from mechanism table
          INTEGER,          INTENT ( IN ) :: NR
          INTEGER,          INTENT ( IN ) :: SPC1RX( MAXSPEC ) ! rx index of 1st occurence of species in mechanism table
          INTEGER,          INTENT ( IN ) :: IP
          CHARACTER( 16 ),  INTENT ( IN ) :: NAMCONSTS( MAXCONSTS )
          REAL( 8 ),        INTENT ( IN ) :: CVAL( MAXCONSTS )
          INTEGER,          INTENT ( IN ) :: SS1RX( MAXNLIST )
          LOGICAL,          INTENT ( IN ) :: LITE               ! option to omit specific write statements
        END SUBROUTINE WREXTS_FORTRAN90
        SUBROUTINE GET_SS_DATA ( LUNOUT, NR ) 
          USE MECHANISM_DATA
          IMPLICIT NONE
          INTEGER, INTENT ( IN )         :: LUNOUT   ! Output unit number
          INTEGER, INTENT ( IN )         :: NR       ! No. of reactions
        END SUBROUTINE GET_SS_DATA
        SUBROUTINE CHECK_SS_SPC ( LUNOUT, NS, SPCLIS, NR, LABEL, SS1RX )
         USE MECHANISM_DATA
         IMPLICIT NONE
         INTEGER, INTENT ( IN )         :: LUNOUT               ! Output unit number
         INTEGER, INTENT ( IN )         ::  NS                  ! No. of species in mechanism
         CHARACTER( 16 ), INTENT ( IN ) ::  SPCLIS( MAXSPEC )   ! List of mechanism species
         INTEGER, INTENT ( IN )         ::  NR                  ! No. of reactions
         CHARACTER( 16 ), INTENT ( IN ) ::  LABEL( MAXRXNUM,2 ) ! Reaction labels
         INTEGER, INTENT ( INOUT )      ::  SS1RX( MAXNLIST )
       END SUBROUTINE CHECK_SS_SPC
       SUBROUTINE WRSS_EXT_FORTRAN90( WRUNIT, NR ) 
         USE MECHANISM_DATA
         IMPLICIT NONE
         INTEGER, INTENT( IN )    ::  WRUNIT     ! logical write unit no.
         INTEGER, INTENT ( IN )   :: NR   ! No. of reactions
       END SUBROUTINE WRSS_EXT_FORTRAN90
      END INTERFACE 
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Initialize module and local mechanism array variables
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      IF( .NOT. ALLOCATED( INDEX_FIXED_SPECIES ) )THEN
          ALLOCATE( INDEX_FIXED_SPECIES( MAXRXNUM, MAXRCTNTS ) )
      END IF
      
      N_SPEC = NS
      N_RXNS = NR   ! loads RBDATA from RXCM.EXT

      MXRCT = MAXRCTNTS
      
      MXRR = 3 * MXRCT
      MXRP = 3 * MXPRD

      MXCOUNT1 = N_SPEC * MAXGL3 * 3
      MXCOUNT2 = N_SPEC * MAXGL3 * 3
      
         ALLOCATE( ISAPORL ( NS ),
     &             ISPARDER( NS, NS ),
     &             LZERO   ( NS, NS ),
     &             IZILCH  ( NS, NCS2 ),
     &             JZILCH  ( NS, NCS2 ), STAT = STATUS )
         IF ( STATUS .NE. 0 ) THEN
            XMSG = '*** Memory allocation failed'
            WRITE( 6, * )'WRT_REACTIONS_MODULE : ' // XMSG
            STOP
         END IF


      ALLOCATE( NKUSERAT( NRXNS,NCS2 ),
     &          NDERIVL ( NRXNS,NCS2 ),
     &          NDERIVP ( NRXNS,NCS2 ), STAT = STATUS )
      IF ( STATUS .NE. 0 ) THEN
         XMSG = 'ERROR allocating NKUSERAT, NDERIVL or NDERIVP'
         WRITE( 6, * )'WRT_REACTIONS_MODULE : ' // XMSG
         STOP
      END IF

      ALLOCATE( IRM2( NRXNS,MXRCT+MXPRD,NCS2 ),
     &          ICOEFF( NRXNS,MXRP,NCS2 ), STAT = STATUS )
      IF ( STATUS .NE. 0 ) THEN
         XMSG = 'ERROR allocating IRM2 or ICOEFF'
         WRITE( 6, * )'WRT_REACTIONS_MODULE : ' // XMSG
         STOP
      END IF

      ALLOCATE( JARRAYPT( NS,NS,NCS2 ), STAT = STATUS )
      IF ( STATUS .NE. 0 ) THEN
         XMSG = 'ERROR allocating JARRAYPT'
         WRITE( 6, * )'WRT_REACTIONS_MODULE : ' // XMSG
         STOP
      END IF

      ALLOCATE( JARRL( NRXNS,MXRR,NCS2 ),
     &          JARRP( NRXNS,MXRP,NCS2 ),
     &          JLIAL( NRXNS,MXRR,NCS2 ),
     &          JPIAL( NRXNS,MXRP,NCS2 ), STAT = STATUS )
      IF ( STATUS .NE. 0 ) THEN
         XMSG = 'ERROR allocating JARRL, JARRP, JLIAL, or JPIAL'
         WRITE( 6, * )'WRT_REACTIONS_MODULE : ' // XMSG
         STOP
      END IF

      ALLOCATE( INEW2OLD( N_SPEC,NCS ),
     &          IOLD2NEW( N_SPEC,NCS ), STAT = STATUS )
      IF ( STATUS .NE. 0 ) THEN
         XMSG = 'ERROR allocating INEW2OLD or IOLD2NEW'
         WRITE( 6, * )'WRT_REACTIONS_MODULE : ' // XMSG
         STOP
      END IF

      ALLOCATE( JZEROA( MXARRAY ),
     &          JZEROB( MXARRAY ), STAT = STATUS )
      IF ( STATUS .NE. 0 ) THEN
         XMSG = 'ERROR allocating JZEROA or JZEROB'
         WRITE( 6, * )'WRT_REACTIONS_MODULE : ' // XMSG
         STOP
      END IF

      ALLOCATE( JZLO( NCS2 ),
     &          IDEC1LO( N_SPEC,NCS2 ),
     &          IDEC1HI( N_SPEC,NCS2 ), STAT = STATUS )
      IF ( STATUS .NE. 0 ) THEN
         XMSG = 'ERROR allocating JZLO, IDEC1LO or IDEC1HI'
         WRITE( 6, * )'WRT_REACTIONS_MODULE : ' // XMSG
         STOP
      END IF

      ALLOCATE( IZEROI( MXCOUNT1 ),
     &          IZEROK( MXCOUNT2 ),
     &          JZERO ( MXCOUNT1 ), STAT = STATUS )
      IF ( STATUS .NE. 0 ) THEN
         XMSG = 'ERROR allocating IZERO, IZEROK, JZERO'
         WRITE( 6, * )'WRT_REACTIONS_MODULE : ' // XMSG
         STOP
      END IF


      ALLOCATE( IJDECA( MXCOUNT2 ),
     &          IJDECB( MXCOUNT2 ),
     &          IKDECA( MXCOUNT2 ),
     &          IKDECB( MXCOUNT2 ),
     &          KJDECA( MXCOUNT2 ),
     &          KJDECB( MXCOUNT2 ), STAT = STATUS )
      IF ( STATUS .NE. 0 ) THEN
         XMSG = 'ERROR allocating IJDECA, IJDECB, IKDECA, IKDECB, KJDECA, or KJDECB'
         WRITE( 6, * )'WRT_REACTIONS_MODULE : ' // XMSG
         STOP
      END IF

      ALLOCATE( JHIZ1( N_SPEC,NCS2 ),
     &          JHIZ2( N_SPEC,NCS2 ),
     &          KZLO1( N_SPEC,NCS2 ),
     &          KZLO2( N_SPEC,NCS2 ),
     &          KZHI0( N_SPEC,NCS2 ),
     &          KZHI1( N_SPEC,NCS2 ), STAT = STATUS )
      IF ( STATUS .NE. 0 ) THEN
         XMSG = 'ERROR allocating JHIZ1, JHIZ2, KZLO1, KZLO2, KZHI0, or KZHI1'
         WRITE( 6, * )'WRT_REACTIONS_MODULE : ' // XMSG
         STOP
      END IF

      ALLOCATE( KZERO( MXARRAY,NCS2 ), STAT = STATUS )
      IF ( STATUS .NE. 0 ) THEN
         XMSG = 'ERROR allocating KZERO'
         WRITE( 6, * )'WRT_REACTIONS_MODULE : ' // XMSG
         STOP
      END IF

      ALLOCATE( KZILCH( N_SPEC,NCS2 ),
     &          MZHI0 ( N_SPEC,NCS2 ),
     &          MZHI1 ( N_SPEC,NCS2 ),
     &          MZILCH( N_SPEC,NCS2 ),
     &          MZLO1 ( N_SPEC,NCS2 ),
     &          MZLO2 ( N_SPEC,NCS2 ), STAT = STATUS )
      IF ( STATUS .NE. 0 ) THEN
         XMSG = 'ERROR allocating KZILCH, MZHI0, MZHI1, MZILCH, MZLO1, or MZLO2'
         WRITE( 6, * )'WRT_REACTIONS_MODULE : ' // XMSG
         STOP
      END IF

!      ALLOCATE( VDIAG( BLKSIZE,NS ), STAT = STATUS )
!      IF ( STATUS .NE. 0 ) THEN
!         XMSG = 'ERROR allocating VDIAG'
!         WRITE( 6, * )'WRT_REACTIONS_MODULE : ' // XMSG
!         STOP
!      END IF

!      ALLOCATE( CC2( BLKSIZE,0:MXARRAY ), STAT = STATUS )
!      IF ( STATUS .NE. 0 ) THEN
!         XMSG = 'ERROR allocating CC2'
!         WRITE( 6, * )'WRT_REACTIONS_MODULE : ' // XMSG
!         STOP
!      END IF

         DO NCSP = 1, NCS2
            DO ISPC = 1, NS
               IZILCH( ISPC,NCSP ) = 0
               JZILCH( ISPC,NCSP ) = 0
               JHIZ1 ( ISPC,NCSP ) = 0
               JHIZ2 ( ISPC,NCSP ) = 0
               KZILCH( ISPC,NCSP ) = 0
               MZILCH( ISPC,NCSP ) = 0
            END DO
         END DO

         DO NCSP = 1, NCS2
            NUSERAT( NCSP ) = 0
            DO NK = 1, NRXNS
               NDERIVL( NK,NCSP ) = 0
               NDERIVP( NK,NCSP ) = 0
            END DO
         END DO
         
         DO NCSP = 1, NCS
            ISCHANG( NCSP ) = 0
         END DO


         DO NCSP = 1, NCS2
            DO ISP = 1, N_SPEC
               IZILCH( ISP,NCSP ) = 0
               JZILCH( ISP,NCSP ) = 0
               JHIZ1 ( ISP,NCSP ) = 0
               JHIZ2 ( ISP,NCSP ) = 0
               KZILCH( ISP,NCSP ) = 0
               MZILCH( ISP,NCSP ) = 0
            END DO
         END DO

         DO NCSP = 1, NCS2
            NUSERAT( NCSP ) = 0
            DO NK = 1, NRXNS
               NDERIVL( NK,NCSP ) = 0
               NDERIVP( NK,NCSP ) = 0
            END DO
         END DO
         
         DO NCSP = 1, NCS
            ISCHANG( NCSP ) = 0
         END DO
         
         JARRAYPT = 0

         IJDECA = 0
         IKDECA = 0
         KJDECA = 0

         IJDECB = 0
         IKDECB = 0
         KJDECB = 0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Determine active reactions for day and then night (i.e., photo
c  reactions determined by BTEST=.TRUE. are not included for nighttime)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               
       DO NRX = 1, NRXNS
          IF ( NREACT( NRX ) .GT. 0 ) THEN
             NUSERAT( NCS ) = NUSERAT( NCS ) + 1
             NU = NUSERAT( NCS )
             NKUSERAT( NU, NCS ) = NRX
             IF ( .NOT. ( BTEST ( IRXBITS( NRX ),1 ) ) ) THEN
                NUSERAT( NCS + 1 ) = NUSERAT( NCS + 1 ) + 1
                NU = NUSERAT( NCS + 1 )
                NKUSERAT( NU, NCS + 1 ) = NRX
             END IF
          END IF
       END DO

 
       NCSP = NCS

       ALLOCATE ( NET_EFFECT( N_SPEC, NRXNS ), 
     &            NET_RCOEFF( N_SPEC, NRXNS ) )


       NET_EFFECT    = 0     ! default setting, species not net reactant or product
       NET_RCOEFF = 0.0D0 ! initialize   
 
       DO NK =  1, NRXNS

C...Set NET_EFFECT for reaction product
            DO NP = 1, NPRDCT( NK )
               ISP2 = IRR( NK,NP+3 )
               
               NET_EFFECT( ISP2, NK )    = 3
               
               NET_RCOEFF( ISP2, NK ) = NET_RCOEFF( ISP2, NK )
     &                                  + REAL(SC( NK,NP ), 8)
            END DO
         
C..Check whether reaction has a species as both reactant and product

         
         DO NRT = 1, NREACT( NK )

            ISP = IRR( NK,NRT )
            NET_RCOEFF( ISP, NK ) = NET_RCOEFF( ISP, NK ) - 1.0D0
            DO NP = 1, NPRDCT( NK )
               ISP2 = IRR( NK,NP+3 )

               IF( ISP .EQ. ISP2 )THEN
                   IF( NET_RCOEFF( ISP, NK ) .EQ. 0.0D0 )THEN ! reaction has no net effect

                       NET_EFFECT( ISP, NK ) = 0

                   ELSE IF( NET_RCOEFF( ISP, NK ) .LT. 0.0D0 )THEN ! net loss

                       IF( NET_RCOEFF( ISP, NK ) .EQ. -1.0D0 )THEN ! only loss process
                           NET_EFFECT( ISP, NK ) = -1
                       ELSE
                           NET_EFFECT( ISP, NK ) = -2
                       END IF

                   ELSE IF( NET_RCOEFF( ISP, NK ) .GT. 0.0D0 )THEN ! loss is not 100% 

                       NET_EFFECT( ISP, NK ) = 2

                   END IF
               END IF
            END DO
            IF( NET_RCOEFF( ISP, NK ) .LT. 0.0D0 )THEN
                IF( NET_RCOEFF( ISP, NK ) .EQ. -1.0D0 )THEN
                     NET_EFFECT( ISP, NK ) = -1
                ELSE
                     NET_EFFECT( ISP, NK ) = -2
                END IF
            END IF 
            
         END DO                      

         WRITE(6,'(5A,I2,A,ES12.4)')'For reactant ', TRIM(MECHANISM_SPC( ISP )),' : reaction ',
     &   RXLABEL( NK ),' NET_EFFECT = ',NET_EFFECT( ISP, NK ),' NET_RCOEFF = ', NET_RCOEFF( ISP, NK )
                       
                   
        END DO               ! END LOOP OVER REACTIONS

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Initialize Prod/loss and PD tabulator arrays
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      NSPECT( NCS ) = N_SPEC
      DO ISP = 1, NSPECT( NCS )
         ISAPORL( ISP ) = 0
      END DO

      DO ISP = 1, NSPECT( NCS )
         DO ISP2 = 1, NSPECT( NCS )
            ISPARDER( ISP,ISP2 ) = 0
         END DO
      END DO
   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set the number of Partial derivative terms in the Jacobian and
c  count the number of terms for each species
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO NRX = 1, NRXNS
        DO NRT = 1, 3
            IREACT = IRR( NRX,NRT )
            IF ( IREACT .EQ. 0 )CYCLE
            DO NRPP = 1, 3 + MXPRD
               IPORR = IRR( NRX,NRPP )
               IF ( IPORR .EQ. 0 )CYCLE
               IF( ABS(NET_RCOEFF( IPORR, NRX )) .GT. 1.0D-6 )ISPARDER( IPORR,IREACT ) = 1               
            END DO
         END DO
      END DO

      DO IREACT = 1, NSPECT( NCS ) 
         DO IPORR = 1, NSPECT( NCS )
            IF ( ISPARDER( IPORR,IREACT ) .EQ. 1 ) 
     &           ISAPORL( IPORR ) = ISAPORL( IPORR ) + 1
         END DO
      END DO
 
      DO IPORR = 1, NSPECT( NCS )
         PRINT*,TRIM(SPCLIS( IPORR )),' has ', ISAPORL( IPORR ),' partial derivatives'
      END DO
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Sort the species, putting all with zero partial derivative 
c  terms at the bottom and those with fewest PD terms at top.
c  Set arrays for species with zero PD terms
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      NOCHANG = NSPECT( NCS )
      ISCHAN   = N_SPEC
      DO JOLD = 1, NSPECT( NCS )
         IF ( ISAPORL( JOLD ) .GT. 0 ) THEN
            ISCHANG( NCS ) = ISCHANG( NCS ) + 1
            JNEW = ISCHANG( NCS )
            INEW2OLD( JNEW,NCS ) = JOLD
            IOLD2NEW( JOLD,NCS ) = JNEW
         ELSE
            INEW2OLD( NOCHANG,NCS ) = JOLD
            IOLD2NEW( JOLD,NCS ) = NOCHANG
            NOCHANG = NOCHANG - 1
         END IF
      END DO
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Now sort by number of PD terms, fewest at position 1, most at
c  the end position. 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO JNEW = 1, ISCHANG( NCS )
c  Uncomment the following three lines to turn off species ordering;
c  not recommended since computational efficiency reduced
!        INEW2OLD( JNEW,NCS ) = JNEW
!        IOLD2NEW( JNEW,NCS ) = JNEW
!        IF ( JNEW .NE. 0 ) CYCLE
         JOLD = INEW2OLD( JNEW,NCS )
         MINVALU = ISAPORL( JOLD )
         IMINOLD = JOLD
         IMINNEW = JNEW

         DO INEW = JNEW + 1, ISCHANG( NCS )
            IOLD = INEW2OLD( INEW,NCS )
            IF ( ISAPORL( IOLD ) .LT. MINVALU ) THEN
               MINVALU = ISAPORL( IOLD )
               IMINOLD = IOLD
               IMINNEW = INEW
            END IF
         END DO

         INEW2OLD( IMINNEW,NCS ) = JOLD
         INEW2OLD( JNEW,NCS )    = IMINOLD
         IOLD2NEW( JOLD,NCS )    = IMINNEW
         IOLD2NEW( IMINOLD,NCS ) = JNEW
      END DO

       DO ISP = 1, NSPECT( NCS )
          PRINT*,'ISP = ',ISP,'SORTED SPECIES( ', INEW2OLD(ISP,NCS),' ) = ',
     &   TRIM( MECHANISM_SPC(INEW2OLD(ISP,NCS)) )
       END DO
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Fill the irm2 array using the new species order developed above.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO NRX = 1, NRXNS
         DO NRT = 1, NREACT( NRX )
            IREACT = IRR( NRX,NRT )
            IRM2( NRX,NRT,NCS ) = IOLD2NEW( IREACT,NCS )
         END DO

         DO NP = 1, NPRDCT( NRX )
            IPROD = IRR( NRX, NP + 3 )
            IRM2( NRX,NP+3,NCS ) = IOLD2NEW( IPROD,NCS )
         END DO
         
!         IF ( NREACT( NRX ) .GT. 0 ) THEN
!            NUSERAT( NCS ) = NUSERAT( NCS ) + 1
!            NU = NUSERAT( NCS )
!           NKUSERAT( NU, NCS ) = NRX
!            IF ( .NOT. ( BTEST ( IRXBITS( NRX ),1 ) ) ) THEN
!               NUSERAT( NCS + 1 ) = NUSERAT( NCS + 1 ) + 1
!               NU = NUSERAT( NCS + 1 )
!               NKUSERAT( NU, NCS + 1 ) = NRX
!            END IF
!         END IF
      END DO



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Do symbolic LU decomposition to determine sparse storage array
c  structures. Done twice, first for day and then for night. An entry
c  of 1 in lzero means a non-negative entry in the Jacobian. First
c  put ones on the diagonal and zeroes everywhere else.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO 700 IFSUN = 1, 2
         NCSP = IFSUN
         DO I = 1, ISCHANG( NCS )
            DO J = 1, ISCHANG( NCS )
               LZERO( J,I ) = 0
            END DO
            LZERO( I,I ) = 1
         END DO
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Fill in the rest of the entries in the Jacobian
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         DO NRX = 1, NUSERAT( NCSP )
            NK = NKUSERAT( NRX,NCSP )
            DO NRT = 1, NREACT( NK )
               IRE  = IRM2( NK,NRT,NCS )
               DO JAL = 1, NREACT( NK )
                  JRE = IRM2( NK,JAL,NCS )
                  JOLD = INEW2OLD( JRE, NCS )
                  IF( ABS(NET_RCOEFF( JOLD, NK )) .GT. 1.0D-6 )THEN
                      LZERO( JRE,IRE ) = 1
                  END IF
               END DO
               DO IAP = 1, NPRDCT( NK )
                  JPR = IRM2( NK,3+IAP,NCS )
                  JOLD = INEW2OLD( JPR, NCS )
                  IF( ABS(NET_RCOEFF( JOLD, NK )) .GT. 1.0D-6 )THEN
                      LZERO( JPR,IRE ) = 1 
                  END IF 
               END DO
           END DO
         END DO
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set up arrays for decomposition / back-substitution of sparse     
c  matrices by removing all calculations involving a zero.          
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( IFNEVER.EQ.0 ) THEN
            IFNEVER = 1
            ICNT    = 0 
            JCNT    = 0 
            ICCOUNT = 0
            JCCOUNT = 0
         END IF
         KOUNT0A = 0
         KOUNT0  = 0
         ICNTA   = 0
         ICNTB   = 0
         JCNTA   = 0
         JCNTB   = 0
         KCNTA   = 0
         MCNTA   = 0
         KCNT    = 0
         MCNT    = 0
         IARRAY2 = 0
         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Count number of entries w/ and w/o sparse matrix storage
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
         DO J = 1, ISCHANG( NCS )
            DO K = 1, ISCHANG( NCS )
               KOUNT0A = KOUNT0A + 1
               IF ( LZERO( J,K ) .EQ. 1 ) KOUNT0 = KOUNT0 + 1
            END DO
         END DO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Do the symbolic decomposition (ludcmp) converting [A] to [L][U] 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ICLO( NCSP ) = ICNT + 1
         JCLO( NCSP ) = JCNT + 1
         DO J = 1, ISCHANG( NCS )
            J1 = J - 1
            
c...  First loop of decomposition
            DO I = 2, ISCHANG( NCS ) 
               I1 = J1 
               IF ( I .LE. J1 ) I1 = I - 1
               DO K = 1, I1
                  ICNTA = ICNTA + 1
                  IF ( LZERO( I,K ) .EQ. 1 .AND. LZERO( K,J ) .EQ. 1 )
     &               THEN
                     IZILCH( J,NCSP ) = IZILCH( J,NCSP ) + 1
                     ICNT             = ICNT + 1
                     ICNTB            = ICNTB + 1
                     IZEROK( ICNT )   = K   
                     IZEROI( ICNT )   = I
                     LZERO( I,J )     = 1 
                  END IF
               END DO
            END DO
c... Second loop of decomposition 
            DO I = J + 1, ISCHANG( NCS ) 
               JCNTA = JCNTA + 1
               IF ( LZERO( I,J ) .EQ. 1 ) THEN
                  JZILCH( J,NCSP ) = JZILCH( J,NCSP ) + 1
                  JCNT             = JCNT  + 1
                  JCNTB            = JCNTB + 1
                  JZERO( JCNT )    = I  
               END IF
            END DO 
         END DO
  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Do symbolic back-substition for solving [L][U]{x}={b}. Store data
c  in sparse matrix pointer jarraypt.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c... First loop of back-substitution
         DO I = 2, ISCHANG( NCS )
            I1 = I - 1
            DO J = 1, I1    
               KCNTA = KCNTA + 1
               IF ( LZERO( I,J ) .EQ. 1 ) THEN 
                  KZILCH( I,NCSP ) = KZILCH( I,NCSP ) + 1
                  KCNT = KCNT + 1
                  IARRAY2 = IARRAY2 + 1
                  KZERO( IARRAY2,NCSP ) = J
                  JARRAYPT( I,J,NCSP ) = IARRAY2 
               END IF
            END DO
         END DO 

c... Second loop of back-substitution 
         DO I = ISCHANG( NCS ) - 1, 1, -1
            I2 = I + 1
            DO J = I + 1, ISCHANG( NCS )
               MCNTA = MCNTA + 1
               IF ( LZERO( I,J ) .EQ. 1 ) THEN 
                  MZILCH( I,NCSP )      = MZILCH( I,NCSP ) + 1
                  MCNT                  = MCNT + 1
                  IARRAY2               = IARRAY2 + 1
                  KZERO( IARRAY2,NCSP ) = J
                  JARRAYPT( I,J,NCSP )  = IARRAY2 
               END IF
            END DO
         END DO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Fill jarraypt with remaining diagonal array points and save counts
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         DO I = 1, ISCHANG( NCS ) 
            IARRAY2 = IARRAY2 + 1
            JARRAYPT( I,I,NCSP ) = IARRAY2 
         END DO
         IARRAY( NCSP ) = IARRAY2 
         PRINT*,'IARRAY( ', NCSP, ' ) = ', IARRAY( NCSP )
         KNTARRAY = KCNTA + MCNTA + ISCHANG( NCS )

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Do decomposition again to change arrays to use jarraypt
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         JCB = JCLO( NCSP ) 
         JZLO( NCSP ) = JCCOUNT
         ICBSUM = ICLO( NCSP ) - 1 
         IJSTEP = 2   

         DO J = 1, ISCHANG( NCS )

c...First loop of decomposition
            IDEC1LO( J,NCSP ) = ICCOUNT + 1
            ICB = ICBSUM  + 1
            ICBSUM = ICBSUM + IZILCH( J, NCSP )

            DO KA = 1, IZILCH( J,NCSP ), IJSTEP
               ICCOUNT = ICCOUNT + 1
               IPA = IZEROI( ICB ) 
               KPA = IZEROK( ICB ) 
               IJDECA( ICCOUNT ) = JARRAYPT( IPA,  J,NCSP ) 
               IKDECA( ICCOUNT ) = JARRAYPT( IPA,KPA,NCSP )
               KJDECA( ICCOUNT ) = JARRAYPT( KPA,  J,NCSP )
               IF ( ICB + 1 .LE. ICBSUM ) THEN
                  IPB = IZEROI( ICB + 1 ) 
                  KPB = IZEROK( ICB + 1 ) 
                  IJDECB( ICCOUNT ) = JARRAYPT( IPB,  J,NCSP ) 
                  IKDECB( ICCOUNT ) = JARRAYPT( IPB,KPB,NCSP )
                  KJDECB( ICCOUNT ) = JARRAYPT( KPB,  J,NCSP )
               END IF
               ICB = ICB + IJSTEP   
            END DO

            IDEC1HI( J,NCSP ) = ICCOUNT  
            
c...Second loop of decomposition
            JZ = JZILCH( J, NCSP )

            DO I = 1, JZ - 1, 2
               JCCOUNT           = JCCOUNT + 1
               JHIZ1( J,NCSP )   = JHIZ1( J,NCSP ) + 1
               IA                = JZERO( JCB )
               IB                = JZERO( JCB + 1 )
               JZEROA( JCCOUNT ) = JARRAYPT( IA,J,NCSP )
               JZEROB( JCCOUNT ) = JARRAYPT( IB,J,NCSP )
               JCB = JCB + 2
            END DO

            IF ( MOD( JZ,2 ) .EQ. 1 ) THEN 
               JCCOUNT           = JCCOUNT + 1
               JHIZ2( J,NCSP )   = JHIZ2( J,NCSP ) + 1
               IA                = JZERO( JCB )
               JZEROA( JCCOUNT ) = JARRAYPT( IA,J,NCSP )
               JCB               = JCB + 1 
            END IF
         END DO

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Group terms to increase efficiency in back-substition
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c... First back-substitution loop
         DO I = 1, ISCHANG( NCS ) 
            KZ              = KZILCH( I,NCSP )
            KZHI0( I,NCSP ) = KZ - 4 
            JZ3             = 0

            DO JZ = 1, KZHI0( I,NCSP ), 5     
               JZ3 = JZ + 4
            END DO  

            KZLO1( I,NCSP ) = JZ3 + 1
            KZHI1( I,NCSP ) = KZ  - 1 
            JZ4 = JZ3 

            DO JZ = JZ3 + 1, KZ - 1, 2    
               JZ4 = JZ + 1
            END DO

            KZLO2( I,NCSP ) = JZ4 + 1
         END DO
 
c... Second loop of back-substitution
         DO I = ISCHANG( NCS ), 1, -1
            MZ = MZILCH( I,NCSP ) 
            MZHI0( I,NCSP ) = MZ - 4  
            JZ3 = 0 

            DO JZ = 1, MZHI0( I,NCSP ), 5  
               JZ3 = JZ + 4 
            END DO

            MZLO1( I,NCSP ) = JZ3 + 1
            MZHI1( I,NCSP ) = MZ  - 1
            JZ4 = JZ3 

            DO JZ = JZ3+1, MZ-1, 2 
               JZ4 = JZ + 1 
            END DO

            MZLO2( I,NCSP ) = JZ4 + 1
         END DO
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Check dimensions 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( ICNT .GT. MXCOUNT2 .OR. JCNT .GT. MXCOUNT1 .OR. 
     &        IARRAY2 .GT. MXARRAY .OR. ICCOUNT .GT. MXCOUNT2 .OR.
     &        JCCOUNT .GT. MXARRAY ) THEN
            WRITE( MSG, 94000 ) 
            WRITE( 6, * ) MSG
            WRITE( MSG, 94020 ) MXCOUNT2, ICNT 
            WRITE( 6, * ) MSG
            WRITE( MSG, 94040 ) MXCOUNT1, JCNT 
            WRITE( 6, * ) MSG
            WRITE( MSG, 94060 ) MXARRAY, IARRAY2 
            WRITE( 6, * ) MSG
            WRITE( MSG, 94080 ) MXARRAY, ICCOUNT 
            WRITE( 6, * ) MSG
            WRITE( MSG, 94100 ) MXARRAY, JCCOUNT 
            WRITE( 6, * ) MSG
            STOP
         END IF           

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set final arrays for partial derivative calculations
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         DO NRX = 1, NUSERAT( NCSP )
            NK = NKUSERAT( NRX,NCSP )
            DO IAL = 1, NREACT( NK )
               IR = IRM2( NK,IAL,NCS )

               DO JAL = 1, NREACT( NK )
                  JR = IRM2( NK,JAL,NCS )
                  IAR = JARRAYPT( JR,IR,NCSP )
                  NDERIVL( NK,NCSP ) = NDERIVL( NK,NCSP ) + 1
                  NLS = NDERIVL( NK,NCSP )
                  JARRL( NK,NLS,NCSP ) = IAR
                  JLIAL( NK,NLS,NCSP ) = IAL
                  NDLMAX = MAX( NLS,NDLMAX )
               END DO
               
               DO IAP = 1, NPRDCT( NK )
                  JP = IRM2( NK,IAP+3,NCS )
                  IAR = JARRAYPT( JP,IR,NCSP )
                  NDERIVP( NK,NCSP ) = NDERIVP( NK,NCSP ) + 1
                  NPR = NDERIVP( NK,NCSP )
                  JARRP(  NK,NPR,NCSP ) = IAR
                  JPIAL(  NK,NPR,NCSP ) = IAL
                  ICOEFF( NK,NPR,NCSP ) = 0
                  IF ( ABS( SC( NK,IAP ) - 1.0 ) .GT. 1.0E-06 ) THEN
                     ICOEFF( NK,NPR,NCSP ) = IAP
                  END IF
                  NDPMAX = MAX( NPR,NDPMAX )
               END DO
            END DO     
         END DO
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Check dimensions of PD arrays
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( NDPMAX .GT. MXRP .OR. NDLMAX .GT. MXRR ) THEN
            WRITE( MSG, 94000 ) 
            WRITE( 6, * ) MSG
            WRITE( MSG, 94200 ) MXRP, NDPMAX 
            WRITE( 6, * ) MSG
            WRITE( MSG, 94220 ) MXRR, NDLMAX 
            WRITE( 6, * ) MSG
            STOP
         END IF
700   CONTINUE

      NPDERIV  = SUM( NREACT )
      
      MXIARRAY  = MAXVAL( IARRAY ) + 1

      ALLOCATE( NDERIVN1( 0:MXIARRAY, NCS2 ),    
     &           NDERIVP1( 0:MXIARRAY, NCS2 ),
     &           NDERIVCO( 0:MXIARRAY, NCS2 ), STAT = STATUS )
      IF ( STATUS .NE. 0 ) THEN
         XMSG = 'ERROR allocating NDERIVN1, NDERIVP1, NDERIVCO  '
         WRITE( 6, * )'WRT_REACTIONS_MODULE : ' // XMSG
         STOP
      END IF
      
      NDERIVN1 = 0
      NDERIVP1 = 0
      NDERIVCO = 0
      
      ALLOCATE( PDERIVN1(  NPDERIV, 0:MXIARRAY, NCS2 ),
     &           PDERIVP1( NPDERIV, 0:MXIARRAY, NCS2 ),
     &           PDERIVCO( NPDERIV, 0:MXIARRAY, NCS2 ), 
     &           PD_COEFF( NPDERIV, 0:MXIARRAY, NCS2 ), STAT = STATUS )
      IF ( STATUS .NE. 0 ) THEN
         XMSG = 'ERROR allocating PDERIVN1, PDERIVP1, PDERIVCO, PD_COEFF '
         WRITE( 6, * )'WRT_REACTIONS_MODULE : ' // XMSG
         STOP
      END IF
      
      PDERIVN1 = 0
      PDERIVP1 = 0
      PDERIVCO = 0
      PD_COEFF = 0.0D0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Find names for output module file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      OPEN ( UNIT = MODULE_UNIT, FILE = FNAME_DATA_MODULE, STATUS = 'UNKNOWN' )

      WRITE( MODULE_UNIT,'(7X,"MODULE RXNS_DATA", 3/ 7X, "IMPLICIT NONE" 3/ )')
       

      LITE = .FALSE.
      
      CALL WREXTS_FORTRAN90 ( MODULE_UNIT, EQUATIONS_MECHFILE,
     &              MECHNAME,
     &              NS, SPCLIS, SPC1RX,
     &              NR,
     &              IP,
     &              NAMCONSTS,
     &              CONST, SS1RX, LITE ) 
      

      CALL WRSPECIAL_EXT_FORTRAN90( MODULE_UNIT )

      IF( N_SS_SPC .GT. 0 ) CALL GET_SS_DATA( LUNOUT, NR ) 


      CALL WRSS_EXT_FORTRAN90( MODULE_UNIT, NR ) 

! set up variables equal to the rate constant of type 10 fall off reactions      
!      DO NXX = 1, NR
!         IF( KTYPE( NXX ) .EQ. 10 )THEN
!             WRITE(MODULE_UNIT, 4505)LABEL(NXX,1)
!         END IF
!      END DO       
! set up variables equal to the rate constant of type 11 reactions      
!      WRITE(MODULE_UNIT,4749)
      IF( NSPECIAL .GT. 0 )THEN
         WRITE(MODULE_UNIT,4750)
      ELSE
         WRITE(MODULE_UNIT,4751)
      END IF
! set up pointers and names for photolysis rate array
      WRITE(MODULE_UNIT,4502)
      DO IPR = 1, NPHOTAB
         WRITE(MODULE_UNIT,4503),PHOTAB(IPR),IPR
      END DO
!      DO IPR = 1, NPHOTAB
!         WRITE(MODULE_UNIT,4557)IPR, PHOTAB(IPR)
!      END DO
      IF( NHETERO .GT. 0 )THEN
 !         WRITE(MODULE_UNIT,5023)NHETERO
          DO IPR = 1, NHETERO
             WRITE(MODULE_UNIT,5024)HETERO(IPR),IPR
          END DO
!          DO IPR = 1, NHETERO
!             WRITE(MODULE_UNIT,5025)IPR, HETERO(IPR)
!          END DO
      ELSE 
!          WRITE(MODULE_UNIT,5026)NHETERO
      END IF

      WRITE( MODULE_UNIT,'(7X, "END MODULE RXNS_DATA")')
     
      CLOSE( MODULE_UNIT )

      OPEN ( UNIT = MODULE_UNIT, FILE = FNAME_FUNC_MODULE, STATUS = 'UNKNOWN' )

      WRITE( MODULE_UNIT,'(7X,"MODULE RXNS_FUNCTION", 3/ 7X, "IMPLICIT NONE" 3/ )')

      WRITE( MODULE_UNIT, 4601)MECHNAME
      WRITE( MODULE_UNIT, 4500) 
      
       ISPC = INDEX(EQN_MECH_KPP,'/mech', BACK= .TRUE.) + 1
       NXX  = INDEX(EQN_MECH_KPP,'.eqn', BACK= .TRUE.)  - 1
       PHRASE = MECHNAME
       CALL CONVERT_CASE ( PHRASE, .FALSE. )
       WRITE(MODULE_UNIT,95050)
       DO ISPC = 1, NSPECIAL
          WRITE(MODULE_UNIT, 4506)SPECIAL(ISPC)
       END DO
       WRITE(MODULE_UNIT,95051)

       DO NXX = 1, NSPECIAL
! count total number of terms in special rates
         SPECIAL_TERMS = 0
         DO IREACT = 1, MAXSPECTERMS
            IF( INDEX_KTERM( NXX, IREACT ) .GT. 0 )THEN
                SPECIAL_TERMS = SPECIAL_TERMS + 1
            END IF
            IF( OPERATORS( NXX, IREACT ) .GT. 0 )THEN
                SPECIAL_TERMS = SPECIAL_TERMS + 1
            END IF
         END DO         
         WRITE(MODULE_UNIT,'(11X, A16)', ADVANCE = 'NO' )SPECIAL( NXX )
         FIRST_TERM = .TRUE.
! first write standard rate constants time concentrations
         COUNT_TERMS = 0
         DO IREACT = 1, MAXSPECTERMS
             IRX  = INDEX_KTERM( NXX, IREACT )
             IF( IRX .LT. 1 )CYCLE
             COUNT_TERMS = COUNT_TERMS + 1 
             IF( FIRST_TERM )THEN
                PHRASE = ' = '
                FIRST_TERM = .FALSE.
                IF(KC_COEFFS( NXX, IREACT ) .LT. 0.0 )PHRASE = ' = ' // ' - '
             ELSE
                WRITE(MODULE_UNIT, 4711, ADVANCE = 'NO' )
                PHRASE = ' + '
                IF(KC_COEFFS( NXX, IREACT ) .LT. 0.0 )PHRASE = ' - '
             END IF
             IF( KC_COEFFS( NXX, IREACT ) .NE. 1.0 )THEN
                 WRITE(MODULE_UNIT, 4708, ADVANCE = 'NO')TRIM(PHRASE),
     &          REAL( ABS( KC_COEFFS( NXX, IREACT ) ), 8), IRX
             ELSE
                 WRITE(MODULE_UNIT, 4706, ADVANCE = 'NO')TRIM(PHRASE),IRX
             END IF
             ISPC = INDEX_CTERM( NXX, IREACT )
             IF( ISPC .LT. 1 )CYCLE
!             WRITE(PHRASE,'(A,I4,A)')' * Y( NCELL, ', IOLD2NEW(ISPC,NCS) , ' ) '
             WRITE(PHRASE,'(A,I4,A)')' * Y( NCELL, IOLD2NEW( ', ISPC, ', NCS) ) '
             WRITE(MODULE_UNIT, 4709, ADVANCE = 'NO')TRIM( PHRASE )
             IF( IREACT .LT. MAXSPECTERMS )THEN
                 IF( COUNT_TERMS .LT. SPECIAL_TERMS )THEN
                     WRITE(MODULE_UNIT, 75006, ADVANCE = 'NO')
                 END IF
             END IF
         END DO
! next write defined operators         
         DO IREACT = 1, MAXSPECTERMS
            IDX = OPERATORS( NXX, IREACT )
            IF( IDX .LT. 1 )CYCLE
             COUNT_TERMS = COUNT_TERMS + 1 
             IF( FIRST_TERM )THEN
                PHRASE = ' = '
                IF(OPERATOR_COEFFS( NXX, IREACT ) .LT. 0.0 )PHRASE = ' = ' // ' - '
                FIRST_TERM = .FALSE.
             ELSE
                WRITE(MODULE_UNIT, 4711, ADVANCE = 'NO' )
                PHRASE = ' + '
                IF(OPERATOR_COEFFS( NXX, IREACT ) .LT. 0.0 )PHRASE = ' - '
             END IF
             IF( OPERATOR_COEFFS( NXX, IREACT ) .NE. 1.0 )THEN
                 WRITE(MODULE_UNIT, 4710, ADVANCE = 'NO')TRIM(PHRASE),
     &           REAL( ABS( OPERATOR_COEFFS( NXX, IREACT ) ), 8), TRIM( SPECIAL( IDX ) )
             ELSE
                 WRITE(MODULE_UNIT, 4712, ADVANCE = 'NO')TRIM(PHRASE),TRIM( SPECIAL( IDX ) )
             END IF
             IF( IREACT .LT. MAXSPECTERMS )THEN
                 IF( COUNT_TERMS .LT. SPECIAL_TERMS )THEN
                     WRITE(MODULE_UNIT, 75006, ADVANCE = 'NO')
                 END IF
             END IF
         END DO 
         WRITE(MODULE_UNIT, * )' '
      END DO
75006 FORMAT(2X, "&")      
      WRITE(MODULE_UNIT,95701)
95701 FORMAT(/ '! define rate constants in terms of special rate operators ' /)
      DO NXX = 1, NSPECIAL_RXN 
         WRITE(MODULE_UNIT,95070)ISPECIAL( NXX,1 ),SPECIAL( ISPECIAL( NXX,2 ) ),
     &   TRIM( LABEL( ISPECIAL( NXX,1 ),1 ) )
      END DO
95070 FORMAT(11X,'RKI( NCELL,',I4,' ) = ',A16,' ! reaction: ',A)
      WRITE(MODULE_UNIT,95060)
      WRITE(MODULE_UNIT,4504)

!      WRITE(MODULE_UNIT,'(A)')' '
! initialize special rate expressions in Update_Rconst subroutine
!      DO NXX = 1, NSPECIAL
!         WRITE(MODULE_UNIT,95100)SPECIAL( NXX ) 
!      END DO

      WRITE(MODULE_UNIT,99890)
      
      IF( KUNITS .EQ. 2 )THEN
          WRITE(MODULE_UNIT,'(3A)')'! All rate constants converted from  molec/cm3 to ppm'
          WRITE(MODULE_UNIT,'(3A)')'! and 1/sec to 1/min'
      ELSE
          WRITE(MODULE_UNIT,'(3A)')'! Only fall off rate constants converted from  molec/cm3 '
          WRITE(MODULE_UNIT,'(3A)')'! and 1/sec to 1/min'
          WRITE(MODULE_UNIT,'(3A)')'! Remainder use ppm and 1/min '
      END IF
      
      IF( IP .GT. 0 )THEN
          WRITE(MODULE_UNIT,5115)
5115   FORMAT(/ 13X,'IF(  LSUNLIGHT )THEN ! set photolysis rates ')
          DO IPR = 1, IP
             NXX = IPH( IPR,1 )
             IF( NXX .LE. 0 )CYCLE
             WRITE(MODULE_UNIT, 5117, ADVANCE= 'NO')LABEL(NXX,1), NXX
             IF ( IPH( IPR,3 ) .NE. 0 )THEN
                IDX = IPH( IPR, 2 )
                IF( RTDAT(1, NXX) .NE. 1.0 )THEN
                   WRITE(MODULE_UNIT,5000, ADVANCE = 'NO')REAL(RTDAT(1, NXX),8),TRIM( PHOTAB(IDX) )
                ELSE
                   WRITE(MODULE_UNIT,5001, ADVANCE = 'NO')TRIM( PHOTAB(IDX) )
                END IF
             ELSE IF( IPH( NXX,3 ) .EQ. 0 )THEN
                IDX = IPH(IPH( NXX,2 ), 2)
                IF( RTDAT(1, NXX) .NE. 1.0 )THEN
                   WRITE(MODULE_UNIT,5100, ADVANCE = 'NO')REAL(RTDAT(1, NXX),8), IDX
                ELSE
                   WRITE(MODULE_UNIT,5101, ADVANCE = 'NO')IDX
                END IF
             END IF
         END DO
          WRITE(MODULE_UNIT,5116)
5116   FORMAT(/13X,'END IF' /)
5117  FORMAT(/    '!  Reaction Label ', A / 16X, 'RKI( NCELL, ', I4, ') = ')
      END IF
      
      DO NXX = 1, NR

         IF( KTYPE( NXX ) .NE. 11 .AND. KTYPE( NXX ) .NE. 0 )THEN
!          WRITE(MODULE_UNIT, 1498 )TRIM(LABEL(NXX,1))
!            CYCLE
!         ELSE
            WRITE(MODULE_UNIT, 1501, ADVANCE= 'NO')LABEL(NXX,1), NXX
         END IF 
         
         SELECT CASE( KTYPE( NXX ) )
          CASE( -1 )
!             IF( KUNITS .EQ. 2 )CALL WRITE_RATE_CONVERT_TIME(MODULE_UNIT, IORDER(NXX))
             DO IPR = 1, NHETERO
                IF ( IHETERO( IPR,1 ) .EQ. NXX )EXIT
             END DO
             IDX = IHETERO( IPR, 2 )
             IF( RTDAT(1, NXX) .NE. 1.0 )THEN
                 WRITE(MODULE_UNIT,5027, ADVANCE = 'NO')REAL(RTDAT(1, NXX),8),TRIM( HETERO(IDX) )
                 PRINT*,REAL(RTDAT(1, NXX),8),TRIM( HETERO(IDX) )
             ELSE
                 WRITE(MODULE_UNIT,5128, ADVANCE = 'NO')TRIM( HETERO(IDX) )
                 PRINT*,TRIM( HETERO(IDX) )
                 WRITE(6,5028)TRIM( HETERO(IDX) )
             END IF
!          CASE(  0 )
!             DO IPR = 1, IP
!                IF ( IPH( IPR,1 ) .EQ. NXX )EXIT
!             END DO
!             IF ( IPH( IPR,3 ) .NE. 0 )THEN
!                IDX = IPH( IPR, 2 )
!                IF( RTDAT(1, NXX) .NE. 1.0 )THEN
!                   WRITE(MODULE_UNIT,5000, ADVANCE = 'NO')REAL(RTDAT(1, NXX),8),TRIM( PHOTAB(IDX) )
!                ELSE
!                   WRITE(MODULE_UNIT,5001, ADVANCE = 'NO')TRIM( PHOTAB(IDX) )
!                END IF
!             ELSE IF( IPH( NXX,3 ) .EQ. 0 )THEN
!                IDX = IPH(IPH( NXX,2 ), 2)
!                IF( RTDAT(1, NXX) .NE. 1.0 )THEN
!                   WRITE(MODULE_UNIT,5100, ADVANCE = 'NO')REAL(RTDAT(1, NXX),8), IDX
!                ELSE
!                   WRITE(MODULE_UNIT,5101, ADVANCE = 'NO')IDX
!                END IF
!             END IF
          CASE( 1 )
             IF( KUNITS .EQ. 2 )CALL WRITE_RATE_CONVERT_LOCAL(MODULE_UNIT, IORDER(NXX))
             WRITE(MODULE_UNIT,'(1PD12.4)', ADVANCE = 'NO')REAL(RTDAT(1, NXX), 8)
          CASE( 2 )
             IF( KUNITS .EQ. 2 )CALL WRITE_RATE_CONVERT_LOCAL(MODULE_UNIT, IORDER(NXX))
             WRITE(MODULE_UNIT,5029, ADVANCE = 'NO')RTDAT(1, NXX), RTDAT(2, NXX)
          CASE( 3 )
             IF( KUNITS .EQ. 2 )CALL WRITE_RATE_CONVERT_LOCAL(MODULE_UNIT, IORDER(NXX))
             WRITE(MODULE_UNIT,5003, ADVANCE = 'NO')RTDAT(1, NXX), RTDAT(3, NXX)
          CASE( 4 )
             IF( KUNITS .EQ. 2 )CALL WRITE_RATE_CONVERT_LOCAL(MODULE_UNIT, IORDER(NXX))
             WRITE(MODULE_UNIT,5004, ADVANCE = 'NO')RTDAT(1, NXX), RTDAT(3, NXX), RTDAT(2, NXX)
          CASE( 5 )
             IRX = INT( RTDAT( 3, NXX) )
             WRITE(MODULE_UNIT,5005, ADVANCE = 'NO')IRX,RTDAT( 1, NXX ), RTDAT(2, NXX )
          CASE( 6 )
             IRX = INT( RTDAT( 2, NXX) )
             IF( RTDAT( 1, NXX ) .NE. 1.0 )THEN
                 WRITE(MODULE_UNIT, 5006, ADVANCE = 'NO')REAL(RTDAT( 1, NXX ), 8), IRX
             ELSE
                 WRITE(MODULE_UNIT, 4706, ADVANCE = 'NO')' ', IRX
             END IF
          CASE( 7 )
             IF( RTDAT(2, NXX) .NE. 0.0 )THEN
                 WRITE(MODULE_UNIT,5014, ADVANCE = 'NO')REAL(RTDAT(1, NXX), 8),REAL(RTDAT(2, NXX), 8)
             ELSE
                 WRITE(MODULE_UNIT,5007, ADVANCE = 'NO')REAL(RTDAT(1, NXX), 8)
             END IF
          CASE( 8 )
             DO IDX = 1, NFALLOFF
                IF( IRRFALL( IDX ) .EQ. NXX )EXIT
             END DO
             CALL WRITE_RATE_CONVERT_LOCAL(MODULE_UNIT, IORDER(NXX))
             WRITE(MODULE_UNIT,5008, ADVANCE = 'NO')RTDAT(1,NXX),(1.0*RTDAT(2,NXX)),RTDAT(3,NXX),
     &      (1.0*RFDAT(1,IDX)),RFDAT(2,IDX),(1.0*RFDAT(3,IDX))
          CASE( 9 )
             DO IDX = 1, NFALLOFF
                IF( IRRFALL( IDX ) .EQ. NXX )EXIT
             END DO
             CALL WRITE_RATE_CONVERT_LOCAL(MODULE_UNIT, IORDER(NXX))
             IF( RFDAT( 2, IDX ) .EQ. 0.0 .AND. RFDAT( 3, IDX ) .EQ. 0.0 )THEN
                 WRITE(MODULE_UNIT,5009, ADVANCE = 'NO')RTDAT(1,NXX),RTDAT(2,NXX),
     &           RTDAT(3,NXX),1.0*RFDAT(1,IDX)
             ELSE 
                 WRITE(MODULE_UNIT,5019, ADVANCE = 'NO')RTDAT(1,NXX),RFDAT(2, IDX),RTDAT(2,NXX),
     &           RTDAT(3,NXX),RFDAT(3, IDX),1.0*RFDAT(1,IDX),RFDAT(4, IDX),RFDAT(5, IDX)
              END IF 
          CASE( 10 )
             DO IDX = 1, NFALLOFF
                IF( IRRFALL( IDX ) .EQ. NXX )EXIT
             END DO
             CALL WRITE_RATE_CONVERT_LOCAL(MODULE_UNIT, IORDER(NXX))
             WRITE(MODULE_UNIT, 5010, ADVANCE = 'NO')RTDAT(1,NXX),RTDAT(3,NXX),RTDAT(2,NXX),
     &      RFDAT(1,IDX),RFDAT(3,IDX),RFDAT(2,IDX),RFDAT(5,IDX),RFDAT(4,IDX)
          CASE( 11 )
             WRITE(MODULE_UNIT, 1498 )TRIM(LABEL(NXX,1))
!             DO IDX = 1, NSPECIAL_RXN
!                IF( ISPECIAL( IDX, 1 ) .EQ. NXX )EXIT
!             END DO
!             I   = ISPECIAL( IDX, 1)
!             IRX = ISPECIAL( IDX, 2)
!             IF( RTDAT( 1, I) .NE. 1.0 )THEN
!                WRITE(MODULE_UNIT,5011, ADVANCE = 'NO')REAL(RTDAT( 1, I),8), TRIM( SPECIAL( IRX ) )
!             ELSE
!                WRITE(MODULE_UNIT,5012, ADVANCE = 'NO')TRIM( SPECIAL( IRX ) )
!             END IF
          END SELECT
!          WRITE( MODULE_UNIT,'(/)')
      END DO

      WRITE(MODULE_UNIT,99991)


      WRITE( MODULE_UNIT,'(7X,"END MODULE RXNS_FUNCTION")')
      CLOSE( MODULE_UNIT )

1498  FORMAT(/ '! RKI for Reaction ', A,' set in SPECIAL_RATES Routine' )

1501  FORMAT(/    '!  Reaction Label ', A / 13X, 'RKI( NCELL, ', I4, ') = ')
1993  FORMAT( / 5X, '*** ERROR: Special label already used'
     &        / 5X, 'Processing for special label number:', I6 )
1994  FORMAT( / 5X, '*** ERROR: Equal sign expected after special label'
     &        / 5X, 'Last line read was:' / A81 )
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
     &        / 5X, 'Processing for reaction number:', I6, 1X, A )
2021  FORMAT( / 5X, '*** ERROR: Label points to currently undefined reaction'
     &        / 5X, 'Processing for reaction number:', I6 )
2031  FORMAT( / 5X, '*** ERROR: Special Rate Coefficient ', A16,
     &              ' uses the unlisted reaction label ', A16 )
2032  FORMAT( / 5X, '*** ERROR: Special Rate Coefficient ', A16,
     &              ' incorrectly uses the reaction ', A16,'.',
     &              ' The reaction order is misinterpreted as 1st or 2nd')
2033  FORMAT( / 5X, '*** ERROR: Special Rate Coefficient ', A16,
     &              ' uses the unlisted species ', A16 )
2034  FORMAT( / 5X, '*** ERROR: Special Rate Coefficient ', A16,
     &              ' incorrectly uses the reaction ', A16,'.',
     &              ' The reaction order is not 2nd.')

3010  FORMAT( / 5X, '*** ERROR: The following steady-state species is also in the ',
     &              'ELIMINATE list' )
3011  FORMAT( 16X, A )

4001  FORMAT( / 5X, '*** ERROR: Number of Steady-state species exceeds max allowable;',
     &              ' increase MAXNLIST' )

4002  FORMAT( / 5X, '*** ERROR: Number of ELIMINATE species exceeds max allowable;',
     &              ' increase MAXNLIST' )
4505  FORMAT('REAL( 8 )  :: RKI_RXN_', A16,' ! rate constant for stated reaction label')        
4506  FORMAT( 7X, 'REAL( 8 )  :: ', A16)        

4500  FORMAT(/7X,'CONTAINS'
     &      2/ 7X,'REAL( 8 ) FUNCTION FALL_T10 ( A0,B0,C0,A1,B1,C1,CE,CF)'
     &      / 9X,'IMPLICIT NONE'
     &      / '! rate constant for CMAQ fall off reaction type 10'
     &      / '! Arguements:'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: A0'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: B0'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: C0'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: A1'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: B1'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: C1'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: CE'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: CF'
     &      / 9X,'! Local:'
     &      / 9X,'REAL( 8 ) K0'
     &      / 9X,'REAL( 8 ) K1'
     &      / 9X,'REAL( 8 ) KEND'
     &      / 9X,'K0 = A0 * CAIR * DEXP(B0*INV_TEMP)* TEMPOT300**C0'
     &      / 9X,'K1 = A1 * DEXP(B1*INV_TEMP) * TEMPOT300**C1'
     &      / 9X,'KEND = ( ( 1.0D0 + ( ( 1.0D0 / CE ) * DLOG10( K0 / K1 ) ) ** 2.0D0 ) )'
     &      / 9X,'KEND = 1.0D0 / KEND'
     &      / 9X,'FALL_T10 = ( K0 / ( 1.0D0 + K0/K1 ) ) * CF ** KEND'
     &      / 9X,'RETURN'
     &      / 7X,'END FUNCTION FALL_T10' 
     &      / 7X,'REAL( 8 ) FUNCTION POWE_T02( A0,B0 )'
     &      / 9X,'IMPLICIT NONE'
     &      / '! rate constant for CMAQ Arrhenius reaction type 2'
     &      / '! Arguements:'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: A0'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: B0'
     &      / 9X,'! Local: None'
     &      / 9X,'POWE_T02 =  A0 * TEMPOT300**B0'
     &      / 9X,'RETURN'
     &      / 7X,'END FUNCTION POWE_T02'
     &      / 7X,'REAL( 8 ) FUNCTION ARRE_T04( A0,B0,C0 )'
     &      / 9X,'IMPLICIT NONE'
     &      / '! rate constant for CMAQ Arrhenius reaction type 4'
     &      / '! Arguements:'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: A0'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: B0'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: C0'
     &      / 9X,'! Local:'
     &      / 9X,'INTRINSIC DEXP'
     &      / 9X,'ARRE_T04 =  A0 * DEXP( B0 * INV_TEMP ) * TEMPOT300**C0'
     &      / 9X,'RETURN'
     &      / 7X,'END FUNCTION ARRE_T04'
     &      / 7X,'REAL( 8 ) FUNCTION ARRE_T03( A0,B0 )'
     &      / '! rate constant for CMAQ Arrhenius reaction type 3'
     &      / 9X,'IMPLICIT NONE'
     &      / '! Arguements:'
     &      / 9X,'REAL( 8 ),     INTENT(IN) ::  A0'
     &      / 9X,'REAL( 8 ),     INTENT(IN) ::  B0'
     &      / 9X,'! Local:'
     &      / 9X,'INTRINSIC DEXP'
     &      / 9X,'ARRE_T03 =  A0 * DEXP( B0 * INV_TEMP )'
     &      / 9X,'RETURN'
     &      / 7X,'END FUNCTION ARRE_T03 '
     &      / 7X,'REAL( 8 ) FUNCTION FALL_T08(A0,C0,A2,C2,A3,C3)'
     &      / '! rate constant for CMAQ fall off reaction type 8'
     &      / 9X,'IMPLICIT NONE'
     &      / '! Arguements:'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: A0'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: C0'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: A2'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: C2'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: A3'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: C3'
     &      / 9X,'! Local:'
     &      / 9X,'REAL( 8 ) K0'
     &      / 9X,'REAL( 8 ) K2'
     &      / 9X,'REAL( 8 ) K3'
     &      / 9X,'INTRINSIC DEXP'
     &      / 9X,'K0 = A0 * DEXP( C0 * INV_TEMP )'
     &      / 9X,'K2 = A2 * DEXP( C2 * INV_TEMP )'
     &      / 9X,'K3 = A3 * DEXP( C3 * INV_TEMP )'
     &      / 9X,'K3 = K3 * CAIR'
     &      / 9X,'FALL_T08 = K0 + K3/( 1.0D0 + K3/K2 )'
     &      / 9X,'RETURN'
     &     2/ 7X,'END FUNCTION FALL_T08'
     &     2/ 7X,'REAL( 8 ) FUNCTION FALL_T11(A1,B1,C1,A2, B2, C2)'
     &      / '! rate constant for CMAQ fall off reaction type 11'
     &      / '! actually expanded form of type 9'
     &      / 9X,'IMPLICIT NONE'
     &      / '! Arguements:'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: A1'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: B1'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: C1'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: A2'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: B2'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: C2'
     &      /9X,'!  Local:'
     &      / 9X,'REAL( 8 ) K1'
     &      / 9X,'REAL( 8 ) K2'
     &      / 9X,'INTRINSIC DEXP'
     &      / 9X,'K1 = A1 * DEXP( C1 * INV_TEMP ) * TEMPOT300**B1'
     &      / 9X,'K2 = A2 * DEXP( C2 * INV_TEMP ) * TEMPOT300**B2'
     &      / 9X,'FALL_T11 = K1 + K2 * CAIR'
     &      / 9X,'RETURN'
     &     2/ 7X,'END FUNCTION FALL_T11'     
     &     2/ 7X,'REAL( 8 ) FUNCTION FALL_T09(A1,C1,A2,C2)'
     &      / '! rate constant for CMAQ fall off reaction type 9'
     &      / 9X,'IMPLICIT NONE'
     &      / '! Arguements:'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: A1'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: C1'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: A2'
     &      / 9X,'REAL( 8 ), INTENT( IN ) :: C2'
     &      /9X,'!  Local:'
     &      / 9X,'REAL( 8 ) K1'
     &      / 9X,'REAL( 8 ) K2'
     &      / 9X,'INTRINSIC DEXP'
     &      / 9X,'K1 = A1 * DEXP( C1 * INV_TEMP )'
     &      / 9X,'K2 = A2 * DEXP( C2 * INV_TEMP )'
     &      / 9X,'FALL_T09 = K1 + K2 * CAIR'
     &      / 9X,'RETURN'
     &     2/ 7X,'END FUNCTION FALL_T09'       )
    
4501   FORMAT( '! Name of Mechanism ', A
     &        / 7X,'PUBLIC             :: CALC_RCONST, SPECIAL_RATES'
     &       2/ 7X,'REAL( 8 ), PRIVATE :: CAIR          ! air number density (wet) [molec/cm^3]'
     &        / 7X,'REAL( 8 ), PRIVATE :: CFACT         ! molec/cc to ppm conversion factor   '/
     &        / 7X,'REAL( 8 ), PRIVATE :: CFACT_SQU     ! molec/cc to ppm conversion factor squared  '/
     &        / 7X,'REAL( 8 ), PRIVATE :: INV_CFACT     ! Reciprocal of molec/cc to ppm conversion factor   '/
     &        / 7X,'REAL( 8 ), PRIVATE :: TEMPOT300     ! temperature divided by 300 K, dimensionaless '/
     &        / 7X,'REAL( 8 ), PRIVATE :: INV_TEMP      ! reciprocal of air temperature, K-1'
     &        / 7X,'REAL( 8 ), PRIVATE :: TEMP          ! air temperature, K'
     &        / 7X,'REAL( 8 ), PRIVATE :: PRESS         ! pressure [Atm] '
     &        / 7X,'REAL( 8 ), PRIVATE :: INV_RFACTOR   ! Convertor: ppm/min to molec/(cm^3*sec)'
     &        / 7X,'REAL( 8 ), PRIVATE :: RFACTOR_SQU   ! Convertor cm^6/(molec^2*sec) to 1/(ppm^2*min)'
     &        / 7X,'REAL( 8 ), PRIVATE :: RFACTOR       ! Convertor cm^3/(molec*sec) to 1/(ppm*min)'
     &        / 7X,'REAL,      PRIVATE :: H2O           ! Cell H2O mixing ratio (ppmV)')

4601   FORMAT( '! Name of Mechanism ', A
     &        / 7X,'PUBLIC             :: CALC_RCONST, SPECIAL_RATES'
     &       2/ 7X,'REAL( 8 ), PRIVATE :: CAIR          ! air number density (wet) [molec/cm^3]'
     &        / 7X,'REAL( 8 ), PRIVATE :: CFACT         ! Convertor cm^3/(molec*sec) to 1/(ppm*min)'/
     &        / 7X,'REAL( 8 ), PRIVATE :: CFACT_SQU     ! Convertor cm^6/(molec^2*sec) to 1/(ppm^2*min)'/
     &        / 7X,'REAL( 8 ), PRIVATE :: INV_CFACT     ! ppm/min to molec/(cm^3*sec)'/
     &        / 7X,'REAL( 8 ), PRIVATE :: TEMPOT300     ! temperature divided by 300 K, dimensionaless '/
     &        / 7X,'REAL( 8 ), PRIVATE :: INV_TEMP      ! reciprocal of air temperature, K-1'
     &        / 7X,'REAL( 8 ), PRIVATE :: TEMP          ! air temperature, K'
     &        / 7X,'REAL( 8 ), PRIVATE :: PRESS         ! pressure [Atm] '
     &        / 7X,'REAL( 8 ), PRIVATE :: INV_RFACT     ! ppm/min to molec/(cm^3*min)'
     &        / 7X,'REAL( 8 ), PRIVATE :: RFACT_SQU     ! cm^6/(molec^2*min) to 1/(ppm^2*min)'
     &        / 7X,'REAL( 8 ), PRIVATE :: RFACT         ! cm^3/(molec*min) to 1/(ppm*min)'
     &        / 7X,'REAL,      PRIVATE :: H2O           ! Cell H2O mixing ratio (ppmV)')
     
4502   FORMAT(  '! pointers and names to specific photolysis rates' )
4503   FORMAT(  7X,'INTEGER, PARAMETER  :: IJ_',A16,' = ', I3 )
4504   FORMAT(' ' )
4555   FORMAT(' ')
4556   FORMAT( 'RFACTOR       = 6.0D-5  * CAIR'
     &       / 'INV_RFACTOR   = 6.0D+7  / CAIR'
     &       / 'RFACTOR_SQU   = 6.0D-11 * CAIR * CAIR'
     &       / 'CFACTOR       = 1.0D0'
     &       / 'INV_TEMP      = 1.0D0 / TEMP'
     &       / 'COEFF_FALLOFF = CAIR ' )
4557   FORMAT('DATA PHOTAB(', I3,' ) / ''',A16,''' /')
4507  FORMAT('RKI_RXN_', A16,A4)        
4706  FORMAT(A,1X,'RKI( NCELL, ', I4,' ) ')
4708  FORMAT(A,1X,1PD12.4,' * RKI( NCELL, ', I4,' ) ')
4709  FORMAT( A )     
4710  FORMAT(A,1X,1PD12.4,' * ', A)
4711  FORMAT( / 5X, '&' 21X)
4712  FORMAT(A, 1X, A)
4713  FORMAT( '!If( .Not. CALC_RCONST )Then'
     &      / '!   Return'
     &      / '!Else'
     &      / '!   CALC_RCONST = .False.'
     &      / '!End If' 
     &      / '! Rate Constant Units produce reaction rates in ppm/min' )
4714  FORMAT('! number mixing ratios of constant atmospheric species, ppmV')     
4749   FORMAT('!Flag to call SPECIAL_RATES rountine in Integrator ')
4750   FORMAT(7X, 'LOGICAL,  PARAMETER :: USE_SPECIAL_RATES = .TRUE. ')
4751   FORMAT(7X, 'LOGICAL,  PARAMETER :: USE_SPECIAL_RATES = .FALSE.')
5000   FORMAT(1PD12.4,' * RJBLK( NCELL, IJ_',A,' )')
5001   FORMAT( 1X, 'RJBLK( NCELL, IJ_',A, ' )' )
5100   FORMAT(1PD12.4,' * RKI( NCELL, ',I4,' )')
5101   FORMAT(  'RKI( NCELL, ',I4,' )')
5002   FORMAT('ARRE_T04( ',1PD12.4,', 0.0000D+0,', 1PD12.4,' )')
5003   FORMAT('ARRE_T03( ',1PD12.4,', ', 1PD12.4,' )')
5004   FORMAT('ARRE_T04( ', 1PD12.4,', ', 1PD12.4,', ', 1PD12.4,' )')
5005   FORMAT('RKI( NCELL, ' I4, ' ) / ARR2( ',1PD12.4,', ',1PD12.4,' )')            
5006   FORMAT(1PD12.4,' * RKI( NCELL, ' I4, ' ) ')   
5007   FORMAT(1PD12.4,' *( 1.0D0 + 0.6D0 * PRESS )')             
5008   FORMAT('FALL_T08( ', 3(1PD12.4,', '), ' & ' / 5X, '&', 47X, 2(1PD12.4,', '), 1PD12.4, ' )' )
5009   FORMAT('FALL_T09( ', 3(1PD12.4,', '), ' & ' / 5X, '&', 47X, 1PD12.4, ' )' )
5010   FORMAT('FALL_T10( ', 3(1PD12.4,', '), ' & ' / 5X,'&', 47X, 3(1PD12.4,', '),  ' & '
     &        / 5X, '&', 47X, 1PD12.4,', ', 1PD12.4,' )')
5011   FORMAT(1PD12.4,' * ',A)             
5012   FORMAT(A)
5014   FORMAT('ARRE_T04( ',1PD12.4,', 0.0000D+0,', 1PD12.4,' )  * PRESS ')             
5019   FORMAT('FALL_T11( ', 3(1PD12.4,', ') / 5X,'&', 47X,  3(1PD12.4,', ')
     &                   / 5X,'&', 47X,  1PD12.4,' )')
5027   FORMAT(1PD12.4,' * KHETERO( NCELL, IK_',A,' )')
5028   FORMAT( 1X, 'KHETERO( NCELL, IK_',A, ' )' )
5128   FORMAT( 1X, 'BLKHET(  NCELL, IK_',A, ' )' )
5023   FORMAT(
     &        / 'INTEGER, PARAMETER  :: NHETERO  = ', I3,'  ! number of heterogeneous rates '
     &        / 'CHARACTER(16), SAVE :: HETERO(  NHETERO )  ! Names of  heterogeneous '
     &        / 'REAL( 8 )           :: KHETERO( NHETERO )  ! grid cell heterogeneous rates ,[min-1]')
5024   FORMAT(7X,'INTEGER, PARAMETER  :: IK_',A16,' = ', I3 )
5025   FORMAT('DATA HETERO(', I3,' ) / ''',A16,''' /')
5026   FORMAT(7X,'INTEGER, PARAMETER  :: NHETERO  = ', I3,'  ! number of heterogeneous rates ')
5029   FORMAT('POWE_T02( ',1PD12.4,', ', 1PD12.4,' )')
94000 FORMAT( 1X,'One of the dimensions below is too small:')
94020 FORMAT( 1X,'DIMENSION: MXCOUNT2 = ',I6,' VARIABLE: ICNT    = ',I6)  
94040 FORMAT( 1X,'DIMENSION: MXCOUNT1 = ',I6,' VARIABLE: JCNT    = ',I6)  
94060 FORMAT( 1X,'DIMENSION: MXARRAY  = ',I6,' VARIABLE: IARRAY2 = ',I6)  
94080 FORMAT( 1X,'DIMENSION: MXARRAY  = ',I6,' VARIABLE: ICCOUNT = ',I6)  
94100 FORMAT( 1X,'DIMENSION: MXARRAY  = ',I6,' VARIABLE: JCCOUNT = ',I6)
94200 FORMAT( 1X,'DIMENSION: MXRP     = ',I6,' VARIABLE: NDPMAX  = ',I6)
94220 FORMAT( 1X,'DIMENSION: MXRR     = ',I6,' VARIABLE: NDLMAX  = ',I6)

95050  FORMAT( 7X,'SUBROUTINE SPECIAL_RATES( NUMCELLS, IOLD2NEW, NCS, Y, RKI )'
     &       /  '! Purpose: calculate special rate operators and update'
     &       /  '!         appropriate rate constants'
     &      //  7X,'USE RXNS_DATA'
     &      /   7X,'IMPLICIT NONE'
     &      //  '! Arguments:'
     &      /   7X,'INTEGER,      INTENT( IN  )   :: NUMCELLS        ! Number of cells in block '
     &      /   7X,'INTEGER,      INTENT( IN  )   :: IOLD2NEW( :,: ) ! species map'
     &      /   7X,'INTEGER,      INTENT( IN  )   :: NCS             ! index for which reaction set'
     &      /   7X,'REAL( 8 ),    INTENT( IN )    :: Y( :, : )       ! species concs'
     &      /   7X,'REAL( 8 ),    INTENT( INOUT ) :: RKI( :, : )     ! reaction rate constant, ppm/min '
     &      /   '! Local:'
     &      /   7X,'INTEGER  NCELL'
     &      /   '! special rate operators listed below' //)
     
     
95051  FORMAT(/ 7X,'DO NCELL = 1, NUMCELLS' 
     &      //  '! define special rate operators' / )
95060  FORMAT(  7X,'END DO'
     &         // 7X,'RETURN'
     &          / 7X,'END SUBROUTINE SPECIAL_RATES')
95100  FORMAT(2X,A16,' = 0.0D0')        


99880 FORMAT(7X,'SUBROUTINE CALC_RCONST( BLKTEMP, BLKPRES, BLKH2O, RJBLK, LSUNLIGHT, RKI, NUMCELLS )' //
     & '!**********************************************************************' //
     & '!  Function: To compute thermal and photolytic reaction rate' /
     & '!            coefficients for each reaction.' //
     & '!  Preconditions: Photolysis rates for individual species must have' /
     & '!                 been calculated and stored in RJPHOT. Expects' /
     & '!                 temperature in deg K, pressure in atm., water' /
     & '!                 vapor in ppmV, and J-values in /min.' /
     & '!  Key Subroutines/Functions Called: None ' /
     & '!***********************************************************************'///
     &      //  7X,'USE RXNS_DATA'
     &       /  7X,'USE AEROSOL_CHEMISTRY     ! rates for heterogeneous reactions' //
     & '        IMPLICIT NONE  ' //
     & '!  Arguements: None ' //
     & '        REAL( 8 ), INTENT( IN  ) :: BLKTEMP( : )      ! temperature, deg K '/
     & '        REAL( 8 ), INTENT( IN  ) :: BLKPRES( : )      ! Reciprocal of temperature, Pa '/
     & '        REAL( 8 ), INTENT( IN  ) :: BLKH2O ( : )      ! water mixing ratio, ppm '/
     & '        REAL( 8 ), INTENT( IN  ) :: RJBLK  ( :, : )   ! photolysis rates, 1/min '/ 
     & '        INTEGER,   INTENT( IN  ) :: NUMCELLS          ! Number of cells in block ' /
     & '        REAL( 8 ), INTENT( OUT ) :: RKI ( :, : )   ! reaction rate constant, ppm/min '/
     & '        LOGICAL,   INTENT( IN  ) :: LSUNLIGHT         ! Is there sunlight? ' /
     & '!..Parameters: ' //
     & '        REAL( 8 ), PARAMETER :: COEF1  = 7.33981D+15     ! Molec/cc to ppm conv factor ' /
     & '        REAL( 8 ), PARAMETER :: CONSTC = 0.6D+0          ! Constant for reaction type 7' /
     & '        REAL( 8 ), PARAMETER :: TI300  = 1.0D+0/300.0D+0 ! reciprocal of 300 deg K' /
     & '!..External Functions: None' //
     & '!..Local Variables:' //
     & '        INTEGER NRT                  ! Loop index for reaction types '/
     & '        INTEGER IRXN                 ! Reaction number'/
     & '        INTEGER JNUM                 ! J-value species # from PHOT)'/
     & '        INTEGER KNUM                 ! Reaction # for a relative rate coeff.'/
     & '        INTEGER N                    ! Loop index for reactions'/
     & '        INTEGER NCELL                ! Loop index for # of cells in the block' 
     & //
     & '          RKI = 0.0 ' /
     & '          DO NCELL = 1, NUMCELLS ' /
     & '!  Set-up conversion factors '/
     & '             INV_TEMP  = 1.0D+00 / BLKTEMP( NCELL ) '/
     & '             CAIR      = 1.0D+06 * COEF1 * BLKPRES( NCELL ) * INV_TEMP '/
     & '             CFACT     = 6.0D-05 * CAIR' / 
     & '             CFACT_SQU = 6.0D-11 * CAIR * CAIR '/
     & '             INV_CFACT = 6.0D+07 / CAIR '/     
     & '             TEMP      = BLKTEMP( NCELL ) '/
     & '             TEMPOT300 = BLKTEMP( NCELL ) * TI300 ' // )

99991  FORMAT(7X // 7X, ' END DO  ' 
     & / '!  Multiply rate constants by [M], [O2], [N2], [H2O], [H2], or [CH4]'
     & / '!  where needed and return'
     & / 7X,'IF ( NWM .GT. 0 ) THEN'
     & / 7X,'   DO NRT = 1, NWM'
     & / 7X,'      IRXN = NRXWM( NRT )'
     & / 7X,'      DO NCELL = 1, NUMCELLS'
     & / 7X,'         RKI( NCELL,IRXN ) = RKI( NCELL,IRXN ) * ATM_AIR'
     & / 7X,'      END DO'
     & / 7X,'   END DO'
     & / 7X,'END IF' 
     & / 7X,'IF ( NWO2 .GT. 0 ) THEN'
     & / 7X,'   DO NRT = 1, NWO2'
     & / 7X,'      IRXN = NRXWO2( NRT )'
     & / 7X,'      DO NCELL = 1, NUMCELLS'
     & / 7X,'         RKI( NCELL,IRXN ) = RKI( NCELL,IRXN ) * ATM_O2'
     & / 7X,'      END DO'
     & / 7X,'   END DO'
     & / 7X,'END IF' 
     & / 7X,'IF ( NWN2 .GT. 0 ) THEN'
     & / 7X,'   DO NRT = 1, NWN2'
     & / 7X,'      IRXN = NRXWN2( NRT )'
     & / 7X,'      DO NCELL = 1, NUMCELLS'
     & / 7X,'         RKI( NCELL,IRXN ) = RKI( NCELL,IRXN ) * ATM_N2'
     & / 7X,'      END DO'
     & / 7X,'   END DO'
     & / 7X,'END IF' 
     & / 7X,'IF ( NWW .GT. 0 ) THEN'
     & / 7X,'   DO NRT = 1, NWW'
     & / 7X,'      IRXN = NRXWW( NRT )'
     & / 7X,'      DO NCELL = 1, NUMCELLS'
     & / 7X,'         RKI( NCELL,IRXN ) = RKI( NCELL,IRXN ) * BLKH2O( NCELL )'
     & / 7X,'      END DO'
     & / 7X,'   END DO'
     & / 7X,'END IF' 
     & / 7X,'IF ( NWH2 .GT. 0 ) THEN'
     & / 7X,'   DO NRT = 1, NWH2'
     & / 7X,'      IRXN = NRXWH2( NRT )'
     & / 7X,'      DO NCELL = 1, NUMCELLS'
     & / 7X,'         RKI( NCELL,IRXN ) = RKI( NCELL,IRXN ) * ATM_H2'
     & / 7X,'      END DO'
     & / 7X,'   END DO'
     & / 7X,'END IF' 
     & / 7X,'IF ( NWCH4 .GT. 0 ) THEN'
     & / 7X,'   DO NRT = 1, NWCH4'
     & / 7X,'      IRXN = NRXWCH4( NRT )'
     & / 7X,'      DO NCELL = 1, NUMCELLS'
     & / 7X,'         RKI( NCELL,IRXN ) = RKI( NCELL,IRXN ) * ATM_CH4'
     & / 7X,'      END DO'
     & / 7X,'   END DO'
     & / 7X,'END IF' 
     & / 7X, 'RETURN' 
     & / 7X,'END SUBROUTINE CALC_RCONST')

99890 FORMAT(7X,'SUBROUTINE CALC_RCONST( BLKTEMP, BLKPRES, BLKH2O, RJBLK, BLKHET, LSUNLIGHT, RKI, NUMCELLS )' //
     & '!**********************************************************************' //
     & '!  Function: To compute thermal and photolytic reaction rate' /
     & '!            coefficients for each reaction.' //
     & '!  Preconditions: Photolysis and heteorogeneous rate must have' /
     & '!                 been calculated and stored in RJPHOT. Expects' /
     & '!                 temperature in deg K, pressure in atm., water' /
     & '!                 vapor in ppmV, and J-values in /min.' /
     & '!  Key Subroutines/Functions Called: None ' /
     & '!***********************************************************************'///
     &      //  7X,'USE RXNS_DATA' //
     & '        IMPLICIT NONE  ' //
     & '!  Arguements: None ' //
     & '        REAL( 8 ), INTENT( IN  ) :: BLKTEMP( : )      ! temperature, deg K '/
     & '        REAL( 8 ), INTENT( IN  ) :: BLKPRES( : )      ! Reciprocal of temperature, Pa '/
     & '        REAL( 8 ), INTENT( IN  ) :: BLKH2O ( : )      ! water mixing ratio, ppm '/
     & '        REAL( 8 ), INTENT( IN  ) :: RJBLK  ( :, : )   ! photolysis rates, 1/min '/ 
     & '        REAL( 8 ), INTENT( IN  ) :: BLKHET ( :, : )   ! heterogeneous rate constants, ???/min'/
     & '        INTEGER,   INTENT( IN  ) :: NUMCELLS          ! Number of cells in block ' /
     & '        LOGICAL,   INTENT( IN  ) :: LSUNLIGHT         ! Is there sunlight? ' /
     & '        REAL( 8 ), INTENT( OUT ) :: RKI ( :, : )      ! reaction rate constant, ppm/min '/
     & '!..Parameters: ' //
     & '        REAL( 8 ), PARAMETER :: COEF1      = 7.33981D+15     ! Molec/cc to ppm conv factor ' /
     & '        REAL( 8 ), PARAMETER :: CONSTC     = 0.6D+0          ! Constant for reaction type 7' /
     & '        REAL( 8 ), PARAMETER :: TI300      = 1.0D+0/300.0D+0 ! reciprocal of 300 deg K' /
     & '!..External Functions: None' //
     & '!..Local Variables:' //
     & '        INTEGER NRT                  ! Loop index for reaction types '/
     & '        INTEGER IRXN                 ! Reaction number'/
     & '        INTEGER JNUM                 ! J-value species # from PHOT)'/
     & '        INTEGER KNUM                 ! Reaction # for a relative rate coeff.'/
     & '        INTEGER N                    ! Loop index for reactions'/
     & '        INTEGER NCELL                ! Loop index for # of cells in the block' 
     & //
     & '          RKI = 0.0 ' /
     & '          DO NCELL = 1, NUMCELLS ' /
     & '!  Set-up conversion factors '/
     & '             INV_TEMP  = 1.0D+00 / BLKTEMP( NCELL ) '/
     & '             CAIR      = 1.0D+06 * COEF1 * BLKPRES( NCELL ) * INV_TEMP '/
     & '!             CFACT     = 6.0D-05 * CAIR' /
     & '!             CFACT_SQU = 6.0D-11 * CAIR * CAIR '/
     & '!             INV_CFACT = 6.0D+07 / CAIR '/     
     & '             RFACT     = 1.0D-06 * CAIR' / 
     & '             RFACT_SQU = 1.0D-12 * CAIR * CAIR '/
     & '             INV_RFACT = 1.0D+06 / CAIR '/     
     & '             CFACT     = 60.0D0  * RFACT' /
     & '             CFACT_SQU = 60.0D0  * RFACT_SQU '/
     & '             INV_CFACT = 60.0D0  * INV_RFACT '/     
     & '             TEMP      = BLKTEMP( NCELL ) '/
     & '             TEMPOT300 = BLKTEMP( NCELL ) * TI300 ' // )

       RETURN
       
       END SUBROUTINE WRT_RATE_CONSTANT
          
       SUBROUTINE  CONVERT_CASE_LOCAL ( BUFFER, UPPER )
C***********************************************************************

C  subroutine body starts at line  41
C
C  FUNCTION:  converts to upcase or lower the text in BUFFER
C             based on values of logic flag UPPER
C
C  PRECONDITIONS REQUIRED:  text is ASCII
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:  prototype 1/91 by CJC
C
C***********************************************************************

      IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:

        CHARACTER*(*)   BUFFER
        LOGICAL         UPPER


C...........   PARAMETER:  ASCII for 'a', 'z', 'A'

        INTEGER       IA, IZ, AADIF

        PARAMETER   ( IA    = 97,
     &                IZ    = 122,
     &                AADIF = 32 )


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER       I, L
        INTEGER       C
        INTEGER       FACTOR
        INTEGER       STRT, FINI
        


C***********************************************************************
C   begin body of subroutine  UPCASE

        L  =  LEN ( BUFFER )
        IF( UPPER )THEN
            FACTOR =  - AADIF
            STRT   =    IA
            FINI   =    IZ
        ELSE
            FACTOR =    AADIF
            STRT   =    IA - AADIF
            FINI   =    IZ - AADIF
        END IF
        
        DO  111  I = 1 , L
            C = ICHAR ( BUFFER ( I:I ) )
            IF ( C .GE. STRT  .AND.  C .LE. FINI ) THEN
                BUFFER ( I:I ) = CHAR ( C + FACTOR )
            END IF
111     CONTINUE        !  end loop on I

        RETURN
        END SUBROUTINE CONVERT_CASE_LOCAL

      SUBROUTINE WRITE_RATE_CONVERT_LOCAL(OUT_UNIT, RXN_ORDER)
        IMPLICIT NONE
        INTEGER, INTENT( IN ) :: OUT_UNIT
        INTEGER, INTENT( IN ) :: RXN_ORDER
        
         SELECT CASE( RXN_ORDER )
           CASE( 0 )
             WRITE(OUT_UNIT, 95000, ADVANCE = 'NO')
           CASE( 1 )
             WRITE(OUT_UNIT, 95001, ADVANCE = 'NO')
           CASE( 2 )
             WRITE(OUT_UNIT, 95002, ADVANCE = 'NO')
           CASE( 3 )
             WRITE(OUT_UNIT, 95003, ADVANCE = 'NO')
        END SELECT
95000   FORMAT(' INV_CFACT * ')                
95001   FORMAT(' 60.0D0 * ')                
95002   FORMAT(' CFACT * ')                
95003   FORMAT(' CFACT_SQU * ')                
        RETURN
      END SUBROUTINE WRITE_RATE_CONVERT_LOCAL
      SUBROUTINE WRITE_RATE_CONVERT_TIME(OUT_UNIT, RXN_ORDER)
        IMPLICIT NONE
        INTEGER, INTENT( IN ) :: OUT_UNIT
        INTEGER, INTENT( IN ) :: RXN_ORDER
        
         SELECT CASE( RXN_ORDER )
           CASE( 0 )
             WRITE(OUT_UNIT, 95000, ADVANCE = 'NO')
           CASE( 1 )
             WRITE(OUT_UNIT, 95001, ADVANCE = 'NO')
           CASE( 2 )
             WRITE(OUT_UNIT, 95002, ADVANCE = 'NO')
           CASE( 3 )
             WRITE(OUT_UNIT, 95003, ADVANCE = 'NO')
        END SELECT
95000   FORMAT(' INV_RFACT * ')                
95001   FORMAT(' ')                
95002   FORMAT(' RFACT * ')                
95003   FORMAT(' RFACT_SQU * ')                
        RETURN
      END SUBROUTINE WRITE_RATE_CONVERT_TIME