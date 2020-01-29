!#=====================================================
!# Computation of the A and B of the hybrid coordinate.
!# The top pressure is assumed to be equal to zero.
!#=====================================================
!# This version is written in F90,
!# and uses 64 bits variables.
!# 13/06/2017 release in free format.
!#
!# Last modifications since jun 2013:
!#  - compile using mpifort
!#=====================================================
!# Documentation : memoeta.tex (P. Benard, 2004)
!#   Design of the hybrid vertical coordinate "eta".
!#
!# Hybridicity strictly follows (Benard, 2004).
!#
!# Stretching function (ZM) is computed according a
!# new algorithm. Atmosphere is split into JPDOM vertical
!# domains (instead of 3):
!#
!# Tunable variables to compute ZM are:
!# - Height of the lower full level.
!# - Pressure of the upper full level.
!# - Heights and height depths of "vertical domains" interfaces.
!#
!# Altitudes at each half-level are computed using a 5th order polynomial
!# function matching tunable parameters.
!#
!# Number of layers for each "vertical domains" are tunable variables,
!# but the program computes recommended numbers ensuring
!# regular variations of height depths.
!#
!# Standard atmosphere is used to do conversions between
!# standard pressure and standard height; in particuliar,
!# once computed half-level heights, routine PRESTA
!# computes half-level pressures.
!#=====================================================
!# Authors: P. Benard and K. Yessad (MF/CNRM/GMAP)
!#=====================================================

PROGRAM MKVGRID

!     ------------------------------------------------------------------

IMPLICIT NONE

! NAMELIST
INTEGER*8 :: NULNAM=4
CHARACTER*200 :: CDNAM
INTEGER*8 :: IOS

! * INTEGER PARAMETER:
! - JPN: Total number of layers.
INTEGER*8 :: JPN=105
INTEGER*8 :: JPNP
! - JPNPRES: Total number of pure pressure layers (minimum 1 included).
INTEGER*8 :: JPNPRES=23
! - JPNSIGM: Total number of pure sigma layers (minimum 1 included).
INTEGER*8 :: JPNSIGM=1
! - JPDOM: Total number of vertical domains (minimum 3).
INTEGER*8 :: JPDOM=9

! * REAL LOCAL:
REAL*8 :: ZP1
REAL*8 :: ZP1H
REAL*8 :: ZVP00
REAL*8 :: ZVP00PR
REAL*8 :: ZALPH
REAL*8 :: ZX
REAL*8,ALLOCATABLE :: ZALTF(:),ZALTH(:)
REAL*8 :: ZALT_TEST,ZALTH_TEST
REAL*8,ALLOCATABLE :: ZPREF(:),ZPREH(:)
REAL*8,ALLOCATABLE :: ZPREHYDF(:),ZPREHYDH(:),ZDELP(:)
REAL*8 :: ZVP200
REAL*8 :: ZAA
REAL*8 :: ZBB
REAL*8 :: ZETAP
REAL*8 :: ZETAS
REAL*8,ALLOCATABLE :: ZM(:)
REAL*8,ALLOCATABLE :: ZH(:)
REAL*8,ALLOCATABLE :: ZAF(:)
REAL*8,ALLOCATABLE :: ZBF(:)
REAL*8,ALLOCATABLE :: ZRST(:)
REAL*8,ALLOCATABLE :: ZSIG(:)
REAL*8 :: ZALT_PROV(1)
REAL*8 :: ZALTH_PROV(1)
REAL*8 :: ZPRE_PROV(1)
REAL*8 :: ZPREH_PROV(1)
REAL*8,ALLOCATABLE :: ZALT_DOM(:)
REAL*8 :: ZALT_BOT
REAL*8,ALLOCATABLE :: ZDALT_DOM(:)
REAL*8 :: ZDALT_BOT
REAL*8,ALLOCATABLE :: ZI_DOM(:)
REAL*8 :: ZI_TOT
REAL*8 :: ZLEV0,ZLEV1,ZLEV2,ZLEV3,ZLEV4,ZLEV
REAL*8 :: ZVAL0,ZVAL1,ZVAL2,ZVAL3,ZVAL4
REAL*8 :: ZC0,ZC1,ZC2,ZC3,ZC4

! * CHARACTER LOCAL:
CHARACTER*100,ALLOCATABLE :: CLNAM_DOM(:)
CHARACTER*4 :: CPN ! character version of JPN

! * INTEGER LOCAL SCALARS:
INTEGER*8 :: IFILE0,IFILE1,IFILE2,IFILE3,ILEV,JK,JLEV,I1
INTEGER*8 :: IN,INP,INPRES,INSIGM,II,IDOM,JDOM
INTEGER*8 :: I_BOT,I_TOT
INTEGER*8,ALLOCATABLE :: I_DOM(:)
INTEGER*8,ALLOCATABLE :: IT_DOM(:)

! * LOGICAL LOCAL SCALARS:
LOGICAL :: LLAPRXPK

! * ADDITIONAL DECLARATIONS FOR "STANDARD Z" CALCULATIONS:
REAL*8,ALLOCATABLE :: ZSTZ(:),ZSTZH(:)

!     =======================================================================
!     * Get namelist

NAMELIST/NAM_DIM/JPN,JPNPRES,JPNSIGM,JPDOM
NAMELIST/NAM_REF/ZP1,ZVP00,ZVP00PR,ZVP200,ZALT_BOT
NAMELIST/NAM_DOM/ZALT_DOM,ZDALT_DOM,IT_DOM,CLNAM_DOM
NAMELIST/NAM_PARAM/LLAPRXPK,ZALPH

IF (IARGC() /= 1) THEN
  WRITE(6,*) 'Error: must provide one and only one argument: path to namelist.'
  STOP 1
ENDIF
CALL GETARG(1, CDNAM)

!     =======================================================================
!     * Defaults

ZP1=10.
! * ZVP00: standard surface pressure (to compute A and B), in pascals:
ZVP00=101325.
! * ZVP00PR: standard surface pressure (to print pressures), in pascals:
ZVP00PR=101325.
! * ZVP200: standard pressure at 200m (mean orography on the Earth) in Pa.:
ZVP200=98945.37974
! * ZALT_BOT: altitude of the bottom full level, in metres.
ZALT_BOT=10.
! * LLAPRXPK:
!   Full layers are assumed to be computed as for the options
!   LVERTFE=F, NDLNPR=0 of ARPEGE/ALADIN.
!   LLAPRXPK=T => pressure(l)=0.5(pressure(lbar-1)+pressure(lbar))
!    ("l" stands for full levels, "lbar" for half levels).
!   LLAPRXPK=F => a more tricky way to compute pressure(l).
!   When using the vertical layers for LVERTFE=F, NDLNPR=0, LAPRXPK=F
!    in the model, it is recommended to use LLAPRXPK=F.
!   When using the vertical layers for LVERTFE=F, NDLNPR=1,
!    in the model, it is recommended to use LLAPRXPK=F.
!   When using the vertical layers for LVERTFE=F, NDLNPR=0, LAPRXPK=T
!    of for LVERTFE=T, it is recommended to use LLAPRXPK=T.
LLAPRXPK=.TRUE.
! * Controls the hybridicity (recommended values between -3 and -1,
!   ZALPH must never be > -1).
ZALPH=-1.6

!     =======================================================================
!     * Read namelist and initializations, allocations


OPEN(NULNAM,FILE=CDNAM,ACTION='READ',IOSTAT=IOS)
IF (IOS /= 0) THEN
  WRITE(6,*) 'Error reading namelist'
  STOP 1
ENDIF
! dimensions
CALL POSNAME(NULNAM, 'NAM_DIM', IOS)
READ(NULNAM, NAM_DIM)
JPNP = JPN + 1
WRITE(CPN,'(I0)') JPN
! for allocations
ALLOCATE(ZALTF(JPN))
ALLOCATE(ZALTH(JPNP))
ALLOCATE(ZPREF(JPN))
ALLOCATE(ZPREH(JPNP))
ALLOCATE(ZPREHYDF(JPN))
ALLOCATE(ZPREHYDH(JPNP))
ALLOCATE(ZDELP(JPN))
ALLOCATE(ZM(0:JPN))
ALLOCATE(ZH(0:JPN))
ALLOCATE(ZAF(JPNP))
ALLOCATE(ZBF(JPNP))
ALLOCATE(ZRST(JPNP))
ALLOCATE(ZSIG(JPNP))
ALLOCATE(ZALT_DOM(JPDOM))
ALLOCATE(ZDALT_DOM(JPDOM))
ALLOCATE(ZI_DOM(JPDOM))
ALLOCATE(CLNAM_DOM(JPDOM))
ALLOCATE(I_DOM(JPDOM))
ALLOCATE(IT_DOM(JPDOM))
ALLOCATE(ZSTZ(JPN))
ALLOCATE(ZSTZH(JPNP))

! Initialize DOMains to ensure namelist consistency
ZALT_DOM(:) = -HUGE(8)
ZDALT_DOM(:) = -HUGE(8)
IT_DOM(:) = -HUGE(8)
CLNAM_DOM(:) = ''

! Then read namelist
CALL POSNAME(NULNAM, 'NAM_REF', IOS)
READ(NULNAM, NAM_REF)
CALL POSNAME(NULNAM, 'NAM_DOM', IOS)
READ(NULNAM, NAM_DOM)
CALL POSNAME(NULNAM, 'NAM_PARAM', IOS)
READ(NULNAM, NAM_PARAM)
CLOSE(NULNAM)

! consistency checks
IF (SUM(IT_DOM(1:JPDOM)) /= JPN) THEN
  WRITE(6,*) 'Error: Inconsistency JPN /= SUM(IT_DOM)'
  STOP 1
ELSEIF (SUM(ZALT_DOM(1:JPDOM-1)) < 0.) THEN
  WRITE(6,*) 'Error: ZALT_DOM must be specified from 1 to JPDOM-1'
  STOP 1
ELSEIF (SUM(ZDALT_DOM(1:JPDOM)) < 0.) THEN
  WRITE(6,*) 'Error: ZDALT_DOM must be specified from 1 to JPDOM'
  STOP 1
ELSEIF (ANY(CLNAM_DOM == '')) THEN
  WRITE(6,*) 'Error: CLNAM_DOM must be specified from 1 to JPDOM'
  STOP 1
ENDIF

!     =======================================================================
! open output files
IFILE0 = 50
IFILE1 = 51
IFILE2 = 52
IFILE3 = 53

OPEN(IFILE0,FILE=TRIM(CDNAM)//'.mapping',ACTION='WRITE',IOSTAT=IOS)
IF (IOS /= 0) THEN
  WRITE(6,*) 'Error opening mapping file'
  STOP 1
ENDIF
OPEN(IFILE1,FILE=TRIM(CDNAM)//'.infosup',ACTION='WRITE',IOSTAT=IOS)
IF (IOS /= 0) THEN
  WRITE(6,*) 'Error opening infosup file'
  STOP 1
ENDIF
OPEN(IFILE2,FILE=TRIM(CDNAM)//'.latex',ACTION='WRITE',IOSTAT=IOS)
IF (IOS /= 0) THEN
  WRITE(6,*) 'Error opening latex file'
  STOP 1
ENDIF
OPEN(IFILE3,FILE=TRIM(CDNAM)//'.namvv1',ACTION='WRITE',IOSTAT=IOS)
IF (IOS /= 0) THEN
  WRITE(6,*) 'Error opening namvv1 file'
  STOP 1
ENDIF
!     =======================================================================

IN=JPN
INP=JPNP
INPRES=JPNPRES
INSIGM=JPNSIGM 
IDOM=JPDOM

I1=1

!     =======================================================================

!     * Computation of mapping function ZM(JLEV)=m(x(JLEV))
!       ZM is equal to pressure/ZVP00 at locations where surface pressure is ZVP00
!       ZM is 0 at the top and 1 at the bottom.

! - Treatment of bottom layer.
ZALT_PROV(1)=ZALT_BOT
CALL PRESTA(I1,ZALT_PROV,ZPRE_PROV)
! exact for LLAPRXPK=T, close to right solution for LLAPRXPK=F
ZPREH_PROV(1)=ZVP00+2.*(ZPRE_PROV(1)-ZVP00)
CALL ALTSTA(I1,ZPREH_PROV,ZALTH_PROV)
ZALTF(IN)=ZALT_BOT
ZALTH(INP)=0.
ZALTH(IN)=ZALTH_PROV(1)
ZPREF(IN)=ZPRE_PROV(1)
ZPREH(INP)=ZVP00
ZPREH(IN)=ZPREH_PROV(1)
ZDALT_BOT=ZALTH(IN)-ZALTH(INP)

! - Treatment of top layer.
IF (LLAPRXPK) THEN
  ZP1H=2.*ZP1
ELSE
  ZP1H=EXP(1.)*ZP1
ENDIF
ZPRE_PROV(1)=ZP1
ZPREH_PROV(1)=ZP1H
CALL ALTSTA(I1,ZPRE_PROV,ZALT_PROV)
CALL ALTSTA(I1,ZPREH_PROV,ZALTH_PROV)
ZALTF(1)=ZALT_PROV(1)
ZALTH(2)=ZALTH_PROV(1)
ZPREF(1)=ZPRE_PROV(1)
ZPREH(2)=ZPREH_PROV(1)
ZALT_DOM(IDOM)=ZALTH(2)

! - Evaluate number of recommended levels in each "vertical domain",
!   in order to have regular variations of altitude depths.

I_BOT=1
ZI_DOM(1)=2.*(ZALT_DOM(1)-ZALTH(IN))/(ZDALT_BOT+ZDALT_DOM(1))
I_DOM(1)=NINT(ZI_DOM(1))
DO JDOM=2,IDOM-1
  ZI_DOM(JDOM)=2.*(ZALT_DOM(JDOM)-ZALT_DOM(JDOM-1))/(ZDALT_DOM(JDOM)+ZDALT_DOM(JDOM-1))
  I_DOM(JDOM)=NINT(ZI_DOM(JDOM))
ENDDO
ZI_DOM(IDOM)=1.+2.*(ZALT_DOM(IDOM)-ZALT_DOM(IDOM-1))/(ZDALT_DOM(IDOM)+ZDALT_DOM(IDOM-1))
I_DOM(IDOM)=NINT(ZI_DOM(IDOM))
ZI_TOT=1.+SUM(ZI_DOM(1:IDOM))
I_TOT=I_BOT+SUM(I_DOM(1:IDOM))

write(6,*) ' For requested levels and levels depths, recommended numbers of layers are:'
write(6,'(A,I2,A,I2)') ' recommended value of IT_DOM( 1) around ',I_DOM(1)+I_BOT, ' // namelist=',IT_DOM(1)
DO JDOM=2,IDOM
  write(6,'(A,I2,A,I2,A,I2)') ' recommended value of IT_DOM(',JDOM,') around ',I_DOM(JDOM), ' // namelist=',IT_DOM(JDOM)
ENDDO
write(6,'(A,I3)') ' recommended value of JPN around ',NINT(ZI_TOT)

DO JDOM=1,IDOM
  IF (I_DOM(JDOM) < 2) THEN
    write(6,'(A,A,A7)') ' Inconsistent reference values for vertical domain ',TRIM(CLNAM_DOM(JDOM))
    write(6,'(A,I2,A,I2,A)') ' Raise ZALT_DOM(',JDOM,') or reduce ZDALT_DOM(',JDOM,')'
    write(6,*) ' Program stops.'
    STOP
  ENDIF
ENDDO

! - Treatment of vertical domain nr 1.
 
ILEV=IN
ZLEV1=REAL(ILEV+1,8)
ZLEV2=REAL(ILEV,8)
ZLEV3=REAL(ILEV-IT_DOM(1)+2,8)
ZLEV4=REAL(ILEV-IT_DOM(1)+1,8)
ZVAL1=ZALTH(ILEV+1)
ZVAL2=ZALTH(ILEV)
ZVAL3=ZALT_DOM(1)-ZDALT_DOM(1)
ZVAL4=ZALT_DOM(1)
DO JLEV=1,IT_DOM(1)-3
  ILEV=ILEV-1
  ZLEV=REAL(ILEV,8)
  ZC1=(ZLEV-ZLEV2)*(ZLEV-ZLEV3)*(ZLEV-ZLEV4)/((ZLEV1-ZLEV2)*(ZLEV1-ZLEV3)*(ZLEV1-ZLEV4))
  ZC2=(ZLEV-ZLEV1)*(ZLEV-ZLEV3)*(ZLEV-ZLEV4)/((ZLEV2-ZLEV1)*(ZLEV2-ZLEV3)*(ZLEV2-ZLEV4))
  ZC3=(ZLEV-ZLEV1)*(ZLEV-ZLEV2)*(ZLEV-ZLEV4)/((ZLEV3-ZLEV1)*(ZLEV3-ZLEV2)*(ZLEV3-ZLEV4))
  ZC4=(ZLEV-ZLEV1)*(ZLEV-ZLEV2)*(ZLEV-ZLEV3)/((ZLEV4-ZLEV1)*(ZLEV4-ZLEV2)*(ZLEV4-ZLEV3))
  ZALTH(ILEV)=ZC1*ZVAL1+ZC2*ZVAL2+ZC3*ZVAL3+ZC4*ZVAL4
ENDDO
ILEV=ILEV-1
ZALTH(ILEV)=ZVAL3
ILEV=ILEV-1
ZALTH(ILEV)=ZVAL4

! - Treatment of vertical domains nr 2 to IDOM-1.

DO JDOM=2,IDOM-1
  ZLEV0=REAL(ILEV+2,8)
  ZLEV1=REAL(ILEV+1,8)
  ZLEV2=REAL(ILEV,8)
  ZLEV3=REAL(ILEV-IT_DOM(JDOM)+1,8)
  ZLEV4=REAL(ILEV-IT_DOM(JDOM),8)
  ZVAL0=ZALTH(ILEV+2)
  ZVAL1=ZALTH(ILEV+1)
  ZVAL2=ZALTH(ILEV)
  ZVAL3=ZALT_DOM(JDOM)-ZDALT_DOM(JDOM)
  ZVAL4=ZALT_DOM(JDOM)
  DO JLEV=1,IT_DOM(JDOM)-2
    ILEV=ILEV-1
    ZLEV=REAL(ILEV,8)
    ZC0=(ZLEV-ZLEV1)*(ZLEV-ZLEV2)*(ZLEV-ZLEV3)*(ZLEV-ZLEV4)/((ZLEV0-ZLEV1)*(ZLEV0-ZLEV2)*(ZLEV0-ZLEV3)*(ZLEV0-ZLEV4))
    ZC1=(ZLEV-ZLEV0)*(ZLEV-ZLEV2)*(ZLEV-ZLEV3)*(ZLEV-ZLEV4)/((ZLEV1-ZLEV0)*(ZLEV1-ZLEV2)*(ZLEV1-ZLEV3)*(ZLEV1-ZLEV4))
    ZC2=(ZLEV-ZLEV0)*(ZLEV-ZLEV1)*(ZLEV-ZLEV3)*(ZLEV-ZLEV4)/((ZLEV2-ZLEV0)*(ZLEV2-ZLEV1)*(ZLEV2-ZLEV3)*(ZLEV2-ZLEV4))
    ZC3=(ZLEV-ZLEV0)*(ZLEV-ZLEV1)*(ZLEV-ZLEV2)*(ZLEV-ZLEV4)/((ZLEV3-ZLEV0)*(ZLEV3-ZLEV1)*(ZLEV3-ZLEV2)*(ZLEV3-ZLEV4))
    ZC4=(ZLEV-ZLEV0)*(ZLEV-ZLEV1)*(ZLEV-ZLEV2)*(ZLEV-ZLEV3)/((ZLEV4-ZLEV0)*(ZLEV4-ZLEV1)*(ZLEV4-ZLEV2)*(ZLEV4-ZLEV3))
    ZALTH(ILEV)=ZC0*ZVAL0+ZC1*ZVAL1+ZC2*ZVAL2+ZC3*ZVAL3+ZC4*ZVAL4
  ENDDO
  ILEV=ILEV-1
  ZALTH(ILEV)=ZVAL3
  ILEV=ILEV-1
  ZALTH(ILEV)=ZVAL4
ENDDO

! - Treatment of vertical domain nr IDOM.

II=IT_DOM(IDOM)-1
ZLEV0=REAL(ILEV+2,8)
ZLEV1=REAL(ILEV+1,8)
ZLEV2=REAL(ILEV,8)
ZLEV3=REAL(ILEV-II+1,8)
ZLEV4=REAL(ILEV-II,8)
ZVAL0=ZALTH(ILEV+2)
ZVAL1=ZALTH(ILEV+1)
ZVAL2=ZALTH(ILEV)
ZVAL3=ZALT_DOM(IDOM)-ZDALT_DOM(IDOM)
ZVAL4=ZALT_DOM(IDOM)
DO JLEV=1,II-2
  ILEV=ILEV-1
  ZLEV=REAL(ILEV,8)
  ZC0=(ZLEV-ZLEV1)*(ZLEV-ZLEV2)*(ZLEV-ZLEV3)*(ZLEV-ZLEV4)/((ZLEV0-ZLEV1)*(ZLEV0-ZLEV2)*(ZLEV0-ZLEV3)*(ZLEV0-ZLEV4))
  ZC1=(ZLEV-ZLEV0)*(ZLEV-ZLEV2)*(ZLEV-ZLEV3)*(ZLEV-ZLEV4)/((ZLEV1-ZLEV0)*(ZLEV1-ZLEV2)*(ZLEV1-ZLEV3)*(ZLEV1-ZLEV4))
  ZC2=(ZLEV-ZLEV0)*(ZLEV-ZLEV1)*(ZLEV-ZLEV3)*(ZLEV-ZLEV4)/((ZLEV2-ZLEV0)*(ZLEV2-ZLEV1)*(ZLEV2-ZLEV3)*(ZLEV2-ZLEV4))
  ZC3=(ZLEV-ZLEV0)*(ZLEV-ZLEV1)*(ZLEV-ZLEV2)*(ZLEV-ZLEV4)/((ZLEV3-ZLEV0)*(ZLEV3-ZLEV1)*(ZLEV3-ZLEV2)*(ZLEV3-ZLEV4))
  ZC4=(ZLEV-ZLEV0)*(ZLEV-ZLEV1)*(ZLEV-ZLEV2)*(ZLEV-ZLEV3)/((ZLEV4-ZLEV0)*(ZLEV4-ZLEV1)*(ZLEV4-ZLEV2)*(ZLEV4-ZLEV3))
  ZALTH(ILEV)=ZC0*ZVAL0+ZC1*ZVAL1+ZC2*ZVAL2+ZC3*ZVAL3+ZC4*ZVAL4
ENDDO
ILEV=ILEV-1
ZALTH(ILEV)=ZVAL3
ILEV=ILEV-2

IF (ILEV /= 1) THEN
  write(6,*) ' Error in ILEV in upper zone'
  write(6,*) ' ILEV must be equal to 1'
  write(6,*) ' ILEV=',ILEV
  write(6,*) ' Program stops.'
  STOP
ENDIF

! - Compute half-level pressures from half-level altitudes using standard atm. 
CALL PRESTA(IN-1,ZALTH(2:IN),ZPREH(2:IN))

! - Fill ZM. 

DO JLEV=0,IN
  ZM(JLEV)=ZPREH(JLEV+1)/ZVP00
  ZM(JLEV)=MAX(0.,MIN(1.,ZM(JLEV)))
ENDDO

WRITE(IFILE0,*) " ZM: "
DO JLEV=0,IN
  WRITE(IFILE0,'(1X,F15.7)') ZM(JLEV)
ENDDO

!     =======================================================================

!     * Computation of *mapped hybridicity* function ZH(JLEV)=h(m(x(JLEV)))

ZETAP=ZM(INPRES)
ZETAS=ZM(IN-INSIGM)

ZAA=ZALPH*ZETAS*ZETAS/(ZETAS-ZETAP)
ZBB=1.+ZAA/ZETAS

WRITE(IFILE0,*) " ZETAP: ",ZETAP
WRITE(IFILE0,*) " ZETAS: ",ZETAS
WRITE(IFILE0,*) " ZAA: ",ZAA
WRITE(IFILE0,*) " ZBB: ",ZBB

DO JLEV=0,INPRES
  ZH(JLEV)=0.
ENDDO

DO JLEV=INPRES+1,(IN-INSIGM)-1
  ZX=ZM(JLEV)
  ZH(JLEV)=ZAA/(ZBB-((ZX-ZETAP)/(ZETAS-ZETAP))**ZALPH)
ENDDO

DO JLEV=(IN-INSIGM),IN
  ZX=ZM(JLEV)
  ZH(JLEV)=ZX
ENDDO

DO JLEV=0,IN
  ZH(JLEV)=MAX(0.,MIN(1.,ZH(JLEV)))
ENDDO

WRITE(IFILE0,*) " ZH: "
DO JLEV=0,IN
  WRITE(IFILE0,'(1X,F15.7)') ZH(JLEV)
ENDDO

!     =======================================================================

!     * A and B functions on half layers (put in ZAF and ZBF):

DO JLEV=0,IN
  ZAF(JLEV+1)=ZVP00*(ZM(JLEV)-ZH(JLEV))
  ZBF(JLEV+1)=ZH(JLEV)
ENDDO

!     =======================================================================

!     * Half level pressures (put in ZPREHYDH):

DO JK=1,INP
  ZPREHYDH(JK)=ZAF(JK)+ZBF(JK)*ZVP00PR
ENDDO

!     =======================================================================

!     * Full level pressures (put in ZPREHYDF):

IF (.NOT.LLAPRXPK) THEN
  ZPREHYDF(1)=ZPREHYDH(2)/EXP(1.)
  DO JLEV=2,IN
    ZPREHYDF(JLEV)=EXP( &
     &    (ZPREHYDH(JLEV+1)*LOG(ZPREHYDH(JLEV+1)) &
     &     -ZPREHYDH(JLEV)*LOG(ZPREHYDH(JLEV)))/ &
     &    (ZPREHYDH(JLEV+1)-ZPREHYDH(JLEV)) - 1. &
     &    )
  ENDDO
ELSE
  ZPREHYDF(1)=0.5*ZPREHYDH(2)
  DO JLEV=2,IN
    ZPREHYDF(JLEV)=0.5*(ZPREHYDH(JLEV+1)+ZPREHYDH(JLEV))
  ENDDO
ENDIF

!     =======================================================================

!     * Pressure depths at full levels (put in ZDELP):

DO JLEV=1,IN
  ZDELP(JLEV)=ZPREHYDH(JLEV+1)-ZPREHYDH(JLEV)
ENDDO

!     =======================================================================

!     * Standard altitude at full levels (put in ZSTZ):

CALL ALTSTA(IN,ZPREHYDF,ZSTZ)

!     =======================================================================

!     * Standard altitude at half levels (put in ZSTZH):

CALL ALTSTA(IN,ZPREHYDH(2:INP),ZSTZH(2:INP))
ZSTZH(1)=ZSTZH(2)+2.*(ZSTZ(1)-ZSTZH(2))

!     =======================================================================

!     * Printings:

DO JK=1,INP
  ZRST(JK)=ZBF(JK)*ZVP200/(ZBF(JK)*ZVP200+ZAF(JK))
  ZSIG(JK)=ZAF(JK)/ZVP200+ZBF(JK)
ENDDO

WRITE(IFILE1,*)
WRITE(IFILE1,'(1X,A,1X,I4)') ' * Number of levels=',IN
WRITE(IFILE1,'(1X,A,1X,I4)') ' * Number of hybrid levels=',IN-INPRES-INSIGM
WRITE(IFILE1,'(1X,A,1X,I4)') ' * Number of pure sigma levels=',INSIGM
WRITE(IFILE1,'(1X,A,1X,I4)') ' * Number of pure pressure levels=',INPRES
WRITE(IFILE1,'(1X,A,1X,I4)') ' * Number of vertical domains=',IDOM

WRITE(IFILE1,*)
DO JDOM=1,IDOM
  WRITE(IFILE1,'(1X,A,A,A,1X,I4)') &
     &   ' * Number of levels in vertical domain ',TRIM(CLNAM_DOM(JDOM)),' =',IT_DOM(JDOM)
ENDDO

WRITE(IFILE1,*)
WRITE(IFILE1,'(1X,A,1X,F15.7,1X,A,1X,F15.7)') &
 & ' * Altitude of bottom full level:              ZALT_BOT    =',ZALT_BOT, &
 & ';  Delta z: ZDALT_BOT    =',ZDALT_BOT
DO JDOM=1,IDOM-1
  WRITE(IFILE1,'(1X,3A,I2,A,1X,F15.7,1X,A,I2,A,1X,F15.7)') &
   &   ' * Altitude of top of vertical domain ', &
   &  TRIM(CLNAM_DOM(JDOM)),': ZALT_DOM(',JDOM,')=',ZALT_DOM(JDOM), &
   &   ';  Delta z: ZDALT_DOM(',JDOM,')=',ZDALT_DOM(JDOM)
ENDDO
WRITE(IFILE1,'(1X,A,29X,F15.7,1X,A,I2,A,1X,F15.7)') &
 & ' * Altitude of half level nr 1: ',ZSTZH(2), &
 & ';  Delta z: ZDALT_DOM(',IDOM,')=',ZDALT_DOM(IDOM)

WRITE(IFILE1,*)
WRITE(IFILE1,'(1X,A,1X,F15.7)') ' * Reference pressure at 200m: ZVP200=',ZVP200
WRITE(IFILE1,'(1X,A,1X,F15.7)') ' * Reference pressure at 0m: ZVP00PR=',ZVP00PR
WRITE(IFILE1,'(1X,A,1X,F15.7)') ' * Pressure of the layer nr 1: ZP1=',ZP1

WRITE(IFILE1,*)
WRITE(IFILE1,'(1X,A,1X,L2)') ' * LLAPRXPK=',LLAPRXPK
WRITE(IFILE1,'(1X,A,1X,F15.7)') ' * ZALPH=',ZALPH

WRITE(IFILE1,*)
WRITE(IFILE1,'(A3,A13,A13,3A10,3A14,3A12)') &
 & 'ILa','      A      ','      B      ', &
 & '   Sigma  ',' 1 - Sigma',' Rap Si-Hy', &
 & '  Prehyd(lbar)','  Prehyd(l)   ', &
 & ' [D Prehyd](l)','   Alti(l)  ',' Alti(lbar) ', &
 & ' [D Alti](l)'
WRITE(IFILE1,*)

DO JK=1,1
  WRITE(IFILE1,'(I3,F13.6,F13.10,3F10.7,F12.4)') &
   & JK-1,ZAF(JK),ZBF(JK),ZSIG(JK),1.-ZSIG(JK),ZRST(JK),ZPREHYDH(JK)
ENDDO
DO JK=2,INP
  WRITE(IFILE1,'(I3,F13.6,F13.10,3F10.7,3F14.6,3F12.4)') &
   & JK-1,ZAF(JK),ZBF(JK),ZSIG(JK),1.-ZSIG(JK),ZRST(JK), &
   & ZPREHYDH(JK),ZPREHYDF(JK-1),ZDELP(JK-1), &
   & ZSTZ(JK-1),ZSTZH(JK),ZSTZH(JK-1)-ZSTZH(JK)
ENDDO

WRITE(IFILE2,*)
DO JK=1,INP
  WRITE(IFILE2,'(A12,F7.3,A28,I3)') &
   & '\\put( 20.00,',100.*(1.-ZSIG(JK)), &
   & '){\\line(1,0){50.00}} % jlev=',JK-1
ENDDO

WRITE(IFILE3,'(A)') ' Namelist obtained with:'

WRITE(IFILE3,*)
WRITE(IFILE3,'(1X,A,1X,I4)') ' * Number of levels=',IN
WRITE(IFILE3,'(1X,A,1X,I4)') ' * Number of hybrid levels=',IN-INPRES-INSIGM
WRITE(IFILE3,'(1X,A,1X,I4)') ' * Number of pure sigma levels=',INSIGM
WRITE(IFILE3,'(1X,A,1X,I4)') ' * Number of pure pressure levels=',INPRES
WRITE(IFILE3,'(1X,A,1X,I4)') ' * Number of vertical domains=',IDOM

WRITE(IFILE3,*)
DO JDOM=1,IDOM
  WRITE(IFILE3,'(1X,A,A,A,1X,I4)') &
   & ' * Number of levels in vertical domain ', &
   & TRIM(CLNAM_DOM(JDOM)),' =',IT_DOM(JDOM)
ENDDO

WRITE(IFILE3,*)
WRITE(IFILE3,'(1X,A,1X,F15.7,1X,A,1X,F15.7)') &
 & ' * Altitude of bottom full level:              ZALT_BOT    =',ZALT_BOT, &
 & ';  Delta z: ZDALT_BOT    =',ZDALT_BOT
DO JDOM=1,IDOM-1
  WRITE(IFILE3,'(1X,3A,I2,A,1X,F15.7,1X,A,I2,A,1X,F15.7)') &
   & ' * Altitude of top of vertical domain ', &
   & TRIM(CLNAM_DOM(JDOM)),': ZALT_DOM(',JDOM,')=',ZALT_DOM(JDOM), &
   & ';  Delta z: ZDALT_DOM(',JDOM,')=',ZDALT_DOM(JDOM)
ENDDO
WRITE(IFILE3,'(1X,A,29X,F15.7,1X,A,I2,A,1X,F15.7)') &
 & ' * Altitude of half level nr 1: ',ZSTZH(2), &
 & ';  Delta z: ZDALT_DOM(',IDOM,')=',ZDALT_DOM(IDOM)

WRITE(IFILE3,*)
WRITE(IFILE3,'(1X,A,1X,F15.7)') ' * Reference pressure at 200m: ZVP200=',ZVP200
WRITE(IFILE3,'(1X,A,1X,F15.7)') ' * Reference pressure at 0m: ZVP00PR=',ZVP00PR
WRITE(IFILE3,'(1X,A,1X,F15.7)') ' * Pressure of the layer nr 1: ZP1=',ZP1

WRITE(IFILE3,*)
WRITE(IFILE3,'(1X,A,1X,L2)') ' * LLAPRXPK=',LLAPRXPK
WRITE(IFILE3,'(1X,A,1X,F15.7)') ' * ZALPH=',ZALPH

WRITE(IFILE3,*)
WRITE(IFILE3,*) ' &NAMVV1'
WRITE(IFILE3,*) '   DVALH(  0)=0.,'
DO JK=2,INP
  ! WRITE(IFILE3,'(5X,F12.6,A)') ZAF(JK),','
  WRITE(IFILE3,'(A,I3,A,F10.4,A)') '    DVALH(',JK-1,')=',ZAF(JK),','
ENDDO
WRITE(IFILE3,*) '   DVBH(  0)=0.,'
DO JK=2,INP
  ! WRITE(IFILE3,'(5X,F12.10,A)') ZBF(JK),','
  WRITE(IFILE3,'(A,I3,A,F10.8,A)') '    DVBH(',JK-1,')=',ZBF(JK),','
ENDDO
WRITE(IFILE3,*) ' /'

!     ------------------------------------------------------------------
DEALLOCATE(ZALTF)
DEALLOCATE(ZALTH)
DEALLOCATE(ZPREF)
DEALLOCATE(ZPREH)
DEALLOCATE(ZPREHYDF)
DEALLOCATE(ZPREHYDH)
DEALLOCATE(ZDELP)
DEALLOCATE(ZM)
DEALLOCATE(ZH)
DEALLOCATE(ZAF)
DEALLOCATE(ZBF)
DEALLOCATE(ZRST)
DEALLOCATE(ZSIG)
DEALLOCATE(ZALT_DOM)
DEALLOCATE(ZDALT_DOM)
DEALLOCATE(ZI_DOM)
DEALLOCATE(CLNAM_DOM)
DEALLOCATE(I_DOM)
DEALLOCATE(IT_DOM)
DEALLOCATE(ZSTZ)

CLOSE(IFILE0)
CLOSE(IFILE1)
CLOSE(IFILE2)
CLOSE(IFILE3)

STOP
END PROGRAM

