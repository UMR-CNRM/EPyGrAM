SUBROUTINE ALTSTA(KLEV,PREHYD,PSTZ)

!     ------------------------------------------------------------------
!     Computes standard altitude from standard pressure.
!     Code taken from SUSTA and PPSTA and adapted.

!     Input:
!      KLEV    : number of levels
!      PREHYD  : standard pressure

!     Output:
!      PSTZ    : standard altitude
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER*8,INTENT(IN) :: KLEV
REAL*8,INTENT(IN)    :: PREHYD(KLEV)
REAL*8,INTENT(OUT)   :: PSTZ(KLEV)

!     ------------------------------------------------------------------

INTEGER*8,PARAMETER :: JPPRO=9

REAL*8 :: ZVP00
LOGICAL :: LLGRZER_ST
INTEGER*8 :: JLEV,IPR,JJPR,IPRO_ST(KLEV)
REAL*8 :: ZSUR,Z_RD,Z_RG,ZROG
REAL*8 :: Z_RZTROP
REAL*8 :: Z_RZSTRA
REAL*8 :: Z_RZSTR2
REAL*8 :: Z_RZSTPO
REAL*8 :: Z_RZMESO
REAL*8 :: Z_RZMES2
REAL*8 :: Z_RZMEPO
REAL*8 :: Z_RZABOV
REAL*8 :: Z_RDTDZ1
REAL*8 :: Z_RDTDZ2
REAL*8 :: Z_RDTDZ3
REAL*8 :: Z_RDTDZ4
REAL*8 :: Z_RDTDZ5
REAL*8 :: Z_RDTDZ6
REAL*8 :: Z_RDTDZ7
REAL*8 :: Z_RDTDZ8
REAL*8 :: Z_RDTDZ9
REAL*8 :: Z_RTSUR
REAL*8 :: ZZDT_ST,ZZDZ_ST
REAL*8 :: ZZ_ST(JPPRO)
REAL*8 :: ZDTDZ_ST(JPPRO)
REAL*8 :: ZT_ST(JPPRO)
REAL*8 :: ZP_ST(JPPRO)
REAL*8 :: ZLNP_ST(JPPRO)
REAL*8 :: ZIP_ST(JPPRO)
REAL*8 :: ZIDTDZ_ST(JPPRO)

!     =======================================================================

!     * Standard altitude (put in PSTZ):
!       These calculations have been adapted from SUSTA and PPSTA

ZSUR =0.

! - constants:
! Z_RD=1000.*6.0221367*1.380658/28.9644
Z_RD=287.0596737
Z_RG=9.80665

! - standard altitude at some referenced levels:
Z_RZTROP=11000.
Z_RZSTRA=20000.
Z_RZSTR2=32100.
Z_RZSTPO=47400.
Z_RZMESO=51400.
Z_RZMES2=71700.
Z_RZMEPO=85700.
Z_RZABOV=100000.

! - standard DT/Dz at some referenced atmosphere layers:
Z_RDTDZ1=-6.5E-3
Z_RDTDZ2=0.
Z_RDTDZ3=1.E-3
Z_RDTDZ4=2.75E-3
Z_RDTDZ5=0.
Z_RDTDZ6=-2.75E-3
Z_RDTDZ7=-1.97E-3
Z_RDTDZ8=0.
Z_RDTDZ9=0.

! - surface standard temperature:
Z_RTSUR=288.15

! - surface standard pressure:
ZVP00=101325.

DO JJPR=1,JPPRO
  IF(JJPR == 1) THEN
    ZZ_ST(1)=ZSUR
    ZDTDZ_ST(1)=Z_RDTDZ1
    ZT_ST(1)=Z_RTSUR
    ZP_ST(1)=ZVP00
  ELSEIF(JJPR == 2) THEN
    ZZ_ST(2)=Z_RZTROP 
    ZDTDZ_ST(2)=Z_RDTDZ2
  ELSEIF(JJPR == 3) THEN
    ZZ_ST(3)=Z_RZSTRA 
    ZDTDZ_ST(3)=Z_RDTDZ3
  ELSEIF(JJPR == 4) THEN
    ZZ_ST(4)=Z_RZSTR2 
    ZDTDZ_ST(4)=Z_RDTDZ4
  ELSEIF(JJPR == 5) THEN
    ZZ_ST(5)=Z_RZSTPO 
    ZDTDZ_ST(5)=Z_RDTDZ5
  ELSEIF(JJPR == 6) THEN
    ZZ_ST(6)=Z_RZMESO
    ZDTDZ_ST(6)=Z_RDTDZ6
  ELSEIF(JJPR == 7) THEN
    ZZ_ST(7)=Z_RZMES2
    ZDTDZ_ST(7)=Z_RDTDZ7
  ELSEIF(JJPR == 8) THEN
    ZZ_ST(8)=Z_RZMEPO
    ZDTDZ_ST(8)=Z_RDTDZ8
  ELSEIF(JJPR == 9) THEN
    ZZ_ST(9)=Z_RZABOV
    ZDTDZ_ST(9)=Z_RDTDZ9
  ENDIF
ENDDO

! - standard temperature and pressure at some referenced levels:
DO JJPR=1,JPPRO-1
  ZZDT_ST=ZDTDZ_ST(JJPR)
  ZZDZ_ST=ZZ_ST(JJPR+1)-ZZ_ST(JJPR)
  IF (ZZDT_ST /= 0.0) THEN
    ZT_ST(JJPR+1)=ZT_ST(JJPR)+ZZDT_ST*ZZDZ_ST
    ZP_ST(JJPR+1)=ZP_ST(JJPR)*(1.0+ZZDT_ST*ZZDZ_ST/ZT_ST(JJPR))**(-Z_RG/Z_RD/ZZDT_ST)
  ELSE
    ZT_ST(JJPR+1)=ZT_ST(JJPR)
    ZP_ST(JJPR+1)=ZP_ST(JJPR)*EXP(-Z_RG/Z_RD/ZT_ST(JJPR)*ZZDZ_ST)
  ENDIF
ENDDO

! LLGRZER_ST : .true. if one or more values of ZDTDZ are zero, .false. otherwise
LLGRZER_ST=.FALSE.
DO JJPR=1,JPPRO
  ZLNP_ST(JJPR)=LOG(ZP_ST(JJPR))
  ZIP_ST(JJPR)=1.0/ZP_ST(JJPR)
  IF (ZDTDZ_ST(JJPR) == 0.0) THEN
    ZIDTDZ_ST(JJPR)=0.0
    LLGRZER_ST=.TRUE.
  ELSE
    ZIDTDZ_ST(JJPR)=1.0/ZDTDZ_ST(JJPR)
  ENDIF
ENDDO

DO JLEV=1,KLEV
  IF(PREHYD(JLEV) >= ZP_ST(2)) THEN
    IPRO_ST(JLEV)=1
  ELSEIF(PREHYD(JLEV) >= ZP_ST(3)) THEN
    IPRO_ST(JLEV)=2
  ELSEIF(PREHYD(JLEV) >= ZP_ST(4)) THEN
    IPRO_ST(JLEV)=3
  ELSEIF(PREHYD(JLEV) >= ZP_ST(5)) THEN
    IPRO_ST(JLEV)=4
  ELSEIF(PREHYD(JLEV) >= ZP_ST(6)) THEN
    IPRO_ST(JLEV)=5
  ELSEIF(PREHYD(JLEV) >= ZP_ST(7)) THEN
    IPRO_ST(JLEV)=6
  ELSEIF(PREHYD(JLEV) >= ZP_ST(8)) THEN
    IPRO_ST(JLEV)=7
  ELSEIF(PREHYD(JLEV) >= ZP_ST(9)) THEN
    IPRO_ST(JLEV)=8
  ELSE
    IPRO_ST(JLEV)=9
  ENDIF
ENDDO

ZROG=Z_RD/Z_RG

IF(.NOT.LLGRZER_ST) THEN
  DO JLEV=1,KLEV
    IPR=IPRO_ST(JLEV)
    PSTZ(JLEV)=ZZ_ST(IPR)+ZT_ST(IPR)*ZIDTDZ_ST(IPR)*((PREHYD(JLEV)*ZIP_ST(IPR))**(-ZDTDZ_ST(IPR)*ZROG)-1.0)
  ENDDO
ELSE
  DO JLEV=1,KLEV
    IPR=IPRO_ST(JLEV)
    IF(ZDTDZ_ST(IPR) /= 0.0) THEN
      PSTZ(JLEV)=ZZ_ST(IPR)+ZT_ST(IPR)*ZIDTDZ_ST(IPR)*((PREHYD(JLEV)*ZIP_ST(IPR))**(-ZDTDZ_ST(IPR)*ZROG)-1.0)
    ELSE
      PSTZ(JLEV)=ZZ_ST(IPR)-ZT_ST(IPR)*ZROG*(LOG(PREHYD(JLEV))-ZLNP_ST(IPR))
    ENDIF
  ENDDO
ENDIF

!     =======================================================================

RETURN
END SUBROUTINE ALTSTA

