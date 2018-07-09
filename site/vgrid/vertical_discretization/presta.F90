SUBROUTINE PRESTA(KLEV,PSTZ,PREHYD)

!     ------------------------------------------------------------------
!     Computes standard pressure from standard altitude.
!     Code taken from SUSTA and PPSTA and adapted.

!     Input:
!      KLEV    : number of levels
!      PSTZ    : standard altitude

!     Output:
!      PREHYD  : standard pressure
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER*8,INTENT(IN) :: KLEV
REAL*8,INTENT(IN)    :: PSTZ(KLEV)
REAL*8,INTENT(OUT)   :: PREHYD(KLEV)

!     ------------------------------------------------------------------

INTEGER*8,PARAMETER :: JPLAY=10
REAL*8 :: ZPP_ST(JPLAY,2)
REAL*8 :: ZEPSZ
INTEGER*8 :: IITER,JITER
INTEGER*8 :: JLEV
REAL*8 :: ZSTZDKM(KLEV)
REAL*8 :: ZSTZ(KLEV)
INTEGER*8 :: ISTZDKM(KLEV)
REAL*8 :: ZSTP1(KLEV),ZSTP2(KLEV)
REAL*8 :: ZSTZ1(KLEV),ZSTZ2(KLEV)
REAL*8 :: ZCOEF(KLEV)

!     =======================================================================

! * ZPP_ST:
!   ZPP_ST(n,1) and ZPP_ST(n,2) are possible pressure bounds,
!   where n=int(z in 10km)+1
!   example:
!    z=1500m: possible pressure bounds are given by ZPP_ST(1,1) and ZPP_ST(1,2)
!    z=88.5km: possible pressure bounds are given by ZPP_ST(9,1) and ZPP_ST(9,2)
!   ZPP_ST(n,1) < ZPP_ST(n,2)

ZPP_ST( 1,1)= 26436.75
ZPP_ST( 2,1)=  5475.25
ZPP_ST( 3,1)=  1171.75
ZPP_ST( 4,1)=   276.75
ZPP_ST( 5,1)=    75.25
ZPP_ST( 6,1)=    20.25
ZPP_ST( 7,1)=     4.25
ZPP_ST( 8,1)=     0.75
ZPP_ST( 9,1)=     0.001
ZPP_ST(10,1)=     0.001

ZPP_ST( 1,2)=101325.
ZPP_ST( 2,2)= 26437.25
ZPP_ST( 3,2)=  5475.75
ZPP_ST( 4,2)=  1172.25
ZPP_ST( 5,2)=   277.25
ZPP_ST( 6,2)=    75.75
ZPP_ST( 7,2)=    20.75
ZPP_ST( 8,2)=     4.75
ZPP_ST( 9,2)=     1.25
ZPP_ST(10,2)=     1.25

!     =======================================================================

!     * Iterations:

IITER=20
ZEPSZ=0.00001

! - first guess:
DO JLEV=1,KLEV
  ZSTZDKM(JLEV) = PSTZ(JLEV)/10000.
  ISTZDKM(JLEV) = INT(ZSTZDKM(JLEV))+1
  IF (ISTZDKM(JLEV) > 10) ISTZDKM(JLEV)=10
  ZSTP1(JLEV) = ZPP_ST(ISTZDKM(JLEV),1)
  ZSTP2(JLEV) = ZPP_ST(ISTZDKM(JLEV),2)
ENDDO
CALL ALTSTA(KLEV,ZSTP1,ZSTZ1)
CALL ALTSTA(KLEV,ZSTP2,ZSTZ2)
DO JLEV=1,KLEV
  ZCOEF(JLEV)=(PSTZ(JLEV)-ZSTZ1(JLEV))/(ZSTZ2(JLEV)-ZSTZ1(JLEV))
  PREHYD(JLEV)=ZSTP1(JLEV)+ZCOEF(JLEV)*(ZSTP2(JLEV)-ZSTP1(JLEV))
ENDDO
CALL ALTSTA(KLEV,PREHYD,ZSTZ)

! - following iterations:
DO JITER=1,IITER
 DO JLEV=1,KLEV
  IF (ZSTZ(JLEV) - PSTZ(JLEV) > ZEPSZ) THEN
    ZSTP2(JLEV) = PREHYD(JLEV)
  ELSEIF (ZSTZ(JLEV) - PSTZ(JLEV) < - ZEPSZ) THEN
    ZSTP1(JLEV) = PREHYD(JLEV)
  ENDIF
 ENDDO
 CALL ALTSTA(KLEV,ZSTP1,ZSTZ1)
 CALL ALTSTA(KLEV,ZSTP2,ZSTZ2)
 DO JLEV=1,KLEV
  IF (ABS(ZSTZ2(JLEV)-ZSTZ1(JLEV)) > ZEPSZ) THEN
   ZCOEF(JLEV)=(PSTZ(JLEV)-ZSTZ1(JLEV))/(ZSTZ2(JLEV)-ZSTZ1(JLEV))
   PREHYD(JLEV)=ZSTP1(JLEV)+ZCOEF(JLEV)*(ZSTP2(JLEV)-ZSTP1(JLEV))
  ELSE
   PREHYD(JLEV)=ZSTP1(JLEV)
  ENDIF
 ENDDO
 CALL ALTSTA(KLEV,PREHYD,ZSTZ)
ENDDO

! - check that algorithm has converged:
DO JLEV=1,KLEV
  IF (ZSTZ(JLEV) - PSTZ(JLEV) > ZEPSZ) THEN
    write(6,*) ' for jlev=',jlev,' PREHYD is too low'
    write(6,*) ' increase iiter'
  ELSEIF (ZSTZ(JLEV) - PSTZ(JLEV) < - ZEPSZ) THEN
    write(6,*) ' for jlev=',jlev,' PREHYD is too high'
    write(6,*) ' increase iiter'
  ENDIF
ENDDO

!     =======================================================================

RETURN
END SUBROUTINE PRESTA

