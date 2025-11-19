MODULE DRVBUF_DAT
!
! DESCRIPTION
!   Holds in-memory driver table lines supplied from Python so that the
!   standard RFM input routines can operate without an on-disk ``rfm.drv``.
!
    USE KIND_DAT
    USE LENREC_DAT
!
  IMPLICIT NONE
!
  LOGICAL          :: DRVBUF_ENABLED = .FALSE.
  INTEGER(I4)      :: DRVBUF_COUNT   = 0
  CHARACTER(LENREC), ALLOCATABLE :: DRVBUF_LINES(:)
!
CONTAINS
!
  SUBROUTINE DRVBUF_SET ( LINES )
    CHARACTER(*), INTENT(IN) :: LINES(:)
    INTEGER(I4) :: I, J, L
    CHARACTER(LENREC) :: CLEAN
!
! Copy Python strings through a LENREC buffer, replacing any embedded NULs.
    CALL DRVBUF_CLEAR()
    IF ( SIZE(LINES) .EQ. 0 ) RETURN
!
    ALLOCATE ( DRVBUF_LINES(SIZE(LINES)) )
    DO I = 1, SIZE(LINES)
      CLEAN = ' '
      L = MIN ( LEN ( LINES(I) ), LENREC )
      DO J = 1, L
        IF ( IACHAR ( LINES(I)(J:J) ) .EQ. 0 ) THEN
          CLEAN(J:J) = ' '
        ELSE
          CLEAN(J:J) = LINES(I)(J:J)
        END IF
      END DO
      DRVBUF_LINES(I) = CLEAN
    END DO
!
    DRVBUF_COUNT   = SIZE ( LINES )
    DRVBUF_ENABLED = .TRUE.
  END SUBROUTINE DRVBUF_SET
!
  SUBROUTINE DRVBUF_CLEAR()
    IF ( ALLOCATED ( DRVBUF_LINES ) ) DEALLOCATE ( DRVBUF_LINES )
    DRVBUF_COUNT   = 0
    DRVBUF_ENABLED = .FALSE.
  END SUBROUTINE DRVBUF_CLEAR
!
END MODULE DRVBUF_DAT
