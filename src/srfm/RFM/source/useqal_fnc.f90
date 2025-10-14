MODULE USEQAL_FNC
CONTAINS
LOGICAL PURE FUNCTION USEQAL ( IDM, ISO, ILS, IUS )
!
! VERSION
!   27JUL24 AD Checked.
!   01MAY17 AD F90 original. Checked.
!
! DESCRIPTION
!   Set TRUE if listed in line molecule qualifiers.
!   Called by GASQAL, REAHIT
!   Values of ISO,ILS,IUS=0 act as wildcards.
!
! VARIABLE KINDS
    USE KIND_DAT
!
! GLOBAL DATA
    USE QALCOM_DAT ! Band/isotope selection qualifiers
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4), INTENT(IN) :: IDM ! HITRAN/RFM index of molecule
    INTEGER(I4), INTENT(IN) :: ISO ! HITRAN isotope#
    INTEGER(I4), INTENT(IN) :: ILS ! Lower vib level
    INTEGER(I4), INTENT(IN) :: IUS ! Upper vib level
!
! LOCAL VARIABLES
    INTEGER(I4) :: IQAL ! Counter for tabulated data
!
! EXECUTABLE CODE --------------------------------------------------------------
!
  DO IQAL = 1, NQAL
    IF ( QAL(IQAL)%IDM .EQ. IDM &
         .AND. &
         ( ISO .EQ. 0 .OR. QAL(IQAL)%ISO .EQ. 0 .OR. ISO .EQ. QAL(IQAL)%ISO ) &
         .AND. &
         ( ILS .EQ. 0 .OR. QAL(IQAL)%ILS .EQ. 0 .OR. ILS .EQ. QAL(IQAL)%ILS ) &
         .AND. &
         ( IUS .EQ. 0 .OR. QAL(IQAL)%IUS .EQ. 0 .OR. IUS .EQ. QAL(IQAL)%IUS ) &
                  ) THEN
      USEQAL = .TRUE.
      RETURN
    END IF
  END DO
!
  USEQAL = .FALSE.
!
END FUNCTION USEQAL
END MODULE USEQAL_FNC


