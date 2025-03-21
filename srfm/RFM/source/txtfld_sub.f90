MODULE TXTFLD_SUB
CONTAINS
PURE SUBROUTINE TXTFLD ( RECORD, NFIELD, ISTA, IEND )
!
! VERSION
!   18DEC23 AD Checked.
!   01MAY17 AD F90 conversion. Checked.
!
! DESCRIPTION
!   Identify start and end points of text field in record
!   General purpose module.
!   Given a character string containing a number of alpha-numeric fields
!   separated by commas or spaces, this function returns the location of the
!   first and last characters of the specified field. 
!   Returns a value 0 if reqd field# is not found before the end of the 
!   string. End value limited to record length.
!
!   The following characters are regarded as spaces:
!     Quote Marks ['] (allows text strings to be optionally enclosed in ' ')
!     Commas [,]      (to conform to free format rules)
!     Tabs ( ICHAR( ) = 9 )
!     <CR> ( ICHAR( ) = 13 )
!   Other punctuation: treated as characters.
!
! VARIABLE KINDS
    USE KIND_DAT
!
  IMPLICIT NONE
!
! ARGUMENTS
    CHARACTER(*), INTENT(IN)  :: RECORD ! ASCII Record to be searched 
    INTEGER(I4),  INTENT(IN)  :: NFIELD ! Field# required
    INTEGER(I4),  INTENT(OUT) :: ISTA   ! Start location of field
    INTEGER(I4),  INTENT(OUT) :: IEND   ! End location of field
!
! LOCAL VARIABLES
    LOGICAL     :: LFIELD ! TRUE if pointer within a field
    LOGICAL     :: LSPACE ! TRUE if pointer in space between fields
    INTEGER(I4) :: I,J    ! Pointers within record
    INTEGER(I4) :: IFIELD ! Field# identified so far
    INTEGER(I4) :: ILEN   ! Effective record length up to ! mark
!
! EXECUTABLE CODE -------------------------------------------------------------
!
! Initialise at start of each call
  IFIELD = 0
  LFIELD = .FALSE.
!
! Check for any exclamation marks and truncate effective record length
  ILEN = INDEX ( RECORD//'!', '!' ) - 1
!
! Now identify text fields
  DO I = 1, ILEN
    LSPACE = ( RECORD(I:I) .EQ. ' '  .OR.  &         ! = space character
               RECORD(I:I) .EQ. ','  .OR.  &         ! = comma
               RECORD(I:I) .EQ. '''' .OR.  &         ! = single quote mark
        ICHAR(RECORD(I:I)) .EQ. 9    .OR.  &         ! = tab
        ICHAR(RECORD(I:I)) .EQ. 13   )               ! = <CR>
    IF ( .NOT. LFIELD .AND. .NOT. LSPACE ) THEN      ! found start of field
      IFIELD = IFIELD + 1
      IF ( IFIELD .EQ. NFIELD ) THEN                 ! start of reqd field#
        ISTA = I
        DO J = I+1, ILEN
          LSPACE = ( RECORD(J:J) .EQ. ' '  .OR. &    ! = space character
                     RECORD(J:J) .EQ. ','  .OR. &    ! = comma
                     RECORD(J:J) .EQ. '''' .OR. &    ! = single quote mark
              ICHAR(RECORD(J:J)) .EQ. 9    .OR. &    ! = tab
              ICHAR(RECORD(J:J)) .EQ. 13 )           ! = <CR>
          IF ( LSPACE ) THEN
            IEND = J - 1
            RETURN          ! found field within RECORD
          END IF
        END DO
        IEND = ILEN          
        RETURN              ! field limited by end of RECORD
      END IF
    END IF
    LFIELD = .NOT. LSPACE 
  END DO
!
! If the program gets this far, the reqd field# was not found in the record
  ISTA = 0
  IEND = 0
!
END SUBROUTINE TXTFLD
END MODULE TXTFLD_SUB

