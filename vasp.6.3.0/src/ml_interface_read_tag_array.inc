#if defined (TYPE_REAL)
#define TYPE_VAL REAL
#define TYPE_SUF REAL(q)
#elif defined(TYPE_INTEGER)
#define TYPE_VAL INTEGER
#define TYPE_SUF INTEGER
#elif defined(TYPE_LOGICAL)
#define TYPE_VAL LOGICAL
#define TYPE_SUF LOGICAL
#endif
#define HELP1(X)    X
#define HELP2(X)    HELP1(X)_
#define HELP3(X, Y) HELP2(X)Y
#define HELP4(X, Y) HELP3(X,Y)
#define PNAME       HELP4(READ_TAG_ARRAY,TYPE_VAL)
      SUBROUTINE PNAME(TAG, ARRAY, NA, DEFAULT_VALUE)
         IMPLICIT NONE

         !> INCAR tag to be processed.
         CHARACTER(len=*), INTENT(IN)    :: TAG
         !> Array where to store default or INCAR tag argument.
         TYPE_SUF,          INTENT(INOUT) :: ARRAY(:)
         !> Number of entries in variable array
         INTEGER,          INTENT(IN)    :: NA
         !> Default value for variable.
         TYPE_SUF,          INTENT(IN)    :: DEFAULT_VALUE

         !> Alternative tags
         CHARACTER(len=:), ALLOCATABLE   :: ALT_TAGS(:)
         INTEGER                         :: I

         ! Set default value.
         ARRAY = DEFAULT_VALUE
         CALL TAG_LIST%SET_STATE(TAG, TL_DEFAULT)

         ALT_TAGS = TAGLIST_GET_ALT_TAGS(TAG)

         DO I = 1, SIZE(ALT_TAGS)
            ! Ignore alternative tags if a higher priority tag was already
            ! found and issue a warning.
            ! Check for per-instance INCAR tag first.
            CALL PROCESS_INCAR(LOPEN,                                &
                               IU0,                                  &
                               IU5,                                  &
                               ADJUSTL(TRIM(PREFIX // ALT_TAGS(I))), &
                               ARRAY,                                &
                               NA,                                   &
                               IERR,                                 &
                               WRITEXMLINCAR,                        &
                               LCONTINUE=.TRUE.,                     &
                               FOUNDNUMBER=N)

            ! If not found, parse INCAR for regular occurence of tag.
            IF (N == 0) THEN
               CALL PROCESS_INCAR(LOPEN,                      &
                                  IU0,                        &
                                  IU5,                        &
                                  ADJUSTL(TRIM(ALT_TAGS(I))), &
                                  ARRAY,                      &
                                  NA,                         &
                                  IERR,                       &
                                  WRITEXMLINCAR,              &
                                  FOUNDNUMBER=N)
            ENDIF

            ! Skip alternative tags if one is already found.
            IF (N > 0) THEN
               EXIT
            ENDIF
         END DO

         ! If found anywhere, set tag list state accordingly.
         IF (N > 0) THEN
            IF (I > 1) THEN
               CALL TAG_LIST%SET_STATE(TAG, TL_INCAR_ALT, I)
            ELSE
               CALL TAG_LIST%SET_STATE(TAG, TL_INCAR, I)
            END IF
         END IF

      END SUBROUTINE PNAME
#undef HELP1
#undef HELP2
#undef HELP3
#undef HELP4
#undef PNAME
#undef TYPE_VAL
#undef TYPE_SUF
