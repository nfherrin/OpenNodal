!a module for string parsing
MODULE string_module
  USE precision_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: parse

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE parse(str,delims,args,nargs)
    CHARACTER(*),INTENT(INOUT) :: str
    CHARACTER(*),INTENT(IN) :: delims
    CHARACTER(*),INTENT(OUT) :: args(:)
    INTEGER,INTENT(OUT) :: nargs

    CHARACTER(LEN_TRIM(str)) :: strsav
    INTEGER :: na,i,lenstr

    strsav=str
    CALL compact(str)
    na=SIZE(args)
    DO i=1,na
      args(i)=' '
    ENDDO
    nargs=0
    lenstr=LEN_TRIM(str)
    IF(lenstr==0)RETURN

    DO
       IF(LEN_TRIM(str) == 0)EXIT
       nargs=nargs+1
       CALL split(str,delims,args(nargs))
       CALL removebksl(args(nargs))
    ENDDO
    str=strsav
  ENDSUBROUTINE parse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE compact(str)
    CHARACTER(*),INTENT(INOUT) :: str

    CHARACTER(1):: ch
    CHARACTER(LEN_TRIM(str)) :: outstr

    INTEGER :: ich,isp,i,k,lenstr

    str=ADJUSTL(str)
    lenstr=LEN_TRIM(str)
    outstr=' '
    isp=0
    k=0

    DO i=1,lenstr
      ch=str(i:i)
      ich=IACHAR(ch)

      SELECTCASE(ich)
        CASE(9,32)     ! space or tab CHARACTER
          IF(isp==0) THEN
            k=k+1
            outstr(k:k)=' '
          ENDIF
          isp=1
        CASE(33:)      ! not a space, quote, or control CHARACTER
          k=k+1
          outstr(k:k)=ch
          isp=0
      ENDSELECT
    ENDDO

    str=ADJUSTL(outstr)
  ENDSUBROUTINE compact

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE split(str,delims,before,sep)
    CHARACTER(*),INTENT(INOUT) :: str
    CHARACTER(*),INTENT(IN) :: delims
    CHARACTER(*),INTENT(OUT) :: before
    CHARACTER,INTENT(INOUT),OPTIONAL :: sep

    LOGICAL :: pres
    CHARACTER :: ch,cha
    INTEGER :: ipos,iposa,i,ibsl,k,lenstr

    pres=PRESENT(sep)
    str=ADJUSTL(str)
    CALL compact(str)
    lenstr=LEN_TRIM(str)
    IF(lenstr == 0) RETURN        ! string str is empty
    k=0
    ibsl=0                        ! backslash initially inactive
    before=' '
    DO i=1,lenstr
       ch=str(i:i)
       IF(ibsl == 1) THEN          ! backslash active
          k=k+1
          before(k:k)=ch
          ibsl=0
          CYCLE
       ENDIF
       IF(ch == '\') THEN          ! backslash with backslash inactive
          k=k+1
          before(k:k)=ch
          ibsl=1
          CYCLE
       ENDIF
       ipos=index(delims,ch)
       IF(ipos == 0) THEN          ! CHARACTER is not a delimiter
          k=k+1
          before(k:k)=ch
          CYCLE
       ENDIF
       IF(ch /= ' ') THEN          ! CHARACTER is a delimiter that is not a space
          str=str(i+1:)
          IF(pres) sep=ch
          EXIT
       ENDIF
       cha=str(i+1:i+1)            ! CHARACTER is a space delimiter
       iposa=index(delims,cha)
       IF(iposa > 0) THEN          ! next CHARACTER is a delimiter
          str=str(i+2:)
          IF(pres) sep=cha
          EXIT
       else
          str=str(i+1:)
          IF(pres) sep=ch
          EXIT
       ENDIF
    ENDDO
    IF(i >= lenstr) str=''
    str=ADJUSTL(str)              ! remove initial spaces
  ENDSUBROUTINE split

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE removebksl(str)
    CHARACTER(*),INTENT(INOUT) :: str

    CHARACTER(1):: ch
    CHARACTER(LEN_TRIM(str))::outstr
    INTEGER :: i,ibsl,k,lenstr

    str=ADJUSTL(str)
    lenstr=LEN_TRIM(str)
    outstr=' '
    k=0
    ibsl=0                        ! backslash initially inactive

    DO i=1,lenstr
      ch=str(i:i)
      IF(ibsl == 1) THEN          ! backslash active
       k=k+1
       outstr(k:k)=ch
       ibsl=0
       CYCLE
      ENDIF
      IF(ch == '\') THEN          ! backslash with backslash inactive
       ibsl=1
       CYCLE
      ENDIF
      k=k+1
      outstr(k:k)=ch              ! non-backslash with backslash inactive
    ENDDO

    str=ADJUSTL(outstr)
  ENDSUBROUTINE removebksl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CHARACTER(64) FUNCTION str(k)
    !   "Convert an INTEGER to string."
    INTEGER, INTENT(IN) :: k

    WRITE (str, *) k
    str = ADJUSTL(str)
  ENDFUNCTION str
ENDMODULE string_module
