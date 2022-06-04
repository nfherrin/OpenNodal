!input functions
MODULE input_module
  USE globals
  USE errors_module
  USE string_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_cmd_args,read_base_input

  !The maximum length of a cardname/blockname
  INTEGER,PARAMETER :: MAX_CARDNAME_LEN=32
  !number of blocks
  INTEGER,PARAMETER :: num_blocks=3
  !The maximum length of a line in the input file
  INTEGER,PARAMETER :: ll_max=200
  !(also need to change in interface IFchanged)
  !The maximum number of params on a given line or inputs for a param
  INTEGER,PARAMETER :: lp_max=100
  !(also need to change in interface IFchanged)

  INTEGER :: in_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE :: cardType
    !name of the card
    CHARACTER(MAX_CARDNAME_LEN) :: cname
    !logical to tell if block has already appeared
    LOGICAL :: found=.FALSE.
    !readin procedure, unique for each card
    PROCEDURE(prototype_wordarg),POINTER :: getcard => NULL()
  ENDTYPE cardType

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE :: blockType
    !name of the block
    CHARACTER(MAX_CARDNAME_LEN) :: bname
    !logical to tell if block has already appeared
    LOGICAL :: found=.FALSE.
    !number of cards in the block
    INTEGER :: num_cards
    !cards in the block
    TYPE(cardType),ALLOCATABLE :: cards(:)
  ENDTYPE blockType

  !actual blocks variables
  TYPE(blockType) :: blocks(num_blocks)

  !Simple abstract interface for an object method subroutine with word array argument
  ABSTRACT INTERFACE
    SUBROUTINE prototype_wordarg(thisCard,twords)
      IMPORT :: cardType
      CLASS(cardType),INTENT(INOUT) :: thisCard
      CHARACTER(200),INTENT(INOUT) :: twords(100)
    ENDSUBROUTINE prototype_wordarg
  ENDINTERFACE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_cmd_args()
    INTEGER :: arg_count

    arg_count = COMMAND_ARGUMENT_COUNT()

    IF(arg_count .NE. 1)STOP 'only give input file argument!'

    CALL GET_COMMAND_ARGUMENT(1,base_in)
  ENDSUBROUTINE read_cmd_args

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_base_input()
    INTEGER :: t_int
    CHARACTER(100) :: t_char

    in_unit=20

    !open input unit
    OPEN(UNIT=in_unit, FILE=base_in, STATUS='OLD', ACTION = "read", IOSTAT=t_int, IOMSG=t_char)
    IF(t_int .NE. 0)THEN
      CALL fatal_error(t_char)
    ENDIF

    !find version ID
    DO
      READ(in_unit,*,IOSTAT=t_int)t_char
      !if we don't find a recognized ID by the end of the file, error out
      IF(t_int .NE. 0)THEN
        CALL fatal_error('No recognized format identifiers in base input')
      ENDIF
      t_char=TRIM(ADJUSTL(t_char))
      IF(t_char .EQ. 'OpenNodal_INP_V1')THEN
        !current OpenNodal_INP_V1 input format
        REWIND(in_unit)
        CALL read_base_inp_v1()
        EXIT
      ENDIF
    ENDDO
  ENDSUBROUTINE read_base_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_base_inp_v1()
    INTEGER :: i

    !data for caseid block
    i=1
    blocks(i)%bname='[CASEID]'

    !data for CORE block
    i=2
    blocks(i)%bname='[CORE]'

    !data for MATERIAL block
    i=3
    blocks(i)%bname='[MATERIAL]'

    DO i=1,num_blocks
      REWIND(in_unit)
      CALL find_block(blocks(i)%bname)
    ENDDO
  ENDSUBROUTINE read_base_inp_v1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE find_block(block_name)
    CHARACTER(MAX_CARDNAME_LEN),INTENT(IN) :: block_name

    CHARACTER(ll_max) :: t_char
    INTEGER :: ios

    !find the next block
    DO
      CALL get_next_line(t_char,ios)
      IF(t_char .EQ. block_name)EXIT
      IF(ios .NE. 0)CALL fatal_error('Could not find block: '//block_name)
    ENDDO
  ENDSUBROUTINE find_block

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_next_line(line,ios)
    CHARACTER(ll_max),INTENT(OUT) :: line
    INTEGER,INTENT(OUT) :: ios

    CHARACTER(ll_max) :: words(100)
    INTEGER :: nwords

    DO
      READ(in_unit,'(A10000)',IOSTAT=ios)line
      IF(ios .NE. 0)EXIT
      line=TRIM(ADJUSTL(line))
      !finding uncommented line that isn't empty
      IF(line(1:1) .NE. '!' .AND. line .NE. '')THEN
        !ignore commented portions of line
        CALL parse(line,'!',words,nwords)
        line=TRIM(ADJUSTL(words(1)))
        EXIT
      ENDIF
    ENDDO
  ENDSUBROUTINE get_next_line
ENDMODULE input_module
