!input functions
MODULE input_module
  USE globals
  USE errors_module
  USE string_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_cmd_args,read_files

  !The maximum length of a cardname/blockname
  INTEGER(ki4),PARAMETER :: MAX_CARDNAME_LEN=32
  !number of blocks
  INTEGER(ki4),PARAMETER :: num_blocks=3
  !The maximum length of a line in the input file
  INTEGER(ki4),PARAMETER :: ll_max=200
  !(also need to change in interface IFchanged)
  !The maximum number of params on a given line or inputs for a param
  INTEGER(ki4),PARAMETER :: lp_max=100
  !(also need to change in interface IFchanged)
  CHARACTER(ll_max) :: xs_in

  INTEGER(ki4) :: in_unit

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
    !number of cards in the block
    INTEGER(ki4) :: num_cards
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
      CHARACTER(200),INTENT(INOUT) :: twords(:)
    ENDSUBROUTINE prototype_wordarg
  ENDINTERFACE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_cmd_args()
    INTEGER(ki4) :: arg_count

    arg_count = COMMAND_ARGUMENT_COUNT()

    IF(arg_count .NE. 1)STOP 'only give input file argument!'

    CALL GET_COMMAND_ARGUMENT(1,base_in)
  ENDSUBROUTINE read_cmd_args

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_files()
    INTEGER(ki4) :: t_int
    CHARACTER(100) :: t_char

    in_unit=20

    !open input unit
    OPEN(UNIT=in_unit, FILE=base_in, STATUS='OLD', ACTION = "READ", IOSTAT=t_int, IOMSG=t_char)
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

    CLOSE(in_unit)

    !open input unit
    OPEN(UNIT=in_unit, FILE=xs_in, STATUS='OLD', ACTION = "READ", IOSTAT=t_int, IOMSG=t_char)
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
      IF(t_char .EQ. 'OpenNodal_ASSY-XS_V1')THEN
        !current OpenNodal_ASSY-XS_V1 input format
        REWIND(in_unit)
        CALL read_xs_inp_v1()
        EXIT
      ENDIF
    ENDDO

    CLOSE(in_unit)
  ENDSUBROUTINE read_files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_base_inp_v1()
    INTEGER(ki4) :: i,j,t_int,nwords
    CHARACTER(ll_max) :: t_char,words(lp_max)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !data for caseid block
    i=1
    blocks(i)%bname='[CASEID]'
    blocks(i)%num_cards=1
    ALLOCATE(blocks(i)%cards(blocks(i)%num_cards))
    j=1
    blocks(i)%cards(j)%cname='title'
    blocks(i)%cards(j)%getcard=>get_title

    !data for CORE block
    i=2
    blocks(i)%bname='[CORE]'
    blocks(i)%num_cards=5
    ALLOCATE(blocks(i)%cards(blocks(i)%num_cards))
    j=1
    blocks(i)%cards(j)%cname='dim'
    blocks(i)%cards(j)%getcard=>get_dim
    j=2
    blocks(i)%cards(j)%cname='size'
    blocks(i)%cards(j)%getcard=>get_size
    j=3
    blocks(i)%cards(j)%cname='apitch'
    blocks(i)%cards(j)%getcard=>get_apitch
    j=4
    blocks(i)%cards(j)%cname='sym'
    blocks(i)%cards(j)%getcard=>get_sym
    j=5
    blocks(i)%cards(j)%cname='assm_map'
    blocks(i)%cards(j)%getcard=>get_assm_map

    !data for MATERIAL block
    i=3
    blocks(i)%bname='[MATERIAL]'
    blocks(i)%num_cards=2
    ALLOCATE(blocks(i)%cards(blocks(i)%num_cards))
    j=1
    blocks(i)%cards(j)%cname='xs_file'
    blocks(i)%cards(j)%getcard=>get_xs_file
    j=2
    blocks(i)%cards(j)%cname='xs_map'
    blocks(i)%cards(j)%getcard=>get_xs_map
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO i=1,num_blocks
      REWIND(in_unit)
      CALL find_block(blocks(i)%bname)
      !loop through cards for this block
      DO j=1,blocks(i)%num_cards
        DO
          CALL get_next_line(t_char,t_int)
          !did not find card name, defaulting. Happens if either EOF or a new block
          IF(t_int .NE. 0)EXIT
          IF(t_char(1:1) .EQ. '[')EXIT
          !parse and check to see if card is found
          CALL parse(t_char,' ',words,nwords)
          words(1)=lowercase(words(1))
          IF(words(1) .EQ. blocks(i)%cards(j)%cname)THEN
            !found the card name
            blocks(i)%cards(j)%found=.TRUE.
            EXIT
          ENDIF
        ENDDO
        !if found, get the data
        IF(blocks(i)%cards(j)%found)THEN
          CALL blocks(i)%cards(j)%getcard(words(1:nwords))
        ENDIF
        !after finding or not finding, get back to the beginning of the block
        REWIND(in_unit)
        CALL find_block(blocks(i)%bname)
      ENDDO
    ENDDO
  ENDSUBROUTINE read_base_inp_v1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_xs_inp_v1()
    use edits_module, only : edit_xs
    CHARACTER(ll_max) :: t_char,words(lp_max)
    INTEGER(ki4) :: ios,nwords,num_local_xs,i,j

    !get info on the number of materials and number of energy groups
    REWIND(in_unit)
    CALL get_next_line(t_char,ios)
    CALL parse(t_char,' ',words,nwords)
    READ(words(2),*)num_local_xs
    READ(words(3),*)num_eg

    DO i=1,num_assm_reg
      ALLOCATE(assm_xs(i)%D(num_eg))
      assm_xs(i)%D=0.0
      ALLOCATE(assm_xs(i)%chi(num_eg))
      assm_xs(i)%chi=0.0
      ALLOCATE(assm_xs(i)%nusigma_f(num_eg))
      assm_xs(i)%nusigma_f=0.0
      ALLOCATE(assm_xs(i)%nu(num_eg))
      assm_xs(i)%nu=0.0
      ALLOCATE(assm_xs(i)%sigma_a(num_eg))
      assm_xs(i)%sigma_a=0.0
      ALLOCATE(assm_xs(i)%sigma_scat(num_eg,num_eg))
      assm_xs(i)%sigma_scat=0.0
    ENDDO

    !get xs data for each of the assembly regions
    DO i=1,num_assm_reg
      REWIND(in_unit)
      DO
        CALL get_next_line(t_char,ios)
        IF(ios .NE. 0)CALL fatal_error('Could not find xs for material "'//TRIM(assm_xs(i)%mat_id)//'"')
        CALL parse(t_char,' ',words,nwords)
        words(1)=lowercase(words(1))
        !check to see if it's a start of a new xs
        IF(words(1) .EQ. 'id')THEN
          !check to see if it's the correct xs
          IF(words(2) .EQ. assm_xs(i)%mat_id)THEN
            !get d
            CALL get_next_line(t_char,ios)
            BACKSPACE(in_unit)
            READ(in_unit,*)assm_xs(i)%D(:)
            !get chi
            CALL get_next_line(t_char,ios)
            BACKSPACE(in_unit)
            READ(in_unit,*)assm_xs(i)%chi(:)
            !get numsigma_f
            CALL get_next_line(t_char,ios)
            BACKSPACE(in_unit)
            READ(in_unit,*)assm_xs(i)%nusigma_f(:)
            !get nu
            CALL get_next_line(t_char,ios)
            BACKSPACE(in_unit)
            READ(in_unit,*)assm_xs(i)%nu(:)
            !get sigma_a
            CALL get_next_line(t_char,ios)
            BACKSPACE(in_unit)
            READ(in_unit,*)assm_xs(i)%sigma_a(:)
            !get sigma_scat
            DO j=1,num_eg
              CALL get_next_line(t_char,ios)
              BACKSPACE(in_unit)
              READ(in_unit,*)assm_xs(i)%sigma_scat(j,:)
            ENDDO
            EXIT
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    call edit_xs()

  ENDSUBROUTINE read_xs_inp_v1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE find_block(block_name)
    CHARACTER(MAX_CARDNAME_LEN),INTENT(IN) :: block_name

    CHARACTER(ll_max) :: t_char
    INTEGER(ki4) :: ios,i

    !find the next block
    DO
      CALL get_next_line(t_char,ios)
      t_char=uppercase(t_char)
      !found the bloc
      IF(t_char .EQ. block_name)EXIT
      IF(ios .NE. 0)CALL fatal_error('Could not find block: '//block_name)
      !make sure all blocks are actual blocks
      IF(t_char(1:1) .EQ. '[')THEN
        DO i=1,num_blocks
          IF(t_char .EQ. blocks(i)%bname)EXIT
        ENDDO
        IF(i .GT. num_blocks)CALL fatal_error('Unrecognized block: "'//TRIM(t_char)//'"')
      ENDIF
    ENDDO
  ENDSUBROUTINE find_block

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_next_line(line,ios)
    CHARACTER(ll_max),INTENT(OUT) :: line
    INTEGER(ki4),INTENT(OUT) :: ios

    CHARACTER(ll_max) :: words(lp_max)
    INTEGER(ki4) :: nwords

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_title(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(:)

    INTEGER(ki4) :: len_str,i

    CALL print_log(TRIM(this_card%cname)//' card found')

    len_str=SIZE(wwords)

    prob_title=''
    DO i=2,len_str
      prob_title=TRIM(prob_title)//' '//TRIM(wwords(i))
    ENDDO
    prob_title=TRIM(ADJUSTL(prob_title))
  ENDSUBROUTINE get_title

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_dim(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(:)

    CALL print_log(TRIM(this_card%cname)//' card found')

    IF(wwords(2) .EQ. '2D')THEN
      prob_dim=2
    ELSE
      CALL fatal_error(TRIM(wwords(2))//' is not a valid dimension.')
    ENDIF
  ENDSUBROUTINE get_dim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_size(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(:)

    CALL print_log(TRIM(this_card%cname)//' card found')

    READ(wwords(2),*)core_size
  ENDSUBROUTINE get_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_apitch(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(:)

    CALL print_log(TRIM(this_card%cname)//' card found')

    READ(wwords(2),*)assm_pitch
  ENDSUBROUTINE get_apitch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_sym(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(:)

    CALL print_log(TRIM(this_card%cname)//' card found')

    prob_sym=TRIM(lowercase(wwords(2)))
  ENDSUBROUTINE get_sym

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_assm_map(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(:)

    INTEGER(ki4) :: nwords
    CHARACTER(ll_max) :: t_char,words(lp_max)
    INTEGER(ki4) :: i,ios,j,oct_sym

    wwords(1)=TRIM(wwords(1))
    CALL print_log(TRIM(this_card%cname)//' card found')

    !allocate the assembly map based upon problem size and symmetry
    !this is the actual problem we will solve, and remember again that the core is assumed square
    ! TODO need to implement ragged core
    SELECTCASE(prob_sym)
      CASE('full')
        core_x_size=core_size
        core_y_size=core_size
      CASE('half')
        core_x_size=CEILING((core_size-1.0D-2)/2.0)
        core_y_size=core_size
      CASE('qtr')
        core_x_size=CEILING((core_size-1.0D-2)/2.0)
        core_y_size=CEILING((core_size-1.0D-2)/2.0)
      CASE DEFAULT
        CALL fatal_error(TRIM(prob_sym)//' is not a valid symmetry')
    ENDSELECT
    ALLOCATE(assm_map(core_x_size,core_y_size))
    assm_map=0

    !read in the core map
    DO i=1,core_y_size
      CALL get_next_line(t_char,ios)
      IF(ios .NE. 0)CALL fatal_error('Given assembly map incomplete!')
      !parse to get all data on line
      CALL parse(t_char,' ',words,nwords)
      DO j=1,nwords
        READ(words(j),*)assm_map(i,j)
      ENDDO
    ENDDO

    !correct if given an octant symmetry
    IF(prob_sym .EQ. 'qtr')THEN
      oct_sym=0
      DO i=1,core_y_size
        DO j=i+1,core_x_size
          oct_sym=oct_sym+ABS(assm_map(i,j))
        ENDDO
      ENDDO
      IF(oct_sym .EQ. 0)THEN
        DO i=1,core_y_size
          DO j=i+1,core_x_size
            assm_map(i,j)=assm_map(j,i)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDSUBROUTINE get_assm_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_xs_file(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(:)

    CALL print_log(TRIM(this_card%cname)//' card found')

    xs_in=wwords(2)
  ENDSUBROUTINE get_xs_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_xs_map(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(:)

    INTEGER(ki4) :: i,ios,nwords
    CHARACTER(ll_max) :: t_char,words(lp_max)

    CALL print_log(TRIM(this_card%cname)//' card found')

    READ(wwords(2),*)num_assm_reg
    IF(MAXVAL(assm_map) .GT. num_assm_reg)CALL fatal_error('Assembly map index outside of range of xs.')
    ALLOCATE(assm_xs(num_assm_reg))

    DO i=1,num_assm_reg
      CALL get_next_line(t_char,ios)
      IF(ios .NE. 0)CALL fatal_error('Given cross section map incomplete!')
      CALL parse(t_char,' ',words,nwords)
      IF(nwords .LT. 3)CALL fatal_error('Each xs map data should have 3 entries!')
      READ(words(1),*)ios
      IF(ios .NE. i)CALL fatal_error('xs map should be in order and complete!')
      IF(words(2) .EQ. 'macro')THEN
        assm_xs(i)%mat_id=TRIM(words(3))
      ELSE
        CALL fatal_error(TRIM(words(2))//' xs not yet supported.')
      ENDIF
    ENDDO
  ENDSUBROUTINE get_xs_map
ENDMODULE input_module
