!OpenNodal is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Module for input routines.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE input_module
  USE globals
  USE errors_module
  USE string_module
  USE edits_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_cmd_args,read_files

  !> The maximum length of a cardname/blockname
  INTEGER(ki4),PARAMETER :: MAX_CARDNAME_LEN=32
  !> number of blocks
  INTEGER(ki4),PARAMETER :: num_blocks=3
  !> The maximum length of a line in the input file
  !> (also need to change in interface IFF changed)
  INTEGER(ki4),PARAMETER :: ll_max=200
  !> The maximum number of params on a given line or inputs for a param
  !> (also need to change in interface IFF changed)
  INTEGER(ki4),PARAMETER :: lp_max=100
  !> cross section input filename
  CHARACTER(ll_max) :: xs_in
  !> In file unit
  INTEGER(ki4) :: in_unit

  !variables to be used to gather data for the variables in main, also the defaults
  INTEGER(ki4) :: p_d=2                           !prob_dim
  INTEGER(ki4) :: c_s=0                           !core_size
  INTEGER(ki4) :: x_s=0                           !core_x_size
  INTEGER(ki4) :: y_s=0                           !core_y_size
  REAL(kr8) :: a_p=0                              !assm_pitch
  CHARACTER(100) :: p_s='full'                    !prob_sym
  INTEGER(ki4), ALLOCATABLE :: a_m(:,:)           !assm_map
  REAL(kr8), ALLOCATABLE :: d_x(:)                !h_x
  REAL(kr8), ALLOCATABLE :: d_y(:)                !h_y
  CHARACTER(100) :: b_o='vacuum'                  !bc_opt
  REAL(kr8), ALLOCATABLE :: alb(:)                !albedos
  INTEGER(ki4) :: n_g                             !num_eg
  INTEGER(ki4) :: ns=1                            !nsplit
  REAL(kr8) :: t_xk = 1d-6                        !tol_xkeff
  REAL(kr8) :: t_xf = 1d-5                        !tol_xflux
  INTEGER(ki4) :: t_m_i = 100                     !tol_max_iter
  CHARACTER(100) :: n_m='fd'                      !nodal_method
  INTEGER(ki4) :: n_a_r=1                         !num_assm_reg
  TYPE(macro_assm_xs_type), ALLOCATABLE :: a_x(:) !assm_xs
  INTEGER(ki4) :: r_m=0                           !refl_mat
  REAL(kr8) :: a_b=0.0D0                          !ax_buckle

  !> This data type stores card information and reads card data
  TYPE :: cardType
    !> name of the card
    CHARACTER(MAX_CARDNAME_LEN) :: cname
    !> logical to tell if block has already appeared
    LOGICAL :: found=.FALSE.
    !> readin procedure, unique for each card
    PROCEDURE(prototype_wordarg), POINTER :: getcard => NULL()
  ENDTYPE cardType

  !> This data type stores card information
  TYPE :: blockType
    !> name of the block
    CHARACTER(MAX_CARDNAME_LEN) :: bname
    !> number of cards in the block
    INTEGER(ki4) :: num_cards
    !> cards in the block
    TYPE(cardType), ALLOCATABLE :: cards(:)
  ENDTYPE blockType

  !> actual blocks variables
  TYPE(blockType) :: blocks(num_blocks)

!---------------------------------------------------------------------------------------------------
!> Explicitly defines the interface for an object method subroutine with word array argument
!> @param thisCard - The card we're retrieving data for
!> @param twords - The string from which the data is being retrieved
!>
  ABSTRACT INTERFACE
    SUBROUTINE prototype_wordarg(thisCard,twords)
      IMPORT :: cardType
      CLASS(cardType), INTENT(INOUT) :: thisCard
      CHARACTER(200), INTENT(INOUT) :: twords(:)
    ENDSUBROUTINE prototype_wordarg
  ENDINTERFACE

CONTAINS

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine reads in command line arguments
!>
  SUBROUTINE read_cmd_args()
    INTEGER(ki4) :: arg_count

    arg_count = COMMAND_ARGUMENT_COUNT()

    IF(arg_count .NE. 1)STOP 'only give input file argument!'

    CALL GET_COMMAND_ARGUMENT(1,base_in)
  ENDSUBROUTINE read_cmd_args

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine reads in all input files to get input data
!> @param prob_dim - problem dimension
!> @param core_x_size - core size in the x direction
!> @param core_y_size - core size in the y direction
!> @param assm_pitch - assembly pitch
!> @param prob_sym - problem symmetry
!> @param assm_map - assembly map
!> @param h_x - node widths in the x direction
!> @param h_y - node widths in the y direction
!> @param bc_opt - boundary condition option
!> @param albedos - albedo boundary conditions
!> @param num_eg - number of energy groups
!> @param nsplit - nsplit value, for decomposing nodes into nsplit number of sub-nodes
!> @param tol_xkeff - keff convergence tolerance
!> @param tol_xflux - flux convergence tolerance
!> @param tol_max_iter - maximum number of iterations
!> @param nodal_method - nodal method option
!> @param num_assm_reg - number of unique assemblies
!> @param assm_xs - assembly level cross sections
!> @param refl_mat - reflector material
!> @param ax_buckle - axial buckling for 2D problems
!>
  SUBROUTINE read_files(prob_dim,core_x_size,core_y_size,assm_pitch,prob_sym,assm_map,h_x,h_y, &
                        bc_opt,albedos,num_eg,nsplit,tol_xkeff,tol_xflux,tol_max_iter,nodal_method, &
                        num_assm_reg,assm_xs,refl_mat,ax_buckle)
    !input/output variables
    INTEGER(ki4), INTENT(OUT) :: prob_dim,core_x_size,core_y_size,num_eg,nsplit,tol_max_iter
    INTEGER(ki4), INTENT(OUT) :: num_assm_reg,refl_mat
    REAL(kr8), INTENT(OUT) :: assm_pitch,tol_xkeff,tol_xflux,ax_buckle
    CHARACTER(100), INTENT(OUT) :: prob_sym,bc_opt,nodal_method
    INTEGER(ki4), INTENT(OUT), ALLOCATABLE :: assm_map(:,:)
    REAL(kr8), INTENT(OUT), ALLOCATABLE :: h_x(:),h_y(:),albedos(:)
    TYPE(macro_assm_xs_type), INTENT(OUT), ALLOCATABLE :: assm_xs(:)
    !local variables
    INTEGER(ki4) :: t_int,i
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

    !open xs input unit
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

    !assign variables to found or default values
    prob_dim=p_d
    core_x_size=x_s
    core_y_size=y_s
    assm_pitch=a_p
    prob_sym=p_s
    ALLOCATE(assm_map(x_s,y_s))
    assm_map=a_m
    ALLOCATE(h_x(x_s),h_y(y_s))
    h_x=d_x
    h_y=d_y
    bc_opt=b_o
    IF(b_o .EQ. 'albedo')THEN
      ALLOCATE(albedos(n_g))
      albedos=alb
    ENDIF
    num_eg=n_g
    nsplit=ns
    tol_xkeff=t_xk
    tol_xflux=t_xf
    tol_max_iter=t_m_i
    nodal_method=n_m
    num_assm_reg=n_a_r
    ALLOCATE(assm_xs(n_a_r))
    refl_mat=r_m
    ax_buckle=a_b
    DO i=1,n_a_r
      ALLOCATE(assm_xs(i)%D(n_g))
      ALLOCATE(assm_xs(i)%chi(n_g))
      ALLOCATE(assm_xs(i)%nusigma_f(n_g))
      ALLOCATE(assm_xs(i)%sigma_t(n_g))
      ALLOCATE(assm_xs(i)%sigma_a(n_g))
      ALLOCATE(assm_xs(i)%nu(n_g))
      ALLOCATE(assm_xs(i)%sigma_r(n_g))
      ALLOCATE(assm_xs(i)%sigma_scat(n_g,n_g))
      assm_xs(i)%mat_id=a_x(i)%mat_id
      assm_xs(i)%D=a_x(i)%D
      assm_xs(i)%chi=a_x(i)%chi
      assm_xs(i)%nusigma_f=a_x(i)%nusigma_f
      assm_xs(i)%nu=a_x(i)%nu
      assm_xs(i)%sigma_a=a_x(i)%sigma_a
      assm_xs(i)%sigma_t=a_x(i)%sigma_t
      assm_xs(i)%sigma_r=a_x(i)%sigma_r
      assm_xs(i)%sigma_scat=a_x(i)%sigma_scat
    ENDDO
  ENDSUBROUTINE read_files

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine reads in the base input file Version 1
!>
  SUBROUTINE read_base_inp_v1()
    INTEGER(ki4) :: i,j,t_int,nwords
    CHARACTER(ll_max) :: t_char,words(lp_max)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !data for caseid block
    i=1
    blocks(i)%bname='[CASE_DETAILS]'
    blocks(i)%num_cards=7
    ALLOCATE(blocks(i)%cards(blocks(i)%num_cards))
    j=1
    blocks(i)%cards(j)%cname='title'
    blocks(i)%cards(j)%getcard=>get_title
    j=2
    blocks(i)%cards(j)%cname='nsplit'
    blocks(i)%cards(j)%getcard=>get_nsplit
    j=3
    blocks(i)%cards(j)%cname='k_eps'
    blocks(i)%cards(j)%getcard=>get_k_eps
    j=4
    blocks(i)%cards(j)%cname='phi_eps'
    blocks(i)%cards(j)%getcard=>get_phi_eps
    j=5
    blocks(i)%cards(j)%cname='max_its'
    blocks(i)%cards(j)%getcard=>get_max_its
    j=6
    blocks(i)%cards(j)%cname='nodal_method'
    blocks(i)%cards(j)%getcard=>get_nodal_method
    j=7
    blocks(i)%cards(j)%cname='weilandt'
    blocks(i)%cards(j)%getcard=>get_weilandt

    !data for CORE block
    i=2
    blocks(i)%bname='[CORE]'
    blocks(i)%num_cards=8
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
    j=6
    blocks(i)%cards(j)%cname='bc'
    blocks(i)%cards(j)%getcard=>get_bc
    j=7
    blocks(i)%cards(j)%cname='refl_mat'
    blocks(i)%cards(j)%getcard=>get_refl_mat
    j=8
    blocks(i)%cards(j)%cname='buckling'
    blocks(i)%cards(j)%getcard=>get_buckling

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

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine reads in the cross sections Version 1
!>
  SUBROUTINE read_xs_inp_v1()
    CHARACTER(ll_max) :: t_char,words(lp_max)
    INTEGER(ki4) :: ios,nwords,num_local_xs,i,j

    !get info on the number of materials and number of energy groups
    REWIND(in_unit)
    CALL get_next_line(t_char,ios)
    CALL parse(t_char,' ',words,nwords)
    READ(words(2),*)num_local_xs
    READ(words(3),*)n_g

    DO i=1,n_a_r
      ALLOCATE(a_x(i)%D(n_g))
      ALLOCATE(a_x(i)%chi(n_g))
      ALLOCATE(a_x(i)%nusigma_f(n_g))
      ALLOCATE(a_x(i)%nu(n_g))
      ALLOCATE(a_x(i)%sigma_a(n_g))
      ALLOCATE(a_x(i)%sigma_t(n_g))
      ALLOCATE(a_x(i)%sigma_r(n_g))
      ALLOCATE(a_x(i)%sigma_scat(n_g,n_g))
      a_x(i)%D=0.0
      a_x(i)%chi=0.0
      a_x(i)%nusigma_f=0.0
      a_x(i)%nu=0.0
      a_x(i)%sigma_a=0.0
      a_x(i)%sigma_t=0.0
      a_x(i)%sigma_r=0.0
      a_x(i)%sigma_scat=0.0
    ENDDO

    !get xs data for each of the assembly regions
    DO i=1,n_a_r
      REWIND(in_unit)
      DO
        CALL get_next_line(t_char,ios)
        IF(ios .NE. 0)CALL fatal_error('Could not find xs for material "'//TRIM(a_x(i)%mat_id)//'"')
        CALL parse(t_char,' ',words,nwords)
        words(1)=lowercase(words(1))
        !check to see if it's a start of a new xs
        IF(words(1) .EQ. 'id')THEN
          !check to see if it's the correct xs
          IF(words(2) .EQ. a_x(i)%mat_id)THEN
            !get d
            CALL get_next_line(t_char,ios)
            BACKSPACE(in_unit)
            READ(in_unit,*)a_x(i)%D(:)
            !get chi
            CALL get_next_line(t_char,ios)
            BACKSPACE(in_unit)
            READ(in_unit,*)a_x(i)%chi(:)
            !get numsigma_f
            CALL get_next_line(t_char,ios)
            BACKSPACE(in_unit)
            READ(in_unit,*)a_x(i)%nusigma_f(:)
            !get nu
            CALL get_next_line(t_char,ios)
            BACKSPACE(in_unit)
            READ(in_unit,*)a_x(i)%nu(:)
            !get sigma_a
            CALL get_next_line(t_char,ios)
            BACKSPACE(in_unit)
            READ(in_unit,*)a_x(i)%sigma_a(:)
            !get sigma_scat
            DO j=1,n_g
              CALL get_next_line(t_char,ios)
              BACKSPACE(in_unit)
              READ(in_unit,*)a_x(i)%sigma_scat(j,:)
            ENDDO
            EXIT
          ENDIF
        ENDIF
      ENDDO
      DO j=1,n_g
        a_x(i)%sigma_t(j)=a_x(i)%sigma_a(j)+SUM(a_x(i)%sigma_scat(:,j))
      ENDDO
    ENDDO

    call edit_xs(a_x,n_a_r,n_g)
  ENDSUBROUTINE read_xs_inp_v1

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine finds a given block in the input file
!> @param block_name - name of the block to find
!>
  SUBROUTINE find_block(block_name)
    CHARACTER(MAX_CARDNAME_LEN), INTENT(IN) :: block_name

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

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine gets the next line in a file, ignoring blanks and commented lines
!> @param line - character string of the line that's been retrieved
!> @param ios - read error indicator
!>
  SUBROUTINE get_next_line(line,ios)
    CHARACTER(ll_max), INTENT(OUT) :: line
    INTEGER(ki4), INTENT(OUT) :: ios

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

!---------------------------------------------------------------------------------------------------
!> @brief Subroutine to read in the problem title
!>
  SUBROUTINE get_title(this_card,wwords)
    CLASS(cardType), INTENT(INOUT) :: this_card
    CHARACTER(ll_max), INTENT(INOUT) :: wwords(:)

    INTEGER(ki4) :: len_str,i

    CALL print_log(TRIM(this_card%cname)//' card found')

    len_str=SIZE(wwords)

    prob_title=''
    DO i=2,len_str
      prob_title=TRIM(prob_title)//' '//TRIM(wwords(i))
    ENDDO
    prob_title=TRIM(ADJUSTL(prob_title))
  ENDSUBROUTINE get_title

!---------------------------------------------------------------------------------------------------
!> @brief Subroutine to read in the problem dimensions
!>
  SUBROUTINE get_dim(this_card,wwords)
    CLASS(cardType), INTENT(INOUT) :: this_card
    CHARACTER(ll_max), INTENT(INOUT) :: wwords(:)

    CALL print_log(TRIM(this_card%cname)//' card found')

    IF(wwords(2) .EQ. '2D')THEN
      p_d=2
    ELSE
      CALL fatal_error(TRIM(wwords(2))//' is not a valid dimension.')
    ENDIF
  ENDSUBROUTINE get_dim

!---------------------------------------------------------------------------------------------------
!> @brief Subroutine to read in the problem size
!>
  SUBROUTINE get_size(this_card,wwords)
    CLASS(cardType), INTENT(INOUT) :: this_card
    CHARACTER(ll_max), INTENT(INOUT) :: wwords(:)

    CALL print_log(TRIM(this_card%cname)//' card found')

    READ(wwords(2),*)c_s
  ENDSUBROUTINE get_size

!---------------------------------------------------------------------------------------------------
!> @brief Subroutine to read in the problem pitch
!>
  SUBROUTINE get_apitch(this_card,wwords)
    CLASS(cardType), INTENT(INOUT) :: this_card
    CHARACTER(ll_max), INTENT(INOUT) :: wwords(:)

    CALL print_log(TRIM(this_card%cname)//' card found')

    READ(wwords(2),*)a_p
  ENDSUBROUTINE get_apitch

!---------------------------------------------------------------------------------------------------
!> @brief Subroutine to read in the core symmetry
!>
  SUBROUTINE get_sym(this_card,wwords)
    CLASS(cardType), INTENT(INOUT) :: this_card
    CHARACTER(ll_max), INTENT(INOUT) :: wwords(:)

    CALL print_log(TRIM(this_card%cname)//' card found')

    p_s=TRIM(lowercase(wwords(2)))
  ENDSUBROUTINE get_sym

!---------------------------------------------------------------------------------------------------
!> @brief Subroutine to read in the core assembly map
!>
  SUBROUTINE get_assm_map(this_card,wwords)
    CLASS(cardType), INTENT(INOUT) :: this_card
    CHARACTER(ll_max), INTENT(INOUT) :: wwords(:)

    INTEGER(ki4) :: nwords
    CHARACTER(ll_max) :: t_char,words(lp_max)
    INTEGER(ki4) :: i,ios,j,oct_sym

    LOGICAL :: odd_prob

    wwords(1)=TRIM(wwords(1))
    CALL print_log(TRIM(this_card%cname)//' card found')

    IF(MOD(c_s,2) .EQ. 1)odd_prob=.TRUE.
    !allocate the assembly map based upon problem size and symmetry
    !this is the actual problem we will solve, and remember again that the core is assumed square
    ! TODO need to implement ragged core
    SELECTCASE(p_s)
      CASE('full')
        x_s=c_s
        y_s=c_s
      CASE('half')
        x_s=CEILING((c_s-1.0D-2)/2.0)
        y_s=c_s
      CASE('qtr')
        x_s=CEILING((c_s-1.0D-2)/2.0)
        y_s=CEILING((c_s-1.0D-2)/2.0)
      CASE DEFAULT
        CALL fatal_error(TRIM(p_s)//' is not a valid symmetry')
    ENDSELECT
    ALLOCATE(a_m(x_s,y_s))
    ALLOCATE(d_x(x_s),d_y(y_s))
    d_x=a_p
    d_y=a_p
    SELECTCASE(p_s)
      CASE('half')
        d_x(1)=a_p*0.5D0
      CASE('qtr')
        d_x(1)=a_p*0.5D0
        d_y(1)=a_p*0.5D0
    ENDSELECT
    a_m=0

    !read in the core map
    DO i=1,y_s
      CALL get_next_line(t_char,ios)
      IF(ios .NE. 0)CALL fatal_error('Given assembly map incomplete!')
      !parse to get all data on line
      CALL parse(t_char,' ',words,nwords)
      DO j=1,nwords
        READ(words(j),*)a_m(j,i)
      ENDDO
    ENDDO

    !correct if given an octant symmetry
    IF(p_s .EQ. 'qtr')THEN
      oct_sym=0
      DO i=1,y_s
        DO j=i+1,x_s
          oct_sym=oct_sym+ABS(a_m(j,i))
        ENDDO
      ENDDO
      IF(oct_sym .EQ. 0)THEN
        DO i=1,y_s
          DO j=i+1,x_s
            a_m(j,i)=a_m(i,j)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDSUBROUTINE get_assm_map

  !---------------------------------------------------------------------------------------------------
  !> @brief Subroutine to read in the boundary conditions option
  !>
    SUBROUTINE get_bc(this_card,wwords)
      CLASS(cardType), INTENT(INOUT) :: this_card
      CHARACTER(ll_max), INTENT(INOUT) :: wwords(:)
      INTEGER(ki4) :: i
      REAL(kr8), ALLOCATABLE :: hxt(:),hyt(:)
      INTEGER(ki4), ALLOCATABLE :: amt(:,:)

      CALL print_log(TRIM(this_card%cname)//' card found')

      b_o=TRIM(ADJUSTL(wwords(2)))

      SELECTCASE(b_o)
        CASE('vac','vacuum','reflective','zero') !nothing to do, supported
        CASE('reflector')
          !add the reflector on the boundary
          ALLOCATE(amt(x_s,y_s))
          ALLOCATE(hxt(x_s),hyt(y_s))
          amt=a_m
          hxt=d_x
          hyt=d_y
          DEALLOCATE(a_m,d_x,d_y)
          SELECT CASE(p_s)
            CASE('full')
              x_s=x_s+2
              y_s=y_s+2
            CASE('half')
              x_s=x_s+1
              y_s=y_s+2
            CASE('qtr')
              x_s=x_s+1
              y_s=y_s+1
          ENDSELECT
          ALLOCATE(a_m(x_s,y_s))
          ALLOCATE(d_x(x_s),d_y(y_s))
          a_m=0
          d_x=a_p
          d_y=a_p
          SELECT CASE(p_s)
            CASE('full')
              a_m(2:x_s-1,2:y_s-1)=amt
              d_x(2:x_s-1)=hxt
              d_y(2:y_s-1)=hyt
            CASE('half')
              a_m(1:x_s-1,2:y_s-1)=amt
              d_x(1:x_s-1)=hxt
              d_y(2:y_s-1)=hyt
            CASE('qtr')
              a_m(1:x_s-1,1:y_s-1)=amt
              d_x(1:x_s-1)=hxt
              d_y(1:y_s-1)=hyt
          ENDSELECT
        CASE('albedo')
          DO i=3,1000000
            IF(wwords(i) .EQ. '')EXIT
          ENDDO
          n_g=i-3
          ALLOCATE(alb(n_g))
          DO i=3,1000000
            IF(wwords(i) .EQ. '')EXIT
            READ(wwords(i),*)alb(i-2)
          ENDDO
        CASE DEFAULT
          CALL fatal_error('Invalid boundary condition given: '//TRIM(b_o))
      ENDSELECT
    ENDSUBROUTINE get_bc

!---------------------------------------------------------------------------------------------------
!> @brief Subroutine to read in the nsplit option
!>
  SUBROUTINE get_nsplit(this_card,wwords)
    CLASS(cardType), INTENT(INOUT) :: this_card
    CHARACTER(ll_max), INTENT(INOUT) :: wwords(:)

    INTEGER(ki4) :: nwords
    CHARACTER(ll_max) :: t_char,words(lp_max)
    INTEGER(ki4) :: i,ios,j,oct_sym

    wwords(1)=TRIM(wwords(1))
    CALL print_log(TRIM(this_card%cname)//' card found')

    READ(wwords(2),*)ns

  ENDSUBROUTINE get_nsplit

!---------------------------------------------------------------------------------------------------
!> @brief Subroutine to read in the keff tolerance
!>
  SUBROUTINE get_k_eps(this_card,wwords)
    CLASS(cardType), INTENT(INOUT) :: this_card
    CHARACTER(ll_max), INTENT(INOUT) :: wwords(:)

    INTEGER(ki4) :: nwords
    CHARACTER(ll_max) :: t_char,words(lp_max)
    INTEGER(ki4) :: i,ios,j,oct_sym

    wwords(1)=TRIM(wwords(1))
    CALL print_log(TRIM(this_card%cname)//' card found')

    READ(wwords(2),*)t_xk

  ENDSUBROUTINE get_k_eps

!---------------------------------------------------------------------------------------------------
!> @brief Subroutine to read in the flux tolerance
!>
  SUBROUTINE get_phi_eps(this_card,wwords)
    CLASS(cardType), INTENT(INOUT) :: this_card
    CHARACTER(ll_max), INTENT(INOUT) :: wwords(:)

    INTEGER(ki4) :: nwords
    CHARACTER(ll_max) :: t_char,words(lp_max)
    INTEGER(ki4) :: i,ios,j,oct_sym

    wwords(1)=TRIM(wwords(1))
    CALL print_log(TRIM(this_card%cname)//' card found')

    READ(wwords(2),*)t_xf

  ENDSUBROUTINE get_phi_eps

!---------------------------------------------------------------------------------------------------
!> @brief Subroutine to read in max number of iterations
!>
  SUBROUTINE get_max_its(this_card,wwords)
    CLASS(cardType), INTENT(INOUT) :: this_card
    CHARACTER(ll_max), INTENT(INOUT) :: wwords(:)

    INTEGER(ki4) :: nwords
    CHARACTER(ll_max) :: t_char,words(lp_max)
    INTEGER(ki4) :: i,ios,j,oct_sym

    wwords(1)=TRIM(wwords(1))
    CALL print_log(TRIM(this_card%cname)//' card found')

    READ(wwords(2),*)t_m_i

  ENDSUBROUTINE get_max_its

!---------------------------------------------------------------------------------------------------
!> @brief Subroutine to read in which nodal method is being used, right now only fd and poly supported
!>
  SUBROUTINE get_nodal_method(this_card,wwords)
    CLASS(cardType), INTENT(INOUT) :: this_card
    CHARACTER(ll_max), INTENT(INOUT) :: wwords(:)

    INTEGER(ki4) :: nwords
    CHARACTER(ll_max) :: t_char,words(lp_max)
    INTEGER(ki4) :: i,ios,j,oct_sym

    wwords(1)=TRIM(wwords(1))
    CALL print_log(TRIM(this_card%cname)//' card found')

    n_m=TRIM(ADJUSTL(wwords(2)))

    IF(n_m .NE. 'poly' .AND. n_m .NE. 'fd')THEN
      CALL fatal_error('Invalid nodal method: '//TRIM(n_m))
    ENDIF

  ENDSUBROUTINE get_nodal_method

!---------------------------------------------------------------------------------------------------
!> @brief Subroutine to read in which nodal method is being used, right now only fd and poly supported
!>
  SUBROUTINE get_weilandt(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(:)

    INTEGER(ki4) :: nwords
    CHARACTER(ll_max) :: t_char,words(lp_max)
    INTEGER(ki4) :: i,ios,j,oct_sym

    CALL print_log(TRIM(this_card%cname)//' card found')

    wwords(2)=TRIM(ADJUSTL(wwords(2)))
    READ(wwords(2),*)dl_weilandt

    IF(dl_weilandt .LT. 0.0D0)THEN
      CALL fatal_error('Invalid Weilandt shift value: '//TRIM(wwords(2)))
    ENDIF

  ENDSUBROUTINE get_weilandt

!---------------------------------------------------------------------------------------------------
!> @brief Subroutine to read in the cross section filename
!>
  SUBROUTINE get_xs_file(this_card,wwords)
    CLASS(cardType), INTENT(INOUT) :: this_card
    CHARACTER(ll_max), INTENT(INOUT) :: wwords(:)

    CALL print_log(TRIM(this_card%cname)//' card found')

    xs_in=wwords(2)
  ENDSUBROUTINE get_xs_file

!---------------------------------------------------------------------------------------------------
!> @brief Subroutine to read in the cross section mapping
!>
  SUBROUTINE get_xs_map(this_card,wwords)
    CLASS(cardType), INTENT(INOUT) :: this_card
    CHARACTER(ll_max), INTENT(INOUT) :: wwords(:)

    INTEGER(ki4) :: i,ios,nwords
    CHARACTER(ll_max) :: t_char,words(lp_max)

    CALL print_log(TRIM(this_card%cname)//' card found')

    READ(wwords(2),*)n_a_r
    IF(MAXVAL(a_m) .GT. n_a_r)CALL fatal_error('Assembly map index outside of range of xs.')
    ALLOCATE(a_x(n_a_r))

    DO i=1,n_a_r
      CALL get_next_line(t_char,ios)
      IF(ios .NE. 0)CALL fatal_error('Given cross section map incomplete!')
      CALL parse(t_char,' ',words,nwords)
      IF(nwords .LT. 3)CALL fatal_error('Each xs map data should have 3 entries!')
      READ(words(1),*)ios
      IF(ios .NE. i)CALL fatal_error('xs map should be in order and complete!')
      IF(words(2) .EQ. 'macro')THEN
        a_x(i)%mat_id=TRIM(ADJUSTL(words(3)))
      ELSE
        CALL fatal_error(TRIM(words(2))//' xs not yet supported.')
      ENDIF
    ENDDO
  ENDSUBROUTINE get_xs_map

!---------------------------------------------------------------------------------------------------
!> @brief Subroutine to read in reflector material to fill in gaps, and/or to compute a reflector
!>    albedo boundary condition
!>
  SUBROUTINE get_refl_mat(this_card,wwords)
    CLASS(cardType), INTENT(INOUT) :: this_card
    CHARACTER(ll_max), INTENT(INOUT) :: wwords(:)

    CALL print_log(TRIM(this_card%cname)//' card found')

    READ(wwords(2),*)r_m
  ENDSUBROUTINE get_refl_mat

!---------------------------------------------------------------------------------------------------
!> @brief Subroutine to read in axial buckling
!>
  SUBROUTINE get_buckling(this_card,wwords)
    CLASS(cardType), INTENT(INOUT) :: this_card
    CHARACTER(ll_max), INTENT(INOUT) :: wwords(:)

    CALL print_log(TRIM(this_card%cname)//' card found')
    IF(wwords(2) .EQ. 'height')THEN
      READ(wwords(3),*)a_b
      a_b=pi**2/a_b**2
    ELSE
      READ(wwords(2),*)a_b
    ENDIF
  ENDSUBROUTINE get_buckling
ENDMODULE input_module
