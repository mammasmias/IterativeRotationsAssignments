module routines
  use ira_mod, only: ira_rp
  implicit none
  public

  type t_settings
     integer :: nfiles = 0
     logical :: log_each = .false.
     logical :: log_file = .false.
     character(len=4096), allocatable :: fnames(:)
   contains
     procedure :: parse_args
  end type t_settings


  type t_sym
     integer :: &
          nbas = 0, &
          n_prin_ax = 0

     character(len=10) :: pg
     real(ira_rp), allocatable :: &
          bas_list(:,:,:), &
          ax_out(:,:), &
          angle_out(:), &
          prin_ax(:,:), &
          dHausdorff_out(:)

     character(len=1), allocatable :: op_out(:)

     integer, allocatable :: &
          perm_list(:,:), &
          n_out(:), &
          p_out(:)
  end type t_sym

contains

  subroutine parse_args( self )
    implicit none
    class( t_settings ), intent(inout) :: self
    integer :: narg, nfiles
    integer :: i, u0, ios
    character(len=4096) :: arg, line
    character(len=4096), allocatable :: tmp(:)
    nfiles = 0
    narg = command_argument_count()

    allocate( tmp(1:max(narg,16)) )
    i = 1
    do while( i <= narg )
       call get_command_argument(i, arg)

       select case( arg )
       case( "--log-file" )
          self%log_file = .true.
          i = i + 1

       case( "--log-each" )
          self%log_each = .true.
          i = i + 1

       case( "--list" )
          ! read filenames from a list file, one per line
          if( i + 1 > narg ) then
             write(*,*) "--list requires a filename"
             error stop 1
          end if
          call get_command_argument(i+1, arg)
          open( newunit=u0, file=trim(arg), status="old", action="read", iostat=ios )
          if( ios /= 0 ) then
             write(*,*) "cannot open list file: "//trim(arg)
             error stop 1
          end if
          do
             read( u0, "(a)", iostat=ios ) line
             if( ios /= 0 ) exit
             line = adjustl(line)
             if( len_trim(line) == 0 ) cycle
             if( line(1:1) == "#" ) cycle
             call push( line )
          end do
          close( u0 )
          i = i + 2

       case default
          call push( arg )
          i = i + 1
       end select
    end do
    ! save from tmp to list
    allocate( self%fnames(1:nfiles) )
    self%fnames = tmp(1:nfiles)
    self%nfiles = nfiles
    deallocate(tmp)

  contains

    subroutine push( fname )
      character(len=*), intent(in) :: fname
      character(len=4096), allocatable :: grow(:)
      if( nfiles == size(tmp) ) then
         allocate( grow(1:2*size(tmp)) )
         grow(1:nfiles) = tmp(1:nfiles)
         call move_alloc( grow, tmp )
      end if
      nfiles = nfiles + 1
      tmp(nfiles) = fname
    end subroutine push

  end subroutine parse_args

  function read_xyz( fname, nat, ctyp, pos )result(ios)
    use ira_mod, only: ira_rp
    implicit none
    character(*), intent(in) :: fname
    integer, intent(out) :: nat
    character(len=4), allocatable, intent(out) :: ctyp(:)
    real(ira_rp), allocatable, intent(out) :: pos(:,:)
    integer :: ios
    integer :: u0, i
    open(newunit=u0, file=fname, status="old", iostat=ios)
    if(ios/=0)return
    read(u0,*) nat
    read(u0,*)
    allocate( pos(1:3,1:nat))
    allocate( ctyp(1:nat))
    do i = 1, nat
       read(u0,*) ctyp(i), pos(:,i)
    end do
    close(u0, status="keep")
  end function read_xyz


  subroutine str_to_int( str_arr1, str_arr2, int_arr1, int_arr2 )

    implicit none
    character(len=*), dimension(:), intent(in) :: str_arr1, str_arr2
    integer, dimension(:), intent(out) :: int_arr1, int_arr2

    integer :: i, j, n1, n2, wlen1, wlen2, wlen
    character(:), allocatable :: old(:)
    logical :: isnew
    integer :: nnew, nnew1, nnew2

    if( size(int_arr1) /= size(str_arr1) ) then
       write(*,*) "unequal size arrays1 in str_to_int. Stopping"
       stop
    end if
    if( size(int_arr2) /= size(str_arr2) ) then
       write(*,*) "unequal size arrays2 in str_to_int. Stopping"
       stop
    end if

    wlen1 = len(str_arr1(1))
    wlen2 = len(str_arr2(1))
    wlen = max( wlen1, wlen2 )

    n1 = size(int_arr1, 1)
    n2 = size(int_arr2, 1)
    allocate( character(len=wlen)::old(n1+n2) )

    !! check on str_arr1
    old(:) = "xcba"
    nnew = 1
    old(nnew) = str_arr1(1)
    do i = 2, n1
       !! check if element is new
       isnew = .not. any( old .eq. str_arr1(i) )
       if( isnew ) then
          nnew = nnew + 1
          old( nnew ) = str_arr1(i)
       end if
    end do
    !! keep how many were new in str_arr1
    nnew1 = nnew

    !! check on str_arr2
    do i = 1, n2
       !! check if element is new
       isnew = .not. any( old .eq. str_arr2(i) )
       if( isnew ) then
          nnew = nnew + 1
          old( nnew ) = str_arr2(i)
       end if
    end do
    !! keep nnew2
    nnew2 = nnew


    !! assign int_arr1
    do i = 1, n1
       do j = 1, nnew1
          if( str_arr1(i) .eq. old(j) ) int_arr1(i) = j
       end do
    end do

    !! assign int_arr2
    do i = 1, n2
       do j = 1, nnew2
          if( str_arr2(i) .eq. old(j) ) int_arr2(i) = j
       end do
    end do

    deallocate(old)
  end subroutine str_to_int

  function sofi_compute( nat, typ, coords, sym_thr, ierr, prescreen_ih )result(sym)
    use ira_mod, only: nmax => sofi_nmax
    implicit none
    integer, intent(in) :: nat
    integer, intent(in) :: typ(:)
    real(ira_rp), intent(in) :: coords(:,:)
    real(ira_rp), intent(in) :: sym_thr
    integer, intent(out) :: ierr
    logical, intent(in), optional :: prescreen_ih
    type( t_sym ) :: sym

    logical :: prescr_ih
    real(ira_rp) :: gc(3)
    integer :: i
    real(ira_rp), allocatable :: pos(:,:)

    allocate( sym%bas_list(1:3, 1:3, 1:nmax))
    allocate( sym%perm_list(1:nat, 1:nmax))
    allocate( sym%op_out(1:nmax))
    allocate( sym%n_out(1:nmax))
    allocate( sym%p_out(1:nmax))
    allocate( sym%ax_out(1:3, 1:nmax))
    allocate( sym%angle_out(1:nmax))
    allocate( sym%dHausdorff_out(1:nmax))
    allocate( sym%prin_ax(1:3, 1:nmax))
    prescr_ih = .false.
    if(present(prescreen_ih))prescr_ih = prescreen_ih

    allocate(pos, mold=coords)
    gc = sum(coords,2)/nat
    do i = 1, nat
       pos(:,i) = coords(:,i) - gc
    end do

    call sofi_compute_all( nat, typ, pos, sym_thr, prescr_ih, &
         sym%nbas, sym%bas_list, sym%perm_list, &
         sym%op_out, sym%n_out, sym%p_out, &
         sym%ax_out, sym%angle_out, sym%dHausdorff_out, sym%pg, &
         sym%n_prin_ax, sym%prin_ax, ierr )

  end function sofi_compute

end module routines
