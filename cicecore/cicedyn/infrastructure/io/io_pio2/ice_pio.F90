!============================================================================
!  Writes netcdf files
!    Created by Mariana Vertenstein, June 2009

  module ice_pio

  use ice_kinds_mod
  use ice_blocks
  use ice_broadcast
  use ice_domain_size
  use ice_communicate, only : my_task, master_task, get_num_procs, MPI_COMM_ICE
  use ice_domain, only : nblocks, blocks_ice
  use ice_fileunits, only : nu_diag
  use ice_exit, only: abort_ice
  use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
  use pio

  implicit none
  private

  interface ice_pio_initdecomp
     module procedure ice_pio_initdecomp_2d
     module procedure ice_pio_initdecomp_3d
     module procedure ice_pio_initdecomp_4d
     module procedure ice_pio_initdecomp_3d_inner
  end interface

  public ice_pio_init
  public ice_pio_initdecomp

#ifdef CESMCOUPLED
  type(iosystem_desc_t), pointer :: ice_pio_subsystem
#else
  type(iosystem_desc_t)          :: ice_pio_subsystem
#endif

!===============================================================================

  contains

!================================================================================
!    Initialize the io subsystem
!    2009-Feb-17 - J. Edwards - initial version

   subroutine ice_pio_init(pio_options, mode, filename, File, clobber, cdf64)

#ifdef CESMCOUPLED
   use shr_pio_mod, only: shr_pio_getiosys, shr_pio_getiotype
#else
#ifdef GPTL
   use perf_mod, only : t_initf
#endif
#endif

   implicit none
   character(len=*)  , intent(in)              :: pio_options(:)  ! pio namelist options
   character(len=*)  , intent(in),    optional :: mode
   character(len=*)  , intent(in),    optional :: filename
   type(file_desc_t) , intent(inout), optional :: File
   logical           , intent(in),    optional :: clobber
   logical           , intent(in),    optional :: cdf64
   ! local variables

   integer (int_kind) :: &
      nml_error          ! namelist read error flag

   integer :: nprocs
   integer :: iostride
   integer :: basetask
   integer :: numiotasks
   integer :: rearranger
   integer :: iotype
   logical :: exists
   logical :: lclobber
   logical :: lcdf64
   integer :: status
   integer :: nmode
   character(len=*), parameter :: subname = '(ice_pio_init)'

#ifdef CESMCOUPLED
   ice_pio_subsystem => shr_pio_getiosys(inst_name)
   iotype =  shr_pio_getiotype(inst_name)
#else

#ifdef GPTL
   !--- initialize gptl
   call t_initf('undefined_NLFileName', LogPrint=.false., mpicom=MPI_COMM_ICE, &
         MasterTask=.true.)
#endif
   nprocs = get_num_procs()
   basetask = min(1,nprocs-1)

   if ((trim(pio_options(1)) == '-99') .or. (trim(pio_options(1)) == 'netcdf')) then
      iotype = PIO_IOTYPE_NETCDF
   else if (trim(pio_options(1)) == 'pnetcdf') then
      iotype = PIO_IOTYPE_PNETCDF
   else if (trim(pio_options(1)) == 'netcdf4c') then
      iotype = PIO_IOTYPE_NETCDF4C
   else if (trim(pio_options(1)) == 'netcdf4p') then
      iotype = PIO_IOTYPE_NETCDF4P
   else
      if (my_task == master_task) then
         write(nu_diag,'(a)') ' no valid iotype set'
      end if
      call abort_ice(subname//'ERROR: aborting with no valid iotype')
      return
   end if

    if ((trim(pio_options(2)) == '-99') .or. (trim(pio_options(2)) == 'box')) then
       rearranger = PIO_REARR_BOX
    else if (trim(pio_options(2)) == 'subset') then
       rearranger = PIO_REARR_SUBSET
    else
       if (my_task == master_task) then
          write(nu_diag,'(a)') ' no valid pio_rearranger set'
       end if
       call abort_ice(subname//'ERROR: aborting with no valid pio_rearranger')
       return
    end if

   if (trim(pio_options(3)) == '-99') then
      iostride = -99
   else
      read(pio_options(3),*)iostride
   end if

   if (trim(pio_options(4)) == '-99') then
      numiotasks = -99
   else
      read(pio_options(4),*)numiotasks
   end if

   ! check for parallel IO, it requires at least two io pes
   if (nprocs > 1 .and. numiotasks == 1 .and. &
        (iotype .eq. PIO_IOTYPE_PNETCDF .or. iotype .eq. PIO_IOTYPE_NETCDF4P)) then
      numiotasks = 2
      iostride = min(iostride, nprocs/2)
      if (my_task == master_task) then
         write(nu_diag,*) ' parallel io requires at least two io pes - following parameters are updated:'
         write(nu_diag,*) trim(subname), ' : iostride = ', iostride
         write(nu_diag,*) trim(subname), ' : numiotasks = ', numiotasks
      end if
   endif

   ! check/set/correct io pio parameters
   if (iostride > 0 .and. numiotasks < 0) then
      numiotasks = max(1, nprocs/iostride)
      if (my_task == master_task ) write(nu_diag,*) trim(subname), ' : update numiotasks = ', numiotasks
    else if(numiotasks > 0 .and. iostride < 0) then
       iostride = max(1, nprocs/numiotasks)
       if (my_task == master_task) write(nu_diag,*) trim(subname), ' : update iostride = ', iostride
    else if(numiotasks < 0 .and. iostride < 0) then
       iostride = max(1,nprocs/4)
       numiotasks = max(1,nprocs/iostride)
       if (my_task == master_task) write(nu_diag,*) trim(subname), ' : update numiotasks = ', numiotasks
       if (my_task == master_task) write(nu_diag,*) trim(subname), ' : update iostride = ', iostride
    end if

    if (basetask + (iostride)*(numiotasks-1) >= nprocs ) then
       if (nprocs < 100) then
          iostride = max(1, nprocs/4)
       else if(nprocs < 1000) then
          iostride = max(1, nprocs/8)
       else
          iostride = max(1, nprocs/16)
       end if
       if(iostride > 1) then
          numiotasks = nprocs/iostride
          basetask = min(1, nprocs-1)
       else
          numiotasks = nprocs
          basetask = 0
       end if
       if (my_task == master_task) then
          write(nu_diag,*) 'iostride, iotasks or root out of bounds - resetting to defaults:'
          write(nu_diag,*) trim(subname), ' : basetask = ', basetask
          write(nu_diag,*) trim(subname), ' : iostride = ', iostride
          write(nu_diag,*) trim(subname), ' : numiotasks = ', numiotasks
       end if
    end if

    if (my_task == master_task) then
       write(nu_diag,'(a,a,i6)') subname,' nprocs     = ',nprocs
       write(nu_diag,'(a,a,i6)') subname,' iostride   = ',iostride
       write(nu_diag,'(a,a,i6)') subname,' root       = ',basetask
       write(nu_diag,'(a,a,i6)') subname,' numiotasks = ',numiotasks
       write(nu_diag,'(a,a,i6,a)') subname,' iotype = ',iotype,' [PNETCDF| NETCDF | NETCDF4C | NETCDF4P]'
       write(nu_diag,'(a,a,i6,a)') subname,' rearranger = ',rearranger,' [BOX | SUBSET]'
    end if
    ! set PIO debug level
    call pio_setdebuglevel(6)

    call pio_init(my_task, MPI_COMM_ICE, numiotasks, master_task, iostride, rearranger, &
                  ice_pio_subsystem, base=basetask)

   ! TODO
   !--- initialize rearranger options
   ! rearranger defaults
   ! pio_rearr_comm_type = PIO_REARR_COMM_P2P
   ! pio_rearr_comm_fcd = PIO_REARR_COMM_FC_2D_ENABLE
   ! pio_rearr_comm_enable_hs_comp2io = .true.
   ! pio_rearr_comm_enable_isend_comp2io = .false.
   ! pio_rearr_comm_max_pend_req_comp2io = 0
   ! pio_rearr_comm_enable_hs_io2comp = .false.
   ! pio_rearr_comm_enable_isend_io2comp = .true.
   ! pio_rearr_comm_max_pend_req_io2comp = 64
   ! ! set PIO rearranger options
   ! if (my_task == master_task) write(nu_diag,*) subname// ' calling pio_set_rearr_opts'
   ! ret = pio_set_rearr_opts(ice_pio_subsystem, pio_rearr_comm_type, &
   !      pio_rearr_comm_fcd, &
   !      pio_rearr_comm_enable_hs_comp2io, &
   !      pio_rearr_comm_enable_isend_comp2io, &
   !      pio_rearr_comm_max_pend_req_comp2io, &
   !      pio_rearr_comm_enable_hs_io2comp, &
   !      pio_rearr_comm_enable_isend_io2comp, &
   !      pio_rearr_comm_max_pend_req_io2comp)
   ! if(ret /= PIO_NOERR) then
   !    call abort_ice(subname//'ERROR: aborting in pio_set_rearr_opts')
   ! end if
#endif
   if (present(mode) .and. present(filename) .and. present(File)) then

      if (trim(mode) == 'write') then
         lclobber = .false.
         if (present(clobber)) lclobber=clobber

         lcdf64 = .false.
         if (present(cdf64)) lcdf64=cdf64

         if (File%fh<0) then
            ! filename not open
            inquire(file=trim(filename),exist=exists)
            if (exists) then
               if (lclobber) then
                  nmode = pio_clobber
                  if (lcdf64) nmode = ior(nmode,PIO_64BIT_OFFSET)
                  status = pio_createfile(ice_pio_subsystem, File, iotype, trim(filename), nmode)
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,' create file ',trim(filename)
                  end if
               else
                  nmode = pio_write
                  status = pio_openfile(ice_pio_subsystem, File, iotype, trim(filename), nmode)
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,' open file ',trim(filename)
                  end if
               endif
            else
               nmode = pio_noclobber
               if (lcdf64) nmode = ior(nmode,PIO_64BIT_OFFSET)
               status = pio_createfile(ice_pio_subsystem, File, iotype, trim(filename), nmode)
               if (my_task == master_task) then
                  write(nu_diag,*) subname,' create file ',trim(filename)
               end if
            endif
         else
            ! filename is already open, just return
         endif
      end if

      if (trim(mode) == 'read') then
         inquire(file=trim(filename),exist=exists)
         if (exists) then
            status = pio_openfile(ice_pio_subsystem, File, iotype, trim(filename), pio_nowrite)
         else
            if(my_task==master_task) then
               write(nu_diag,*) 'ice_pio_ropen ERROR: file invalid ',trim(filename)
            end if
            call abort_ice(subname//'ERROR: aborting with invalid file')
         endif
      end if

   end if

   end subroutine ice_pio_init

!================================================================================

   subroutine ice_pio_initdecomp_2d(iodesc, precision)

      type(io_desc_t), intent(out) :: iodesc
      integer(kind=int_kind), optional, intent(in) :: precision

      integer (kind=int_kind) :: &
          iblk,ilo,ihi,jlo,jhi,lon,lat,i,j,n,k

      type(block) :: this_block

      integer(kind=int_kind), pointer :: dof2d(:)
      integer(kind=int_kind) :: lprecision
      character(len=*), parameter :: subname = '(ice_pio_initdecomp_2d)'

      lprecision = 8
      if (present(precision)) lprecision = precision

      allocate(dof2d(nx_block*ny_block*nblocks))

      n=0
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j=1,ny_block
         do i=1,nx_block
            n = n+1
            if (j < jlo .or. j>jhi) then
               dof2d(n) = 0
            else if (i < ilo .or. i > ihi) then
               dof2d(n) = 0
            else
               lon = this_block%i_glob(i)
               lat = this_block%j_glob(j)
               dof2d(n) = (lat-1)*nx_global + lon
            endif
         enddo !i
         enddo !j
      end do

      if (lprecision == 8) then
         call pio_initdecomp(ice_pio_subsystem, pio_double, (/nx_global,ny_global/), &
              dof2d, iodesc)
      else
         call pio_initdecomp(ice_pio_subsystem, pio_real, (/nx_global,ny_global/), &
              dof2d, iodesc)
      endif

      deallocate(dof2d)

   end subroutine ice_pio_initdecomp_2d

!================================================================================

   subroutine ice_pio_initdecomp_3d (ndim3, iodesc, remap, precision)

      integer(kind=int_kind), intent(in) :: ndim3
      type(io_desc_t), intent(out) :: iodesc
      logical, optional :: remap
      integer(kind=int_kind), optional, intent(in) :: precision
      integer (kind=int_kind) :: &
          iblk,ilo,ihi,jlo,jhi,lon,lat,i,j,n,k

      type(block) :: this_block
      logical :: lremap
      integer(kind=int_kind), pointer :: dof3d(:)
      integer(kind=int_kind) :: lprecision
      character(len=*), parameter :: subname = '(ice_pio_initdecomp_3d)'

      lprecision = 8
      if (present(precision)) lprecision = precision

      allocate(dof3d(nx_block*ny_block*nblocks*ndim3))
      lremap=.false.
      if (present(remap)) lremap=remap
      if (lremap) then
         ! Reorder the ndim3 and nblocks loops to avoid a temporary array in restart read/write
         n=0
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            do k=1,ndim3
               do j=1,ny_block
                  do i=1,nx_block
                     n = n+1
                     if (j < jlo .or. j>jhi) then
                        dof3d(n)=0
                     else if (i < ilo .or. i > ihi) then
                        dof3d(n) = 0
                     else
                        lon = this_block%i_glob(i)
                        lat = this_block%j_glob(j)
                        dof3d(n) = ((lat-1)*nx_global + lon) + (k-1)*nx_global*ny_global
                     endif
                  enddo !i
               enddo !j
            enddo !ndim3
         enddo ! iblk
   else
         n=0
         do k=1,ndim3
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j=1,ny_block
                  do i=1,nx_block
                     n = n+1
                     if (j < jlo .or. j>jhi) then
                        dof3d(n)=0
                     else if (i < ilo .or. i > ihi) then
                        dof3d(n) = 0
                     else
                        lon = this_block%i_glob(i)
                        lat = this_block%j_glob(j)
                        dof3d(n) = ((lat-1)*nx_global + lon) + (k-1)*nx_global*ny_global
                     endif
                  enddo !i
               enddo !j
            enddo ! iblk
         enddo !ndim3
      endif

      if (lprecision == 8) then
         call pio_initdecomp(ice_pio_subsystem, pio_double, (/nx_global,ny_global,ndim3/), &
              dof3d, iodesc)
      else
         call pio_initdecomp(ice_pio_subsystem, pio_real, (/nx_global,ny_global,ndim3/), &
              dof3d, iodesc)
      endif

      deallocate(dof3d)

   end subroutine ice_pio_initdecomp_3d

!================================================================================

   subroutine ice_pio_initdecomp_3d_inner(ndim3, inner_dim, iodesc, precision)

      integer(kind=int_kind), intent(in) :: ndim3
      logical, intent(in) :: inner_dim
      type(io_desc_t), intent(out) :: iodesc
      integer(kind=int_kind), optional, intent(in) :: precision

      integer (kind=int_kind) :: &
          iblk,ilo,ihi,jlo,jhi,lon,lat,i,j,n,k

      type(block) :: this_block

      integer(kind=int_kind), pointer :: dof3d(:)
      integer(kind=int_kind) :: lprecision
      character(len=*), parameter :: subname = '(ice_pio_initdecomp_3d_inner)'

      lprecision = 8
      if (present(precision)) lprecision = precision

      allocate(dof3d(nx_block*ny_block*nblocks*ndim3))

      n=0
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j=1,ny_block
         do i=1,nx_block
         do k=1,ndim3
            n = n+1
            if (j < jlo .or. j>jhi) then
               dof3d(n) = 0
            else if (i < ilo .or. i > ihi) then
               dof3d(n) = 0
            else
               lon = this_block%i_glob(i)
               lat = this_block%j_glob(j)
               dof3d(n) = k + ((lon-1) + (lat-1)*nx_global)*ndim3
            endif
         end do !ndim3
         enddo  !i
         enddo  !j
      end do    !iblk

      if (lprecision == 8) then
         call pio_initdecomp(ice_pio_subsystem, pio_double, (/ndim3,nx_global,ny_global/), &
              dof3d, iodesc)
      else
         call pio_initdecomp(ice_pio_subsystem, pio_real, (/ndim3,nx_global,ny_global/), &
              dof3d, iodesc)
      endif

      deallocate(dof3d)

   end subroutine ice_pio_initdecomp_3d_inner

!================================================================================

   subroutine ice_pio_initdecomp_4d (ndim3, ndim4, iodesc, precision)

      integer(kind=int_kind), intent(in) :: ndim3, ndim4
      type(io_desc_t), intent(out) :: iodesc
      integer(kind=int_kind), optional, intent(in) :: precision

      integer (kind=int_kind) :: &
          iblk,ilo,ihi,jlo,jhi,lon,lat,i,j,n,k,l

      type(block) :: this_block

      integer(kind=int_kind), pointer :: dof4d(:)
      integer(kind=int_kind) :: lprecision
      character(len=*), parameter :: subname = '(ice_pio_initdecomp_4d)'

      lprecision = 8
      if (present(precision)) lprecision = precision

      allocate(dof4d(nx_block*ny_block*nblocks*ndim3*ndim4))

      n=0
      do l=1,ndim4
      do k=1,ndim3
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j=1,ny_block
         do i=1,nx_block
            n = n+1
            if (j < jlo .or. j>jhi) then
               dof4d(n)=0
            else if (i < ilo .or. i > ihi) then
               dof4d(n) = 0
            else
               lon = this_block%i_glob(i)
               lat = this_block%j_glob(j)
               dof4d(n) = ((lat-1)*nx_global + lon) &
                        + (k-1)*nx_global*ny_global &
                        + (l-1)*nx_global*ny_global*ndim3
            endif
         enddo !i
         enddo !j
      enddo ! iblk
      enddo !ndim3
      enddo !ndim4

      if (lprecision == 8) then
         call pio_initdecomp(ice_pio_subsystem, pio_double, &
             (/nx_global,ny_global,ndim3,ndim4/), dof4d, iodesc)
      else
         call pio_initdecomp(ice_pio_subsystem, pio_real, &
             (/nx_global,ny_global,ndim3,ndim4/), dof4d, iodesc)
      endif

      deallocate(dof4d)

   end subroutine ice_pio_initdecomp_4d

!================================================================================

  end module ice_pio

!================================================================================
