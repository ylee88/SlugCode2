module io

#include "definition.h"

  use grid_data, only : gr_V,                 &
                        gr_xCoord, gr_yCoord, &
                        gr_glb_nx, gr_glb_ny, &
                        gr_nx, gr_ny,         &
                        gr_dx, gr_dy,         &
                        gr_ibeg, gr_iend,     &
                        gr_i0, gr_imax,       &
                        gr_ngc,               &
                        gr_xbeg

  use sim_data, only : sim_name, sim_hdf5, sim_pIO
  use block_data
  use mpi

  implicit none

  integer, save :: nCounter

contains

  subroutine io_writeOutput(init_file_name, t, nstep, ioCounter, eTime, commit_hash)

    implicit none
    character(len=MAX_STRING_LENGTH), intent(IN) :: init_file_name, commit_hash
    real, intent(IN) :: t, eTime
    integer, intent(IN) :: nstep, ioCounter

    if (sim_hdf5) then
      if (sim_pIO) then
        call io_writeHDF5_p(init_file_name, t, nstep, ioCounter, eTime, commit_hash)
      else
        call io_writeHDF5(init_file_name, t, nstep, ioCounter, eTime, commit_hash)
      end if
    else
      call io_writeASCII(ioCounter)
    end if
  end subroutine io_writeOutput

  subroutine io_writeHDF5_p(init_file_name, t, nstep, ioCounter, eTime, commit_hash)

    use HDF5
    use mpi, ONLY: MPI_COMM_WORLD, MPI_INFO_NULL
    use io_hdf5_misc, ONLY: io_writeHDF5_simInfo

    implicit none

    character(len=200) :: ofile
    character(len=5)  :: cCounter
    character(len=50) :: dset_prim, dset_x, dset_y

    character(len=MAX_STRING_LENGTH), intent(IN) :: init_file_name, commit_hash
    real,    intent(IN) :: t, eTime
    integer, intent(IN) :: nstep, ioCounter

    real, allocatable, dimension(:,:,:)   :: img_V
    real, allocatable, dimension(:)       :: yCoord, xCoord

    integer(HSIZE_T), dimension(1) :: dims_XYZ
    integer(HSIZE_T), dimension(2) :: xyoffset
    integer(HSIZE_T), dimension(3) :: dims_V, glb_dims_V, hs_start, hs_stride, hs_count
    integer(HID_T) :: file_id, dspace_id, dset_id, mspace_id
    integer(HID_T) :: plist_file, plist_dset, plist_dxfer

    integer :: i, j, comm, info, error
    integer :: rank_V, rank_XYZ

    allocate(img_V(NSYS_VAR, gr_nx, gr_ny))
    img_V(1:NSYS_VAR, :, :) = gr_V(1:NSYS_VAR, gr_ibeg(XDIM):gr_iend(XDIM), gr_ibeg(YDIM):gr_iend(YDIM))

    write(cCounter,910) ioCounter + 10000
    ofile = trim(sim_name)//'_'//cCounter//'.slug'

    dset_prim = "prim_vars"
    dset_x    = "xCoord"
    dset_y    = "yCoord"

    rank_V     = 3
    dims_V     = (/ NSYS_VAR, gr_nx,     gr_ny     /)
    glb_dims_V = (/ NSYS_VAR, gr_glb_nx, gr_glb_ny /)

    ! get communicator from CAF
    comm = MPI_COMM_WORLD
    info = MPI_INFO_NULL
    ! call MPI_info_create(info, error)

    ! init hdf5 fortran interface
    call H5open_f(error)

    ! set up file access property, plist_file
    call H5Pcreate_f(H5P_FILE_ACCESS_F, plist_file, error)
    call H5Pset_fapl_mpio_f(plist_file, comm, info, error)
    ! Create the file collectively.
    call H5Fcreate_f(ofile, H5F_ACC_TRUNC_F, file_id, error, access_prp=plist_file)
    ! close parallel for fapl
    call H5Pclose_f(plist_file, error)


    ! write primitive variables in parallel
    ! create datatspace for primitive vars
    call H5Screate_simple_f(rank_V, glb_dims_V, dspace_id, error)!, NULL)
    ! create memspace (datachunk (local)) for primitive vars
    call H5Screate_simple_f(rank_V, dims_V, mspace_id, error)!, NULL)

    ! create chunked dataset
    ! init plist for data transfer
    call H5Pcreate_f(H5P_DATASET_CREATE_F, plist_dset, error)
    call H5Pset_chunk_f(plist_dset, rank_V, dims_V, error)
    ! create dataset in collective mode
    call H5Dcreate_f(file_id, dset_prim, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error, plist_dset)


    ! select hyperslab
    xyoffset = (/ bl_i, bl_j /)
    xyoffset = xyoffset-1
    hs_start(1) = 0
    hs_start(2) = xyoffset(1)*gr_nx
    hs_start(3) = xyoffset(2)*gr_ny
    hs_count  = (/ 1, 1, 1 /)
    hs_stride = (/ 1, 1, 1 /)

    call H5Sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hs_start, hs_count, error, hs_stride, dims_V)!, NULL);

    ! collective writing
    call H5Pcreate_f(H5P_DATASET_XFER_F, plist_dxfer, error)
    call H5Pset_dxpl_mpio_f(plist_dxfer, H5FD_MPIO_COLLECTIVE_F, error)

    call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, img_V(1:NSYS_VAR,:,:), &
                    glb_dims_V, error, &
                    file_space_id=dspace_id, mem_space_id=mspace_id, xfer_prp = plist_dxfer)

    call H5Pclose_f(plist_dxfer, error)
    call H5Pclose_f(plist_dset, error)

    call H5Sclose_f(mspace_id, error)
    call H5Sclose_f(dspace_id, error)

    !close dataset
    call H5Dclose_f(dset_id, error)


    ! TODO: below block should be written in serial.
    ! if (bl_ID == 1) then

      allocate(xCoord(gr_glb_nx))
      allocate(yCoord(gr_glb_ny))

      do i = 1, gr_glb_nx
        xCoord(i) = (real(i)-0.5)*gr_dx + gr_xbeg(XDIM)
      end do
      do j = 1, gr_glb_ny
        yCoord(j) = (real(j)-0.5)*gr_dy + gr_xbeg(YDIM)
      end do

      rank_XYZ = 1
      ! dataspace for x-coordinates
      dims_XYZ = gr_glb_nx
      call H5Screate_simple_f(rank_XYZ, dims_XYZ, dspace_id, error)
      call H5Dcreate_f(file_id, dset_x, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
      ! write data
      call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xCoord, dims_XYZ, error)
      ! close dataspace && dataset for x-coordinates
      call H5Dclose_f(dset_id,error)
      call H5Sclose_f(dspace_id,error)

      ! dataspace for y-coordinates
      dims_XYZ = gr_glb_ny
      call H5Screate_simple_f(rank_XYZ, dims_XYZ, dspace_id, error)
      call H5Dcreate_f(file_id, dset_y, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
      ! write data
      call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, yCoord, dims_XYZ, error)
      ! close dataspace && dataset for y-coordinates
      call H5Dclose_f(dset_id,error)
      call H5Sclose_f(dspace_id,error)

      deallocate(xCoord)
      deallocate(yCoord)

      ! write simulation infos
      call io_writeHDF5_simInfo(file_id, init_file_name, t, nstep, eTime, commit_hash)

    ! end if

    ! done!
    ! close file and fortran interface
    call h5fclose_f(file_id, error)
    call h5close_f(error)

    deallocate(img_V)

910 format(i5)

  end subroutine io_writeHDF5_p



  subroutine io_writeHDF5(init_file_name, t, nstep, ioCounter, eTime, commit_hash)

    ! TODO: I should rewrite this using MPI communications

    use HDF5
    use io_hdf5_misc, ONLY: io_writeHDF5_simInfo

    implicit none

    character(len=MAX_STRING_LENGTH), intent(IN) :: init_file_name, commit_hash
    real, intent(IN) :: t, eTime
    integer, intent(IN) :: nstep, ioCounter

    integer :: i, j, dest
    character(len=200) :: ofile
    character(len=5)  :: cCounter

    integer :: ibeg, jbeg, iend, jend
    real, dimension(gr_glb_nx) :: xCoord
    real, dimension(gr_glb_ny) :: yCoord
    real, dimension(NUMB_VAR, gr_glb_nx, gr_glb_ny) :: V

    character(len=50) :: dset_prim, dset_x, dset_y
    integer(HID_T)    :: file_id, dspace_id, dset_id
    integer           :: error, rank_V, rank_XYZ
    integer(HSIZE_T),allocatable, dimension(:) :: dims_V, dims_XYZ, dims_B

    integer:: ierr, tag
    integer, dimension(MPI_STATUS_SIZE) :: stat

    tag = 1
    if (bl_ID == 0) then
      ! I am root; receiving data
      dset_prim = "prim_vars"
      dset_x    = "xCoord"
      dset_y    = "yCoord"
      do i = 1, gr_glb_nx
        xCoord(i) = (real(i)-0.5)*gr_dx + gr_xbeg(XDIM)
      end do
      do j = 1, gr_glb_ny
        yCoord(j) = (real(j)-0.5)*gr_dy + gr_xbeg(YDIM)
      end do

      V(:,1:gr_nx,1:gr_ny) = gr_V(:,gr_ibeg(XDIM):gr_iend(XDIM), gr_ibeg(YDIM):gr_iend(YDIM))

      if (bl_nProcs > 1) then
        do i = 1, bl_iProcs
          do j = 1, bl_jProcs
            dest = bl_grid(i,j)
            if (dest /= 0) then

              ibeg = (i-1)*gr_nx + 1
              jbeg = (j-1)*gr_ny + 1
              iend = ibeg-1 + gr_nx
              jend = jbeg-1 + gr_ny

              ! receive from other processors
              call MPI_Recv( V(:,ibeg:iend, jbeg:jend), NUMB_VAR*gr_nx*gr_ny, MPI_DOUBLE_PRECISION, &
                             dest, tag, MPI_COMM_WORLD, stat, ierr )
            end if
          end do
        end do
      end if

      ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
      !
      ! convert conter number to character
      write(cCounter,910) ioCounter + 10000

      ! file name for ascii output
      !ofile = 'slug_'//trim(sim_name)//'_'//cCounter//'.dat'
      ofile = trim(sim_name)//'_'//cCounter//'.slug'
      allocate(dims_V(3))
      allocate(dims_B(4))
      allocate(dims_XYZ(1))
      rank_XYZ = 1
      dims_V = (/ NSYS_VAR,  gr_glb_nx, gr_glb_ny/)
      rank_V = 3

      !hdf5 fortran interface
      call h5open_f(error)
      !create data file
      call h5fcreate_f(ofile, H5F_ACC_TRUNC_F, file_id, error)

      !datatspace for primitive vars
      call h5screate_simple_f(rank_V, dims_V, dspace_id, error)
      !dataset for primitive vars
      call h5dcreate_f(file_id, dset_prim, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
      !write prim vars
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, V(DENS_VAR:NSYS_VAR,:,:), dims_V, error)

      !close dataset
      call h5dclose_f(dset_id, error)
      !close dataspace
      call h5sclose_f(dspace_id, error)


      !gp betas
      ! rank_V = 4
      ! dims_B = (/ NSYS_VAR, gr_glb_nx, gr_glb_ny, NDIM /)
      ! dset_prim = 'betas'
      ! call h5screate_simple_f(rank_V, dims_B, dspace_id, error)
      ! call h5dcreate_f(file_id, dset_prim, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
      ! call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, gr_betas(:,gr_ibeg(XDIM):gr_iend(XDIM),gr_ibeg(YDIM):gr_iend(YDIM),:), dims_B, error)
      !  !close dataset
      ! call h5dclose_f(dset_id, error)
      ! !close dataspace
      ! call h5sclose_f(dspace_id, error)

      !dataspace for x-coordinates
      dims_XYZ = (/gr_glb_nx/)
      call h5screate_simple_f(rank_XYZ, dims_XYZ, dspace_id, error)
      call h5dcreate_f(file_id, dset_x, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
      !write data
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xCoord, dims_XYZ, error)

      call h5dclose_f(dset_id,error)
      call h5sclose_f(dspace_id,error)

      !dataspace for x-coordinates
      dims_XYZ = gr_glb_ny
      call h5screate_simple_f(rank_XYZ, dims_XYZ, dspace_id, error)
      call h5dcreate_f(file_id, dset_y, h5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
      !write data
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, yCoord, dims_XYZ, error)

      call h5dclose_f(dset_id,error)
      call h5sclose_f(dspace_id,error)

      ! write simulation infos
      call io_writeHDF5_simInfo(file_id, init_file_name, t, nstep, eTime, commit_hash)

      !close hdf5 file and interface
      call h5fclose_f(file_id,error)
      call h5close_f(error)

      deallocate(dims_V)
      deallocate(dims_B)
      deallocate(dims_XYZ)

    else
      ! I am NOT root; sending data
      call MPI_SEND( gr_V(:,gr_ibeg(XDIM):gr_iend(XDIM), gr_ibeg(YDIM):gr_iend(YDIM)), &
                     NUMB_VAR*gr_nx*gr_ny, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, ierr)

    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

910 format(i5)

  end subroutine io_writeHDF5


  subroutine io_writeASCII(ioCounter)
    implicit none

    integer, intent(IN) :: ioCounter

    integer :: i, j, nVar, dest
    character(len=200) :: ofile
    character(len=5)  :: cCounter

    integer :: ibeg, jbeg, iend, jend
    real, dimension(gr_glb_nx) :: xCoord
    real, dimension(gr_glb_ny) :: yCoord
    real, dimension(NUMB_VAR, gr_glb_nx, gr_glb_ny) :: V

    integer:: ierr, tag
    integer, dimension(MPI_STATUS_SIZE) :: stat

    tag = 1
    if (bl_ID == 0) then
      !I am root; receiving data
      do i = 1, gr_glb_nx
        xCoord(i) = (real(i)-0.5)*gr_dx + gr_xbeg(XDIM)
      end do
      do j = 1, gr_glb_ny
        yCoord(j) = (real(j)-0.5)*gr_dy + gr_xbeg(YDIM)
      end do

      V(:,1:gr_nx,1:gr_ny) = gr_V(:,gr_ibeg(XDIM):gr_iend(XDIM), gr_ibeg(YDIM):gr_iend(YDIM))

      if (bl_nProcs > 1) then
        do i = 1, bl_iProcs
          do j = 1, bl_jProcs
            dest = bl_grid(i,j)
            if (dest /= 0) then

              ibeg = (i-1)*gr_nx + 1
              jbeg = (j-1)*gr_ny + 1
              iend = ibeg-1 + gr_nx
              jend = jbeg-1 + gr_ny

              ! receive from other processors
              call MPI_Recv( V(:,ibeg:iend, jbeg:jend), NUMB_VAR*gr_nx*gr_ny, MPI_DOUBLE_PRECISION, &
                             dest, tag, MPI_COMM_WORLD, stat, ierr )
            end if
          end do
        end do
      end if

      ! convert conter number to character
      write(cCounter,910) ioCounter + 10000

      ! file name for ascii output
      !ofile = 'slug_'//trim(sim_name)//'_'//cCounter//'.dat'
      ofile = trim(sim_name)//'_'//cCounter//'.dat'

      open(unit=20,file=ofile,status='unknown')
      do i=1, gr_glb_nx
        do j = 1, gr_glb_ny
          write(20,920)xCoord(i),ycoord(j),(V(nVar,i,j),nVar=1,EINT_VAR)
        end do
        write(20,100)
      end do
      close(20)

    else
      ! I am NOT root; sending data
      call MPI_SEND( gr_V(:,gr_ibeg(XDIM):gr_iend(XDIM), gr_ibeg(YDIM):gr_iend(YDIM)), &
                     NUMB_VAR*gr_nx*gr_ny, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, ierr)

    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

100 format()
910 format(i5)
920 format(1x,f16.8,1x,f16.8,1x,NUMB_VAR f32.16)

  end subroutine io_writeASCII


end module io
