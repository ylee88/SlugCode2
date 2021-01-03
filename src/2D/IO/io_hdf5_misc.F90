module io_hdf5_misc

#include "definition.h"

  implicit none

contains

  subroutine io_writeHDF5_simInfo(file_id, init_file_name, t, nstep, eTime, commit_hash)

    use HDF5

    implicit none

    integer(HID_T), intent(IN) :: file_id
    real,    intent(IN) :: t, eTime
    integer, intent(IN) :: nstep
    character(len=MAX_STRING_LENGTH), intent(IN) :: init_file_name, commit_hash

    integer(HID_T) :: dspace_id, dset_id
    integer(HSIZE_T), dimension(1) :: dimT
    integer(HSIZE_T), dimension(1) :: dims_txt

    integer(HID_T) :: strtype
    integer(SIZE_T) :: strlen

    character(len=MAX_STRING_LENGTH), dimension(:), allocatable :: init_params
    integer :: file_length
    integer :: rankT, error

    !dataspace for the single element
    rankT = 1
    dimT  = 1
    call H5Screate_simple_f(rankT, dimT, dspace_id, error)

    ! elapsed (CPU) time
    call H5Dcreate_f(file_id, 'eTime', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, eTime, dimT, error)
    call H5Dclose_f(dset_id,error)
    ! simulation time
    call H5Dcreate_f(file_id, 'time', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, t, dimT, error)
    call H5Dclose_f(dset_id,error)
    ! current step number
    call H5Dcreate_f(file_id, 'nStep', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER, nStep, dimT, error)
    call H5Dclose_f(dset_id,error)

    ! close dataspace
    call H5Sclose_f(dspace_id,error)


    ! dataspace for init parameters
    ! allocate contents of init_file_name to init_params
    call text_to_array(init_file_name, file_length, init_params)
    ! extend string datatype
    strlen = MAX_STRING_LENGTH
    call H5Tcopy_f(H5T_FORTRAN_S1, strtype, error)
    call H5Tset_size_f(strtype, strlen, error)
    ! datatspace
    dims_txt = (/ file_length /)
    call H5Screate_simple_f(1, dims_txt, dspace_id, error)
    call H5Dcreate_f(file_id, 'init_params', strtype, dspace_id, dset_id, error)
    ! write data
    call H5Dwrite_f(dset_id, strtype, init_params, dims_txt, error)

    ! call H5Tclose_f(strtype,error)   !! strtype will be used for writing version, save
    call H5Dclose_f(dset_id,error)
    call H5Sclose_f(dspace_id,error)
    deallocate(init_params)


    !dataspace for git version
    dims_txt = (/ 1 /)
    call h5Screate_simple_f(1, dims_txt, dspace_id, error)
    call h5Dcreate_f(file_id, 'version', strtype, dspace_id, dset_id, error)
    !write data
    call h5Dwrite_f(dset_id, strtype, commit_hash, dims_txt, error)

    call H5Tclose_f(strtype,error)
    call H5Dclose_f(dset_id,error)
    call H5Sclose_f(dspace_id,error)

  end subroutine io_writeHDF5_simInfo


  subroutine text_to_array(ofile, line_num, text_array)

    ! open text file, store the contents into string array

    implicit none

    character(len=MAX_STRING_LENGTH), intent(IN) :: ofile
    character(len=MAX_STRING_LENGTH), allocatable, dimension(:), intent(INOUT) :: text_array
    integer, intent(OUT) :: line_num

    character(len=MAX_STRING_LENGTH) :: ctmp
    integer :: i, uid, ierr

    uid = 50

    open(unit=uid, file=ofile, status='unknown')

    ! get line number
    ierr = 0
    line_num = 0
    do while(ierr == 0)
      line_num = line_num + 1
      read(uid, '(A)', iostat=ierr) ctmp
    end do
    line_num = line_num - 1

    allocate(text_array(line_num))

    rewind(uid)
    do i = 1, line_num
      read(uid, '(A)') text_array(i)
    end do

    close(uid)

  end subroutine text_to_array


end module io_hdf5_misc
