subroutine read_version(commit_hash)

#include "definition.h"

  character(len=MAX_STRING_LENGTH), intent(IN OUT) :: commit_hash

  integer :: uid, pos
  character(len=MAX_STRING_LENGTH) :: exec_cmd, path

  uid = 10

  ! get exec path
  call get_command_argument(0, exec_cmd)
  pos = scan(exec_cmd, '/', .true.)
  path = exec_cmd(1:pos)

  open(unit=uid, file=TRIM(ADJUSTL(path))//'GIT_VERSION', status='unknown', action='read')
  read(uid, '(A)') commit_hash
  close(uid)

end subroutine read_version
