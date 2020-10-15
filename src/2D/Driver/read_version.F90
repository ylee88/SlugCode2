subroutine read_version(commit_hash)

#include "definition.h"

  character(len=MAX_STRING_LENGTH), intent(IN OUT) :: commit_hash

  integer :: uid

  uid = 10

  open(unit=uid, file='GIT_VERSION', status='unknown', action='read')
  read(uid, '(A)') commit_hash
  close(uid)

end subroutine read_version
