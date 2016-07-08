module error

contains
 subroutine neo_abort(message)
  
  use memory, only : dealloc_all
  implicit none

  character (len=*),intent(in) :: message

  write(*,*) 'NEOART STOPPED: ', message
  call dealloc_all()
  stop 1

 end subroutine
 
 subroutine neo_warn(message)
  
  implicit none
  character (len=*),intent(in) :: message
  write(*,*) '-------------WARNING!!!---------------'
  write(*,*) message
  write(*,*) '--------------------------------------'


 end subroutine

end module
