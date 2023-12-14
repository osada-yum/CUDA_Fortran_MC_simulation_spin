module GPU_simulation
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, GPU_simulation!"
  end subroutine say_hello
end module GPU_simulation
