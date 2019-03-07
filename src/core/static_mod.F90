module static_mod

  use params_mod
  use mesh_mod

  implicit none

  type static_on_cell_type
    real(real_kind), allocatable :: ghs(:)
  end type static_on_cell_type

  type static_type
    type(static_on_cell_type) cell
  end type static_type

  type(static_type) static

contains

  subroutine static_init()

    if (.not. allocated(static%cell%ghs)) allocate(static%cell%ghs(nCells))

  end subroutine static_init

  subroutine static_final()

    if (allocated(static%cell%ghs)) deallocate(static%cell%ghs)

  end subroutine static_final

end module static_mod