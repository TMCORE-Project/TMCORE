program gen_mesh

  use log_mod
  use io_mod
  use string_mod
  use delaunay_voronoi_mod

  implicit none

  integer, parameter :: nx = 40962
  real(8), allocatable :: x(:,:)

  allocate(x(nx,3))
  
  call io_init()
  call io_create_dataset(name='grids', file_path='./point.'//trim(to_string(nx))//'.nc', mode='input')
  call io_start_input('grids')
  call io_input('vtx_p', x, dataset_name='grids')
  call io_end_input('grids')

  call delaunay_voronoi_init(nx, x=x(:,1), y=x(:,2), z=x(:,3))
  call delaunay_triangulation()

  deallocate(x)

end program gen_mesh
