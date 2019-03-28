program gen_mesh

  use const_mod
  use log_mod
  use io_mod
  use string_mod
  use array_mod
  use linked_list_mod
  use delaunay_voronoi_mod

  implicit none

  integer, parameter :: nx = 40962
  real(8), allocatable :: x(:,:)

  allocate(x(nx,3))
  
  call io_init()
  call io_create_dataset(name='grids', file_path='./point.' // trim(to_string(nx)) // '.nc', mode='input')
  call io_start_input('grids')
  call io_input('vtx_p', x, dataset_name='grids')
  call io_end_input('grids')

  call delaunay_voronoi_init(nx, x=x(:,1), y=x(:,2), z=x(:,3))
  call delaunay_triangulation(output_delaunay)

  deallocate(x)

contains

  subroutine output_delaunay(DVT_array, DT_list, tag)

    type(array_type), intent(in) :: DVT_array
    type(linked_list_type), intent(in) :: DT_list
    character(*), intent(in), optional :: tag

    integer i, j
    real(8), allocatable :: lon(:)
    real(8), allocatable :: lat(:)
    integer, allocatable :: idx(:,:)
    class(*), pointer :: DVT, DT
    type(linked_list_iterator_type) iterator

    allocate(lon(DVT_array%size))
    allocate(lat(DVT_array%size))
    allocate(idx(3,DT_list%size))

    do i = 1, DVT_array%size
      select type (DVT => DVT_array%value_at(i))
      type is (delaunay_vertex_type)
        lon(i) = DVT%lon * deg
        lat(i) = DVT%lat * deg
      end select
    end do

    i = 1
    iterator = linked_list_iterator(DT_list)
    do while (.not. iterator%ended())
      select type (DT => iterator%value)
      type is (delaunay_triangle_type)
        do j = 1, 3
          select type (DVT => DT%DVT%value_at(j))
          type is (delaunay_vertex_type)
            idx(j,i) = DVT%id
          end select
        end do
      end select
      i = i + 1
      call iterator%next()
    end do

    call io_create_dataset(name='delaunay', file_path='mesh.' // trim(to_string(nx)) // '.nc', mode='output')
    call io_add_dim('num_DVT', size=DVT_array%size, dataset_name='delaunay')
    call io_add_dim('num_DT',  size=DT_list%size,   dataset_name='delaunay')
    call io_add_dim('THREE',   size=3,              dataset_name='delaunay')
    call io_add_var('lon_DVT',    long_name='longitude of DVT',  units='degree_east',  dim_names=['num_DVT'],            data_type='real(8)', dataset_name='delaunay')
    call io_add_var('lat_DVT',    long_name='latitude of DVT',   units='degree_north', dim_names=['num_DVT'],            data_type='real(8)', dataset_name='delaunay')
    call io_add_var('DT_DVT_idx', long_name='DVT indices of DT', units='1',            dim_names=['THREE  ', 'num_DT '], data_type='integer', dataset_name='delaunay')
    call io_start_output(dataset_name='delaunay', tag=tag)
    call io_output('lon_DVT',    lon, dataset_name='delaunay')
    call io_output('lat_DVT',    lat, dataset_name='delaunay')
    call io_output('DT_DVT_idx', idx, dataset_name='delaunay')
    call io_end_output('delaunay')

  end subroutine output_delaunay

end program gen_mesh
