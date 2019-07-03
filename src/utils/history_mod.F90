module history_mod

  use params_mod
  use mesh_mod
  use time_mod
  use io_mod
  use log_mod
  use string_mod
  use static_mod
  use state_mod
  use operators_mod

  implicit none

  private

  public history_init
  public history_final
  public history_write

  interface history_write
    module procedure history_write_state
  end interface history_write

contains

  subroutine history_init()

    character(10) time_value, time_units
    real(real_kind) seconds

    time_value = string_split(history_interval(1), 1)
    time_units = string_split(history_interval(1), 2)
    read(time_value, *) seconds
    select case (time_units)
    case ('days')
      seconds = seconds * 86400
    case ('hours')
      seconds = seconds * 3600
    case ('minutes')
      seconds = seconds * 60
    case ('seconds')
      seconds = seconds
    case default
      call log_error('Invalid history interval ' // trim(history_interval(1)) // '!')
    end select

    if (output_file_prefix /= 'N/A') then
      call io_create_dataset('h0', desc=case_name, file_prefix=output_file_prefix)
    else
      call io_create_dataset('h0', desc=case_name, file_prefix=case_name)
    end if

    call io_add_att('h0', 'source',         'TMCORE')
    call io_add_att('h0', 'dt',             dt)
    call io_add_att('h0', 'time_scheme',    time_scheme)
    call io_add_att('h0', 'author',         'N/A')
    call io_add_att('h0', 'on_a_sphere',    'YES')
    call io_add_att('h0', 'sphere_radius',  radius)
    ! Dimensions
    call io_add_dim('h0', 'Time',           add_var=.true.)
    call io_add_dim('h0', 'nCells',         size=nCells)
    call io_add_dim('h0', 'nEdges',         size=nEdges)
    call io_add_dim('h0', 'nVertices',      size=nVertices)
    call io_add_dim('h0', 'TWO',            size=2)
    call io_add_dim('h0', 'vertexDegree',   size=vertexDegree)
    call io_add_dim('h0', 'maxEdges',       size=maxEdges)
    call io_add_dim('h0', 'maxEdges2',      size=maxEdges2)
    ! Mesh parameters
    call io_add_var('h0', 'lonCell',        long_name='Longitude on the cell',                       units='radian', dim_names=['nCells      '],                 data_type='real(8)')
    call io_add_var('h0', 'latCell',        long_name='Latitude on the cell',                        units='radian', dim_names=['nCells      '],                 data_type='real(8)')
    call io_add_var('h0', 'xCell',          long_name='Cartesian X on the cell',                     units='m',      dim_names=['nCells      '],                 data_type='real(8)')
    call io_add_var('h0', 'yCell',          long_name='Cartesian Y on the cell',                     units='m',      dim_names=['nCells      '],                 data_type='real(8)')
    call io_add_var('h0', 'zCell',          long_name='Cartesian Z on the cell',                     units='m',      dim_names=['nCells      '],                 data_type='real(8)')
    call io_add_var('h0', 'indexToCellID',  long_name='Global cell ID',                              units='1',      dim_names=['nCells      '],                 data_type='integer')
    call io_add_var('h0', 'lonEdge',        long_name='Longitude on the edge',                       units='radian', dim_names=['nEdges      '],                 data_type='real(8)')
    call io_add_var('h0', 'latEdge',        long_name='Latitude on the edge',                        units='radian', dim_names=['nEdges      '],                 data_type='real(8)')
    call io_add_var('h0', 'xEdge',          long_name='Cartesian X on the edge',                     units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
    call io_add_var('h0', 'yEdge',          long_name='Cartesian Y on the edge',                     units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
    call io_add_var('h0', 'zEdge',          long_name='Cartesian Z on the edge',                     units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
    call io_add_var('h0', 'indexToEdgeID',  long_name='Global edge ID',                              units='1',      dim_names=['nEdges      '],                 data_type='integer')
    call io_add_var('h0', 'lonVertex',      long_name='Longitude on the vertex',                     units='radian', dim_names=['nVertices   '],                 data_type='real(8)')
    call io_add_var('h0', 'latVertex',      long_name='Latitude on the vertex',                      units='radian', dim_names=['nVertices   '],                 data_type='real(8)')
    call io_add_var('h0', 'xVertex',        long_name='Cartesian X on the vertex',                   units='m',      dim_names=['nVertices   '],                 data_type='real(8)')
    call io_add_var('h0', 'yVertex',        long_name='Cartesian Y on the vertex',                   units='m',      dim_names=['nVertices   '],                 data_type='real(8)')
    call io_add_var('h0', 'zVertex',        long_name='Cartesian Z on the vertex',                   units='m',      dim_names=['nVertices   '],                 data_type='real(8)')
    call io_add_var('h0', 'indexToVertexID',long_name='Global vertex ID',                            units='1',      dim_names=['nVertices   '],                 data_type='integer')
    call io_add_var('h0', 'areaCell',       long_name='Primary cell area',                           units='m2',     dim_names=['nCells      '],                 data_type='real(8)')
    call io_add_var('h0', 'areaTriangle',   long_name='Dual cell area',                              units='m2',     dim_names=['nVertices   '],                 data_type='real(8)')
    call io_add_var('h0', 'areaEdge',       long_name='Defined edge area',                           units='m2',     dim_names=['nEdges      '],                 data_type='real(8)')
    call io_add_var('h0', 'nEdgesOnCell',   long_name='Edge number on the cell',                     units='1',      dim_names=['nCells      '],                 data_type='integer')
    call io_add_var('h0', 'nEdgesOnEdge',   long_name='Edge number to reconstruct tangent velocity', units='1',      dim_names=['nEdges      '],                 data_type='integer')
    call io_add_var('h0', 'cellsOnCell',    long_name='Cell indices that surround cell',             units='1',      dim_names=['maxEdges    ', 'nCells      '], data_type='integer')
    call io_add_var('h0', 'cellsOnEdge',    long_name='Cell indices that saddle cell',               units='1',      dim_names=['TWO         ', 'nEdges      '], data_type='integer')
    call io_add_var('h0', 'cellsOnVertex',  long_name='Cell indices that surround vertex',           units='1',      dim_names=['vertexDegree', 'nVertices   '], data_type='integer')
    call io_add_var('h0', 'edgesOnCell',    long_name='Edge indices on the cell',                    units='1',      dim_names=['maxEdges    ', 'nCells      '], data_type='integer')
    call io_add_var('h0', 'edgesOnEdge',    long_name='Edge indices to reconstruct tangent velocity',units='1',      dim_names=['maxEdges2   ', 'nEdges      '], data_type='integer')
    call io_add_var('h0', 'edgesOnVertex',  long_name='Edge indices on the vertex',                  units='1',      dim_names=['vertexDegree', 'nVertices   '], data_type='integer')
    call io_add_var('h0', 'verticesOnCell', long_name='Vertex indices on the cell',                  units='1',      dim_names=['maxEdges    ', 'nCells      '], data_type='integer')
    call io_add_var('h0', 'verticesOnEdge', long_name='Vertex indices on the edge',                  units='1',      dim_names=['TWO         ', 'nEdges      '], data_type='integer')
    ! Dynamical variables
    call io_add_var('h0', 'u',              long_name='Normal wind on the edge',                     units='m s-1',  dim_names=['nEdges   ', 'Time     '],       data_type='real(8)')
    call io_add_var('h0', 'zonalWind',      long_name='Zonal Wind',                                  units='m s-1',  dim_names=['nCells   ', 'Time     '],       data_type='real(8)')
    call io_add_var('h0', 'merdionalWind',  long_name='Merdional Wind',                              units='m s-1',  dim_names=['nCells   ', 'Time     '],       data_type='real(8)')
    call io_add_var('h0', 'h',              long_name='Geopotential height on the cell',             units='m',      dim_names=['nCells   ', 'Time     '],       data_type='real(8)')
    call io_add_var('h0', 'pv',             long_name='Potential vorticity on the cell',             units='',       dim_names=['nVertices', 'Time     '],       data_type='real(8)')
    call io_add_var('h0', 'div',            long_name='Divergence',                                  units='s-1',    dim_names=['nCells   ', 'Time     '],       data_type='real(8)')
    call io_add_var('h0', 'tm',             long_name='total mass',                                  units='m2 s-2', dim_names=['Time     '],                    data_type='real(8)')
    call io_add_var('h0', 'te',             long_name='total energy',                                units='m4 s-4', dim_names=['Time     '],                    data_type='real(8)')

    call time_add_alert('history_write', seconds=seconds)

    call log_notice('History module is initialized.')

  end subroutine history_init

  subroutine history_final()

    call log_notice('History module is finalized.')

  end subroutine history_final

  subroutine history_write_state(state, static)

    type(state_type),  intent(inout) :: state
    type(static_type), intent(in   ) :: static

    call div_operator(state%edge%u, state%cell%div)

    call io_start_output('h0', time_elapsed_seconds(), new_file=.false.)
    call io_output('h0', 'lonCell',         lonCell)
    call io_output('h0', 'latCell',         latCell)
    call io_output('h0', 'xCell',           xCell)
    call io_output('h0', 'yCell',           yCell)
    call io_output('h0', 'zCell',           zCell)
    call io_output('h0', 'indexToCellID',   indexToCellID)
    call io_output('h0', 'lonEdge',         lonEdge)
    call io_output('h0', 'latEdge',         latEdge)
    call io_output('h0', 'xEdge',           xEdge)
    call io_output('h0', 'yEdge',           yEdge)
    call io_output('h0', 'zEdge',           zEdge)
    call io_output('h0', 'indexToEdgeID',   indexToEdgeID)
    call io_output('h0', 'lonVertex',       lonVertex)
    call io_output('h0', 'latVertex',       latVertex)
    call io_output('h0', 'xVertex',         xVertex)
    call io_output('h0', 'yVertex',         yVertex)
    call io_output('h0', 'zVertex',         zVertex)
    call io_output('h0', 'indexToVertexID', indexToVertexID)
    call io_output('h0', 'areaCell',        areaCell)
    call io_output('h0', 'areaTriangle',    areaTriangle)
    call io_output('h0', 'areaEdge',        areaEdge)
    call io_output('h0', 'nEdgesOnCell',    nEdgesOnCell)
    call io_output('h0', 'nEdgesOnEdge',    nEdgesOnEdge)
    call io_output('h0', 'cellsOnCell',     cellsOnCell)
    call io_output('h0', 'cellsOnEdge',     cellsOnEdge)
    call io_output('h0', 'cellsOnVertex',   cellsOnVertex)
    call io_output('h0', 'edgesOnCell',     edgesOnCell)
    call io_output('h0', 'edgesOnEdge',     edgesOnEdge)
    call io_output('h0', 'edgesOnVertex',   edgesOnVertex)
    call io_output('h0', 'verticesOnCell',  verticesOnCell)
    call io_output('h0', 'verticesOnEdge',  verticesOnEdge)
    call io_output('h0', 'u',               state%edge%u)
    call io_output('h0', 'zonalWind',       state%cell%u)
    call io_output('h0', 'merdionalWind',   state%cell%v)
    call io_output('h0', 'h',               (state%cell%gd + static%cell%ghs) / g)
    call io_output('h0', 'pv',              state%vertex%pv)
    call io_output('h0', 'div',             state%cell%div)
    call io_output('h0', 'tm',              state%total_mass)
    call io_output('h0', 'te',              state%total_energy)
    call io_end_output('h0')

  end subroutine history_write_state

end module history_mod
