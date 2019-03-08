module history_mod

  use params_mod
  use mesh_mod
  use io_mod
  use log_mod
  use string_mod
  use static_mod
  use state_mod

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

    call io_init()

    if (output_file_prefix /= 'N/A') then
      call io_create_dataset(desc='TMCORE dataset', file_prefix=output_file_prefix, frames_per_file=frames_per_file)
    else
      call io_create_dataset(desc='TMCORE dataset', file_prefix=case_name, frames_per_file=frames_per_file)
    end if

    call io_add_meta('dt', dt)
    call io_add_meta('time_scheme', time_scheme)
    call io_add_meta('author', 'N/A')
    ! Dimensions
    call io_add_dim('time',      add_var=.true.)
    call io_add_dim('nCells',    size=nCells)
    call io_add_dim('nEdges',    size=nEdges)
    call io_add_dim('nVertices', size=nVertices)
    call io_add_dim('two',       size=2)
    call io_add_dim('three',     size=3)
    call io_add_dim('maxEdges',  size=maxEdges)
    ! Mesh parameters
    call io_add_var('lonCell',        long_name='Longitude on the cell',      units='degrees_east',  dim_names=['nCells   '])
    call io_add_var('latCell',        long_name='Latitude on the cell',       units='degrees_north', dim_names=['nCells   '])
    call io_add_var('lonEdge',        long_name='Longitude on the edge',      units='degrees_east',  dim_names=['nEdges   '])
    call io_add_var('latEdge',        long_name='Latitude on the edge',       units='degrees_north', dim_names=['nEdges   '])
    call io_add_var('lonVertex',      long_name='Longitude on the vertex',    units='degrees_north', dim_names=['nVertices'])
    call io_add_var('latVertex',      long_name='Latitude on the vertex',     units='degrees_north', dim_names=['nVertices'])
    call io_add_var('areaCell',       long_name='Primary cell area',          units='m2',            dim_names=['nCells   '])
    call io_add_var('areaTriangle',   long_name='Dual cell area',             units='m2',            dim_names=['nVertices'])
    call io_add_var('areaEdge',       long_name='Defined edge area',          units='m2',            dim_names=['nEdges   '])
    call io_add_var('nEdgesOnCell',   long_name='Edge number on the cell',    units='1',             dim_names=['nCells   '])
    call io_add_var('edgesOnCell',    long_name='Edge indices on the cell',   units='1',             dim_names=['maxEdges ', 'nCells   '])
    call io_add_var('verticesOnEdge', long_name='Vertex indices on the edge', units='1',             dim_names=['two      ', 'nEdges   '])
    call io_add_var('edgesOnVertex',  long_name='Edge indices on the vertex', units='1',             dim_names=['three    ', 'nVertices'])
    ! Dynamical variables
    call io_add_var('u',  long_name='Normal wind on the edge',         units='m s-1',  dim_names=['nEdges   ', 'time     '])
    call io_add_var('h',  long_name='Geopotential height on the cell', units='m',      dim_names=['nCells   ', 'time     '])
    call io_add_var('pv', long_name='Potential vorticity on the cell', units='',       dim_names=['nVertices', 'time     '])
    call io_add_var('tm', long_name='total mass',                      units='m2 s-2', dim_names=['time     '])
    call io_add_var('te', long_name='total energy',                    units='m4 s-4', dim_names=['time     '])

    call log_notice('History module is initialized.')

  end subroutine history_init

  subroutine history_final()

    call log_notice('History module is finalized.')

  end subroutine history_final

  subroutine history_write_state(state, static)

    type(state_type), intent(in) :: state
    type(static_type), intent(in) :: static

    call io_start_output()
    call io_output('lonCell',        lonCell * deg)
    call io_output('latCell',        latCell * deg)
    call io_output('lonEdge',        lonEdge * deg)
    call io_output('latEdge',        latEdge * deg)
    call io_output('lonVertex',      lonVertex * deg)
    call io_output('latVertex',      latVertex * deg)
    call io_output('areaCell',       areaCell)
    call io_output('areaTriangle',   areaTriangle)
    call io_output('areaEdge',       areaEdge)
    call io_output('nEdgesOnCell',   nEdgesOnCell)
    call io_output('edgesOnCell',    edgesOnCell)
    call io_output('verticesOnEdge', verticesOnEdge)
    call io_output('edgesOnVertex',  edgesOnVertex)
    call io_output('u',              state%edge%u)
    call io_output('h',              (state%cell%gd + static%cell%ghs) / g)
    call io_output('pv',             state%vertex%pv)
    call io_output('tm',             state%total_mass)
    call io_output('te',             state%total_energy)
    call io_end_output()

  end subroutine history_write_state

end module history_mod
