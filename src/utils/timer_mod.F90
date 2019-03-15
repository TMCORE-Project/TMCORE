module timer_mod

  use hash_table_mod
  use log_mod

  implicit none

  private

  public timer_init
  public timer_start
  public timer_end
  public timer_duration

  type timer_type
    real(8) start_time
    real(8) end_time
  end type timer_type

  type(hash_table_type) timers

contains

  subroutine timer_init()

    timers = hash_table()

  end subroutine timer_init

  subroutine timer_start(name)

    character(*), intent(in) :: name

    type(timer_type) timer

    call cpu_time(timer%start_time)

    print *, timer%start_time
    call timers%insert(name, timer)

  end subroutine timer_start

  subroutine timer_end(name)

    character(*), intent(in) :: name

    type(timer_type), pointer :: timer

    timer => get_timer(name)
    call cpu_time(timer%end_time)

  end subroutine timer_end

  real(8) function timer_duration(name) result(res)

    character(*), intent(in) :: name

    type(timer_type) timer

    timer = get_timer(name)
    res = timer%end_time - timer%start_time

  end function timer_duration

  function get_timer(name) result(res)

    character(*), intent(in) :: name
    type(timer_type), pointer :: res

    if (timers%hashed(name)) then
      select type (value => timers%value(name))
      type is (timer_type)
        res => value
        return
      class default
        call log_error('Unknown timer ' // trim(name) // '!')
      end select
    end if

  end function get_timer

end module timer_mod