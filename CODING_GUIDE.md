# Directory structure

The directory structure should conform to normal standards and be easy to be familiar with.

```
- src:   source codes
  - core
  - utils
  - test_cases
    - sw_model
- lib:   external libraries
- build: out source build directory
- run:   test run directory
```

# Fortran style:

Indentation width should be two spaces and no TAB!

`enddo`, `endif` and similar types should not be used, use `end do` and `end if` instead.

Underscore naming style is preferred, but for compatibility with MPAS, we can tolerate with `edgesOnCell` and so on.
