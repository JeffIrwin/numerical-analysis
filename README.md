
![](https://github.com/JeffIrwin/numerical-analysis/workflows/CI/badge.svg)

# numerical-analysis

Exercises from the textbook "Introduction to numerical analysis" by Stoer and
Bulirsch.  The textbook can be purchased here:
https://link.springer.com/book/10.1007/978-0-387-21738-3

<!--
Note to self:

Also see local dir "/mnt/c/Users/jirwi/cpp/math523/" with my original code from
when I took Math 523 at Penn State in 2014
-->

## Dependencies
- `fpm`: https://fpm.fortran-lang.org/
- `gfortran` or another Fortran compiler
- Check the [Dockerfile](Dockerfile) for an up to date list of dependencies

## Run the test application

Run the test suite using `fpm`, the [Fortran package
manager](https://fpm.fortran-lang.org/):
```
fpm test
```

## Build the library

```
fpm build
```

