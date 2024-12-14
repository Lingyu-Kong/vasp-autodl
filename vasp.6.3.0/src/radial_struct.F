#include "symbol.inc"
MODULE radial_struct_def
    USE prec
    ! structure which is used for the logarithmic grid
    ! the grid points are given by R(i) = RSTART * exp [H (i-1)]
    TYPE rgrid
       REAL(q)  :: RSTART                   ! starting point
       REAL(q)  :: REND                     ! endpoint
       REAL(q)  :: RMAX                     ! radius of augmentation sphere
       REAL(q)  :: D                        ! R(N+1)/R(N) = exp(H)
       REAL(q)  :: H                        !
       REAL(q),POINTER :: R(:)  => NULL()   ! radial grid (r-grid)
       REAL(q),POINTER :: SI(:) => NULL()   ! integration prefactors on r-grid
       INTEGER  :: NMAX                     ! number of grid points
    END TYPE rgrid

    ! This parameter determines at which magnetization the aspherical contributions
    ! to the one center magnetization are truncated in the non collinear case
    !   Without any truncation the aspherical terms for non magnetic atoms
    ! tend to yield spurious but meaningless contributions to the potential
    ! so that convergence to the groundstate can not be achieved
    ! for details see the routines RAD_MAG_DIRECTION and RAD_MAG_DENSITY
    REAL(q), PARAMETER :: MAGMIN=1E-2

    ! for non collinear calculations, setting
    ! the parameter USE_AVERAGE_MAGNETISATION  means that the aspherical 
    ! contributions to the one center magnetisation are projected onto the
    ! average magnetization direction in the PAW sphere instead of the
    ! local moment of the spherical magnetization density at
    ! each grid-point
    ! USE_AVERAGE_MAGNETISATION improves the numerical stability significantly
    ! and must be set
    LOGICAL :: USE_AVERAGE_MAGNETISATION=.TRUE.
END MODULE radial_struct_def
