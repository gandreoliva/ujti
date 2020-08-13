module geq_sphax_base
  use, intrinsic :: iso_fortran_env, dp=>real64
  implicit none
  real(dp) :: asp_ratio ! rpolar / requatorial
  real(dp) :: surf_r_eq
  real(dp) :: surf_bl_a
  interface
    function pot_generic(pos) result(pot)
      import dp
      real(dp), dimension(0:3) :: pos
      real(dp), dimension(1:5) :: pot
      ! Index convention:
      !   pot: {1: V, 2: W, 3: X, 4: Y, 5: Z}
    end function
    function dpot_generic(pos) result(dpot)
      import dp
      real(dp), dimension(0:3) :: pos
      real(dp), dimension(1:5,1:2) :: dpot
      ! Index convention:
      !   dpot, first index:
      !     {1: V, 2: W, 3: X, 4: Y, 5: Z}
      !   dpot, second index: (variable w.r.t. which we differentiate)
      !     {1: r, 2: th}
    end function
  end interface
  procedure(pot_generic), pointer :: selected_pot_func => null()
  procedure(dpot_generic), pointer :: selected_dpot_func => null()
  save
contains
  function get_metric(pos) result(metric)
    real(dp), dimension(0:3,0:3) :: metric
    real(dp), dimension(0:3) :: pos
    real(dp), dimension(1:5) :: pot
    pot = selected_pot_func(pos)
    include "sphax-metric.f90"
  end function
  function geqrhs(pos, vel)
    real(dp), dimension(0:3) :: geqrhs
    real(dp), dimension(0:3) :: pos
    real(dp), dimension(0:3) :: vel
    real(dp), dimension(1:5) :: pot
    real(dp), dimension(1:5,1:2) :: dpot
    pot = selected_pot_func(pos)
    dpot = selected_dpot_func(pos)
    include "sphax-geq.f90"
  end function
  function surface_sph(pos,padding)
    real(dp) :: padding
    real(dp), dimension(0:3) :: pos
    logical :: surface_sph
    if (pos(1)/surf_r_eq - sin(pos(2))**2 - asp_ratio*cos(pos(2))**2 < padding) then
      surface_sph = .true.
    else
      surface_sph = .false.
    end if
  end function
  function surface_bl(pos,padding)
    real(dp) :: padding
    real(dp), dimension(0:3) :: pos
    logical :: surface_bl
    real(dp) :: r_sph, th_sph
    associate( r_bl => pos(1), th_bl => pos(2), a=> surf_bl_a )
      r_sph = sqrt( r_bl**2 + a**2*sin(th_bl)**2 )
      th_sph = acos( r_bl/r_sph * cos(th_bl) )
      if (r_sph/surf_r_eq - sin(th_sph)**2 - asp_ratio*cos(th_sph)**2 < padding ) then
        surface_bl = .true.
      else
        surface_bl = .false.
      end if
    end associate
  end function
end module
