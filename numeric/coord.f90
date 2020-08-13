module coord
  use, intrinsic :: iso_fortran_env, dp => real64
  implicit none
  real(dp) :: coord_bl_a ! Boyer Lindquist parameter 'a'
  
  ! Coordinate identifiers
  integer, parameter :: coord_sph = 1
  integer, parameter :: coord_bl = 2

  ! For coordinate transformations selection
  abstract interface
    function pos3d_conv(r)
      import dp
      real(dp), dimension(1:3) :: r, pos3d_conv
    end function
    function vel3d_conv(v,r)
      import dp
      real(dp), dimension(1:3) :: v, r, vel3d_conv
    end function
  end interface


  procedure(pos3d_conv), pointer :: pos3d_cart_to_selected => null()
  procedure(pos3d_conv), pointer :: pos3d_selected_to_cart => null()
  procedure(vel3d_conv), pointer :: vel3d_cart_to_selected => null()
  procedure(vel3d_conv), pointer :: vel3d_selected_to_cart => null()

contains

  ! Coordinate selector
  subroutine set_coordinates(coord_syst)
    integer, intent(in) :: coord_syst
    select case(coord_syst)
    case(coord_sph)
      pos3d_cart_to_selected => pos3d_cart_to_sph
      pos3d_selected_to_cart => pos3d_sph_to_cart
      vel3d_cart_to_selected => vel3d_cart_to_sph
      vel3d_selected_to_cart => vel3d_sph_to_cart
    case(coord_bl)
      pos3d_cart_to_selected => pos3d_cart_to_bl
      pos3d_selected_to_cart => pos3d_bl_to_cart
      vel3d_cart_to_selected => vel3d_cart_to_bl
      vel3d_selected_to_cart => vel3d_bl_to_cart
    end select
  end subroutine

  ! Spherical coordinates
  function pos3d_cart_to_sph(r_cart) result(r_sph)
    real(dp), dimension(1:3) :: r_cart, r_sph
    associate(&
      & x=>r_cart(1), y=>r_cart(2), z=>r_cart(3), &
      & r=> r_sph(1), th=> r_sph(2), ph=>r_sph(3) &
      &)
      r = sqrt(x*x + y*y + z*z)
      th = acos(z/r)
      ph = atan2(y,x)
    end associate
  end function

  function pos3d_sph_to_cart(r_sph) result(r_cart)
    real(dp), dimension(1:3) :: r_cart, r_sph
    associate(&
      & x=>r_cart(1), y=>r_cart(2), z=>r_cart(3), &
      & r=> r_sph(1), th=> r_sph(2), ph=>r_sph(3) &
      &)
      x = r*sin(th)*cos(ph)
      y = r*sin(th)*sin(ph)
      z = r*cos(th)
    end associate
  end function

  function vel3d_cart_to_sph(v_cart,r_sph) result(v_sph)
    real(dp), dimension(1:3) :: v_cart, r_sph, v_sph
    associate(&
      & vx=>v_cart(1), vy=>v_cart(2), vz=>v_cart(3), &
      & r=> r_sph(1), th=> r_sph(2), ph=>r_sph(3), &
      & vr=> v_sph(1), vth=>v_sph(2), vph=>v_sph(3) &
      &)
      vr = (vy*sin(ph)+vx*cos(ph))*sin(th)+vz*cos(th)
      vth = -(vz*sin(th)+((-vy*sin(ph))-vx*cos(ph))*cos(th))/r
      vph = -(vx*sin(ph)-vy*cos(ph))/(r*sin(th))
    end associate
  end function

  function vel3d_sph_to_cart(v_sph,r_sph) result(v_cart)
    real(dp), dimension(1:3) :: v_cart, v_sph, r_sph
    associate(&
      & vx=>v_cart(1), vy=>v_cart(2), vz=>v_cart(3), &
      & r=> r_sph(1), th=> r_sph(2), ph=>r_sph(3), &
      & vr=> v_sph(1), vth=>v_sph(2), vph=>v_sph(3) &
      &)
      vx = cos(ph)*r*cos(th)*vth+cos(ph)*sin(th)*vr-sin(ph)*r*sin(th)*vph
      vy = sin(ph)*r*cos(th)*vth+sin(ph)*sin(th)*vr+cos(ph)*r*sin(th)*vph
      vz = cos(th)*vr-r*sin(th)*vth
    end associate
  end function

  ! Boyer-Lindquist coordinates
  function pos3d_cart_to_bl(r_cart) result(r_bl)
    real(dp), dimension(1:3) :: r_cart, r_bl
    real(dp) :: b
    associate(&
      & x=>r_cart(1), y=>r_cart(2), z=>r_cart(3), &
      & r=> r_bl(1), th=> r_bl(2), ph=>r_bl(3), &
      & a=>coord_bl_a &
      &)
      b = x**2 + y**2 + z**2 - a**2
      r = sqrt( (b+sqrt(b**2 + 4*a**2*z**2))/2 )
      th = acos(z/r)
      ph = atan2(y,x)
    end associate
  end function

  function pos3d_bl_to_cart(r_bl) result(r_cart)
    real(dp), dimension(1:3) :: r_bl, r_cart
    associate(&
      & x=>r_cart(1), y=>r_cart(2), z=>r_cart(3), &
      & r=> r_bl(1), th=> r_bl(2), ph=>r_bl(3), &
      & a=>coord_bl_a &
      &)
      x = sqrt(r**2 + a**2)*sin(th)*cos(ph)
      y = sqrt(r**2 + a**2)*sin(th)*sin(ph)
      z = r*cos(th)
    end associate
  end function

  function vel3d_cart_to_bl(v_cart,r_bl) result(v_bl)
    real(dp), dimension(1:3) :: v_cart, r_bl, v_bl
    associate(&
      & vx=>v_cart(1), vy=>v_cart(2), vz=>v_cart(3), &
      & r=> r_bl(1), th=> r_bl(2), ph=>r_bl(3), &
      & vr=> v_bl(1), vth=>v_bl(2), vph=>v_bl(3), &
      & a=>coord_bl_a &
      &)
      vr = (cos(th)*(r**2*vz+a**2*vz)+r*sqrt(r**2+a**2)*sin(th)*(sin(ph)*vy+&
        &cos(ph)*vx))/(a**2*cos(th)**2+r**2)
      vth = (sqrt(r**2+a**2)*cos(th)*(sin(ph)*vy+cos(ph)*vx)-r*sin(th)*vz)/(a&
        &**2*cos(th)**2+r**2)
      vph = -(sin(ph)*vx-cos(ph)*vy)/(sqrt(r**2+a**2)*sin(th))
    end associate
  end function

  function vel3d_bl_to_cart(v_bl,r_bl) result(v_cart)
    real(dp), dimension(1:3) :: v_cart, r_bl, v_bl
    associate(&
      & vx=>v_cart(1), vy=>v_cart(2), vz=>v_cart(3), &
      & r=> r_bl(1), th=> r_bl(2), ph=>r_bl(3), &
      & vr=> v_bl(1), vth=>v_bl(2), vph=>v_bl(3), &
      & a=>coord_bl_a &
      &)
      vx = (cos(ph)*(r**2+a**2)*cos(th)*vth+cos(ph)*r*sin(th)*vr-sin(ph)*(r*&
        &*2+a**2)*sin(th)*vph)/sqrt(r**2+a**2)
      vy = (sin(ph)*(r**2+a**2)*cos(th)*vth+sin(ph)*r*sin(th)*vr+cos(ph)*(r*&
        &*2+a**2)*sin(th)*vph)/sqrt(r**2+a**2)
      vz = cos(th)*vr-r*sin(th)*vth
    end associate
  end function

end module
