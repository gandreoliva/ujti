!! Use the preprocessor! (compile with ztoa_make.sh)

!! Calculates the time of arrival + gravitational redshift
!! for an arbitrary inclination angle 'incl' w.r.t. the horizontal.
!! The 'x' axis is used as a reference non-inclined position

#define BWFX 1
#define SHFT 2
#define KAFT 3
#define KBLN 4

#define BL_COORD 1
#define SPH_COORD 2


program toa
  use, intrinsic :: iso_fortran_env, dp=>real64
  use, intrinsic :: ieee_arithmetic
  use geq_sphax_base
  use METRIC
  use nsl_solver
  implicit none
  real(dp), parameter :: pi = 3.1415192d0
  character(len=128) :: dataid
  real(dp), dimension(0:3,0:3) :: g
  real(dp), dimension(1:3) :: pos, k
  real(dp), dimension(0:3) :: x,u
  real(dp), dimension(1:3,1:3) :: rotmaty
  real(dp) :: incldeg, ns_rotfreq, ns_re

  integer :: outputfile,metafile,logfile,i,j,imax,imin
  real(dp) :: z, rot, incl
  logical :: on_surf
  real(dp) :: nan
  !! uses nan to fill redshifts where there is no image
  nan = ieee_value(1d0,ieee_quiet_nan)

  !! Binding of functions for the specific metric into the numeric modules
  selected_pot_func => get_potential
  selected_dpot_func => get_deriv_potential
  selected_get_metric => get_metric
  selected_geqrhs => geqrhs

#if USECOORD == BL_COORD
  selected_surf => surface_bl
#elif USECOORD == SPH_COORD
  selected_surf => surface_sph
#endif


  read(*,'(a)') dataid
  read(*,*) incldeg

  !! Logfile and binary output file
  open(newunit=nsl_solver__logfile, file="log.txt", status="replace", action="write")
  open(newunit=outputfile, file="data/zt/"//trim(dataid)//".dat", form="unformatted", access="stream", status="replace")
  open(newunit=metafile, file="data/zt/meta/"//trim(dataid)//".txt", status="replace")


  !! Binary output format:
  !! columns 1,2,3 : 3D Cartesian position of photons as seen from the observer
  !! column 4: "z" (redshift)
  !! columns 5,6,7: 3D Cartesian position of the photon at the surface of the NS
  !! "nan" is used to fill columns 4-7 when the photon does not end up at the
  !! surface of the NS
  !! column 8: final x(0) coordinate


  !! ====================================
  !! Parameters
  !! Precomputed with RNS


#if OBJ == SHFT
  ! ! Config SHFT
   m = 0.185445
   a = 0.0612915
   q = 0.00272009
   asp_ratio = 0.897217
   ns_re = 11.4355
   ns_rotfreq = 716
   rot = 2*pi*ns_rotfreq/3d8*ns_re*1d3

#elif OBJ == BWFX
  ! Config BWFX
    m = 0.283632
    a = 0.0542427
    q = 0.00124465
    asp_ratio = 0.967529
    ns_re = 9.48716
    ns_rotfreq = 622
    rot = 2*pi*ns_rotfreq/3d8*ns_re*1d3

#elif OBJ == KAFT
  ! Config KAFT
  m = 0.164782
  a = 0.0866435
  q = 0.00509572
  asp_ratio = 0.756836
  ns_rotfreq = 1000
  ns_re = 12.6457
  rot = 2*pi*ns_rotfreq/3d8*ns_re*1d3

#elif OBJ == KALN
  ! Config KALN
  m = 0.205571
  a = 0.142883
  q = 0.0092869
  asp_ratio = 0.570313
  ns_rotfreq = 1000
  ns_re = 20.0559
  rot = 2*pi*ns_rotfreq/3d8*ns_re*1d3

#endif



#if TR == 1
  !! Frutos
   q = -q
#elif TR == 2
  !! Mink
   m = 0d0; a = 0d0; q = 0d0; rot = 0d0
#elif TR == 3
  ! GB, HT (no transf)
#endif


  !! Inclination of the observer w.r.t. the horizontal plane
  !! (using the 'x' axis as a reference, this is simply
  !! the angle between the normal to the observer's plane
  !! and the 'x' axis)
  incl = incldeg*pi/180

  !! 8< - - - - - - - - - - - - - - - - -


#if USECOORD == BL_COORD
  !! Binding of "a" for coordinate transformations with Boyer-Lindquist
  coord_bl_a = a
  !! use Boyer-Lindquist coordinates
  call set_coordinates(coord_bl)
  surf_r_eq = 1d0
  surf_bl_a = a
#elif USECOORD == SPH_COORD
  call set_coordinates(coord_sph)
  surf_r_eq = 1d0
#endif

  !! Computational domain
  imax = 640/8
  imin = -640/8


  !! Rotation matrix in 'y' for the initial position and direction
  rotmaty(1,1) = cos(-incl)
  rotmaty(1,3) = sin(-incl)
  rotmaty(2,2) = 1
  rotmaty(3,1) = -sin(-incl)
  rotmaty(3,3) = cos(-incl)

  !! Cartesian 3D initial direction of photons
  k = [-1d0,0d0,0d0]
  k = matmul(rotmaty,k)


  include "print_art.f90"
  print*, "Gravitational redshift calculation (all domain)"
  print*, "data id:", dataid

  do i = imin,imax
    !! Progress output
    write(*,'(a,a,f5.1,a,$)') char(13), "computing...", (imax+i)/(1.0*abs(imin)+abs(imax))*100.0, "%"
    do j = imin,imax
      !! Cartesian 3D position of photons as seen from the observer (start of geodesics)
      pos = [100d0,8*0.003125*real(j,dp),8*0.003125*real(i,dp)]
      pos = matmul(rotmaty,pos)

      x(0) = 0.0
      x(1:3) = pos3d_cart_to_selected(pos)

      u(0) = 1d0
      u(1:3) = vel3d_cart_to_selected(k,x(1:3))
      call normalize_u(x,u)

      ! call grsh_cart_rkf(x,u,1d0,rot,150d0,1d-8, 5d0,1d-10,z,on_surf)
      !            x,u,phot_freq,rot,lmax,dlmin,dlmax,tol,z,on_surf
      call grsh_cart_rkf_tl(x,u,1d0,rot,150d0,1d-8, 5d0,1d-10,z,on_surf,600.)

      !call grsh_cart_eul(x,u,1d0,rot,5d0,1d-3,z,on_surf)
      !grsh_cart_eul(x,u,om,rot,lmax,dl,z,on_surf)

      if (i == imin .and. j == imin) then
        write(*,'(a,$)') '_'
      end if

      !! If the surface is reached at the end of the calculation,
      !! output the corresponding "z", if not, output nan
      if ( on_surf .eqv. .true. ) then

#if USECOORD == BL_COORD
        write(outputfile) pos,z,pos3d_bl_to_cart(x(1:3)),x(0)
#elif USECOORD == SPH_COORD
        write(outputfile) pos,z,pos3d_sph_to_cart(x(1:3)),x(0)
#endif
      else
#if USECOORD == BL_COORD
        write(outputfile) pos,nan,pos3d_bl_to_cart(x(1:3)),nan
#elif USECOORD == SPH_COORD
        write(outputfile) pos,nan,pos3d_sph_to_cart(x(1:3)),nan
#endif
      end if


    end do
  end do
  print*, ""

  print '(a)', "[Report]"
  print*, " dataid = ", dataid
  print*, " m = ", m
  print*, " a = ", a
  print*, " q = ", q
  print*, " asp_ratio", asp_ratio
  print*, " rot = ", rot
  print*, " nx,ny = ", abs(imax - imin) + 1
  print*, "8< - - - - - - - - - - - - - - -"


  write(metafile,*) abs(imax - imin) + 1
  write(metafile,*) m
  write(metafile,*) a
  write(metafile,*) q
  write(metafile,*) asp_ratio
  write(metafile,*) ns_re
  write(metafile,*) ns_rotfreq

end program
