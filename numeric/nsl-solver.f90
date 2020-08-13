module nsl_solver
  use, intrinsic :: iso_fortran_env, dp=>real64
  use coord
  implicit none

  ! Warning:
  !   - "vel" refers to the derivative of position w.r.t. the affine parameter,
  !     it is not, in general, a physical velocity!

  integer :: nsl_solver__logfile ! Logfile for errors

  private :: get_metric, geqrhs, surface, pos3d_conv, vel3d_conv
  ! Conection to functions provided by a specific metric
  abstract interface
    function get_metric(pos) result(metric)
      import dp
      real(dp), dimension(0:3,0:3) :: metric
      real(dp), dimension(0:3) :: pos
    end function
    function geqrhs(pos, vel)
      import dp
      real(dp), dimension(0:3) :: geqrhs
      real(dp), dimension(0:3) :: pos
      real(dp), dimension(0:3) :: vel
    end function
    function surface(pos,padding)
      import dp
      real(dp) :: padding
      real(dp), dimension(0:3) :: pos
      logical :: surface
    end function
  end interface

  ! Function selection pointers: the user can later select which metric
  !   to use by pointing these to other module.
  procedure(get_metric), pointer :: selected_get_metric => null()
  procedure(geqrhs), pointer :: selected_geqrhs => null()
  procedure(surface), pointer :: selected_surf => null()

contains

  subroutine grsh_cart_rkf_tl(x,u,om,rot,lmax,dlmin,dlmax,tol,z,on_surf,tl)
    ! GR frequency shift for 1 geodesic cartesian -> specific coords
    !
    ! Arguments
    ! --------
    !   x : initial position (specif. coord) of the photon
    !   u : initial dx/dl (specif. coord) of the photon
    !   om : angular frequency of the photon (can be safely set to 1d0)
    !   lmax : maximum value of the affine parameter
    !   dlmin : minimum step size for the affine parameter
    !   dlmax : maximum step size for the affine parameter
    !   tol : tolerance (RKF method)
    !   rot : angular velocity of the NS. t[g]/R = c*t[SI]/R
    !       => rot[g] = rot[SI]/c = 2*pi*f/c. In units of R = 1, f -> f*R
    !       => (rot[g]*R) = 2*pi*f*R/c
    !   tl : timeout for the calculation
    !
    ! Returns
    ! -------
    !   x : the last position (specif. coord) of the photon
    !   u : the last der. of position (specif. coord)
    !   z : redshift
    !   on_surf : true if the photon reached the surface of the NS

    real(dp), dimension(0:3), intent(inout) :: x,u
    real(dp), intent(in) :: lmax,om,dlmin,dlmax,tol,rot
    real(dp), intent(out) :: z
    logical, intent(out) :: on_surf
    real(dp) :: resc,dl,l,delta,om_i,om_f
    real(dp), dimension(0:3) :: k1,k2,k3,k4,k5,k6,rvec
    real(dp), dimension(0:3) :: ku1,ku2,ku3,ku4,ku5,ku6,ruvec
    logical :: flag
    integer :: iterations
    real(dp), dimension(0:3,0:3) :: g
    real(dp), dimension(0:3) :: vem ! velocity of the emitter
    real :: tl
    real :: starttime, endtime

    ! i: initial, f: final in the integration
    ! RKF coefficients
    real(dp), parameter :: c21=1/4d0, c31=3/32d0, c32=9d0/32, c41=1932/2197d0,&
   & c42=-7200/2197d0, c43=7296/2197d0, c51=439/216d0, c52=-8d0, c53=3680/513d0,&
   & c54=-845/4104d0, c61=-8/27d0, c62=2d0, c63=-3544/2565d0, c64=1859/4104d0,&
   & c65=-11/40d0
    ! -----------------------------------------------------------

    call cpu_time(starttime)

    ! initial values
    l = 0d0
    dl = dlmax
    flag = .true.
    iterations = 0

    ! initial values of position and velocity
    ! Cartesian to spherical coordinates

    g = selected_get_metric(x)

    om_i = (g(0,0)*u(0) + g(0,1)*u(1) + g(0,2)*u(2) + g(0,3)*u(3))*sqrt(-1/g(0,0))
    on_surf = .false.
    !write(outputfile) logical(.false.,1), r(1), r(2), r(3)

    do while(flag .eqv. .true.)
      call cpu_time(endtime)
      if ((endtime-starttime) > tl ) then
        flag = .false.
      end if

        ! k's of RKF method
        k1=dl*selected_geqrhs(x,u)
        k2=dl*selected_geqrhs(x+c21*k1,u+c21*k1)
        k3=dl*selected_geqrhs(x+c31*k1+c32*k2,u+c31*k1+c32*k2)
        k4=dl*selected_geqrhs(x+c41*k1+c42*k2+c43*k3,&
            &u+c41*k1+c42*k2+c43*k3)
        k5=dl*selected_geqrhs(x+c51*k1+c52*k2+c53*k3+c54*k4,&
            &u+c51*k1+c52*k2+c53*k3+c54*k4)
        k6=dl*selected_geqrhs(x+c61*k1+c62*k2+c63*k3+c64*k4+c65*k5,&
            &u+c61*k1+c62*k2+c63*k3+c64*k4+c65*k5)
        ku1=dl*u
        ku2=dl*(u+c21*k1)
        ku3=dl*(u+c31*k1+c32*k2)
        ku4=dl*(u+c41*k1+c42*k2+c43*k3)
        ku5=dl*(u+c51*k1+c52*k2+c53*k3+c54*k4)
        ku6=dl*(u+c61*k1+c62*k2+c63*k3+c64*k4+c65*k5)
        ! used to determine if the integration continues
        rvec = (1/dl)*(1/360d0*k1-128/4275d0*k3-2197/75240d0*k4&
            &+1/50d0*k5+2/55d0*k6)
        ruvec = (1/dl)*(1/360d0*ku1-128/4275d0*ku3-2197/75240d0*ku4&
            &+1/50d0*ku5+2/55d0*ku6)
        resc = sqrt(dot_product(rvec,rvec)+dot_product(ruvec,ruvec))
        if (resc <= tol) then
            ! r <= tol means it's OK to output
            l = l + dl
            u = u + (25/216d0)*k1+(1408/2565d0)*k3+(2197/4104d0)*k4-(1/5d0)*k5
            x = x + (25/216d0)*ku1+(1408/2565d0)*ku3+(2197/4104d0)*ku4-(1/5d0)*ku5

            ! if the surface is reached...
            if (selected_surf(x,0.00d0) .eqv. .true. ) then
              flag = .false.
              on_surf = .true.
              g = selected_get_metric(x)
              vem = [sqrt(-1/(g(0,0)+2*g(0,3)*rot+rot*rot*g(3,3))),0d0,0d0,&
                  & rot*sqrt(-1/(g(0,0)+2*g(0,3)*rot+rot*rot*g(3,3))) ]
              om_f = (g(0,0)*u(0)+g(0,3)*u(3))*vem(0) + (g(3,0)*u(0)+g(3,3)*u(3))*vem(3)
              !---
              ! r_f = [x(1)*sin(x(2))*cos(x(3)),x(1)*sin(x(2))*sin(x(3)),x(1)*cos(x(2))]
              z = om_f/om_i - 1
              !write(outputfile) logical(.true.,1), -1d0, 0d0, 0d0
              return
            end if

            ! output of result
            ! r = pos3d_selected_to_cart(x(1:3))
            ! k = vel3d_selected_to_cart(u(1:3),x(1:3))
          !  write(outputfile) logical(.false.,1), r(1), r(2), r(3)

        end if


        ! Adjustment of the next step
        delta = 0.84d0*(tol/resc)**(1/4d0)
        if (delta <= 0.1) then
            dl = 0.1d0*dl
        else if (delta >= 4) then
            dl = 4*dl
        else
            dl = delta*dl
        end if
        if (dl > dlmax) then
            dl = dlmax
        end if

        if ( l >= lmax ) then
            flag = .false.
        else if (dl<dlmin) then

            flag = .false.
        end if
        iterations = iterations + 1


    end do
    ! write(outputfile) logical(.true.,1), 1d0, 0d0, 0d0
  end subroutine

  ! subroutine grsh_cart_eul(x,u,om,rot,lmax,dl,z,on_surf)
  !
  !   real(dp), dimension(0:3), intent(inout) :: x,u
  !   real(dp), intent(in) :: lmax,dl,om,rot
  !   real(dp), intent(out) :: z
  !   logical, intent(out) :: on_surf
  !   real(dp) :: om_i,om_f,l
  !   integer :: iterations
  !   real(dp), dimension(0:3,0:3) :: g
  !   real(dp), dimension(0:3) :: vem ! velocity of the emitter
  !   ! i: initial, f: final in the integration
  !   ! -----------------------------------------------------------
  !
  !   ! initial values
  !   l = 0.0
  !   iterations = 0
  !
  !   ! initial values of position and velocity
  !   ! Cartesian to spherical coordinates
  !
  !   g = selected_get_metric(x)
  !
  !   om_i = (g(0,0)*u(0) + g(0,1)*u(1) + g(0,2)*u(2) + g(0,3)*u(3))*sqrt(-1/g(0,0))
  !   on_surf = .false.
  !   !write(outputfile) logical(.false.,1), r(1), r(2), r(3)
  !
  !   do while(l < lmax)
  !       ! k's of RKF method
  !           l = l + dl
  !           u = u + selected_geqrhs(x,u)*dl
  !           x = x + u*dl
  !
  !           ! if the surface is reached...
  !           if (selected_surf(x,0.001d0) .eqv. .true. ) then
  !             on_surf = .true.
  !             g = selected_get_metric(x)
  !             vem = [sqrt(-1/(g(0,0)+2*g(0,3)*rot+rot*rot*g(3,3))),0d0,0d0,&
  !                 & rot*sqrt(-1/(g(0,0)+2*g(0,3)*rot+rot*rot*g(3,3))) ]
  !             om_f = (g(0,0)*u(0)+g(0,3)*u(3))*vem(0) + (g(3,0)*u(0)+g(3,3)*u(3))*vem(3)
  !             z = om_f/om_i - 1
  !             return
  !           end if
  !
  !   end do
  ! end subroutine



  

  subroutine normalize_u(x,u)
    ! Normalization of velocity for null geodesics (ds^2 == 0)
    real(dp), dimension(0:3), intent(inout) :: x,u
    real(dp), dimension(0:3,0:3) :: g
    real(dp) :: A, g0iui, gijuiuj
    g = selected_get_metric(x)
    g0iui = g(0,1)*u(1) + g(0,2)*u(2) + g(0,3)*u(3)
    gijuiuj = g(1,1)*u(1)*u(1)+g(1,2)*u(1)*u(2)+g(1,3)*u(1)*u(3)+&
        & g(2,2)*u(2)*u(2)+g(2,1)*u(1)*u(2)+g(2,3)*u(2)*u(3)+&
        & g(3,1)*u(3)*u(1)+g(3,2)*u(3)*u(2)+g(3,3)*u(3)*u(3)
    A = (-g0iui*u(0)+sqrt((g0iui*u(0))**2-4*(gijuiuj)*(g(0,0)*&
        &u(0)*u(0)))) / (2d0*gijuiuj)
    u(1) = A*u(1)
    u(2) = A*u(2)
    u(3) = A*u(3)
  end subroutine normalize_u
end module
