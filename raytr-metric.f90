program raytr
  use, intrinsic :: iso_fortran_env, dp=>real64
  use geq_sphax_base
  include "gensrc/raytr-metricname.f90"
  use indiv_ngeod_solver
  implicit none
  real(dp), parameter :: pi = 3.1415192d0
  real(dp), dimension(0:3,0:3) :: g
  real(dp), dimension(1:3) :: pos, k, k0
  integer :: outputfile,logfile,metafile
  real(dp) :: dlmin, tol, sepgeod
  character(len=128) :: dataid, metricname
  integer :: view,i,imax,imin,ifactor

  selected_pot_func => get_potential
  selected_dpot_func => get_deriv_potential
  selected_get_metric => get_metric
  selected_geqrhs => geqrhs
  selected_surf => surface_sph

  indiv_ngeod_solver__logfile = logfile

  read(*,'(a)') dataid
  read(*,'(a)') metricname

  read(*,*) m
  read(*,*) a
  read(*,*) q
  read(*,*) asp_ratio
  read(*,*) view
  read(*,*) tol
  read(*,*) dlmin


  include "gensrc/raytr-paramtransf.f90"


  open(newunit=logfile, file="log.txt", status="replace", action="write")
  open(newunit=outputfile, file="data/raytr/"//trim(dataid)//"/"//trim(metricname)//".dat",&
    & form="unformatted", access="stream", status="replace")
  open(newunit=metafile, file="data/raytr/"//trim(dataid)//"/meta/"//trim(metricname)//".txt",&
    & status="replace", action="write")


  if (asp_ratio < 0 .and. asp_ratio >= -1) then
    !! Calculate aspect ratio, assuming solid body =====
    !!   >>> J2 = (Iz - Ix)/M*R^2 = Q/(M_g*R^2);  R=Req
    !!  If solid ellipsoid, Iz = 2/5*M*R^2, Ix = 1/5*M*(R^2 + Rp^2)
    !!      => (M_g/R) - (M_g/R)*Rp^2/R^2 = (Q/R^3)*5/2
    !!  /!\ The abs is added assuming that's always oblate!
    asp_ratio = sqrt(1-5*abs(q)/(2*m))
  end if

  surf_r_eq = 1d0
  coord_bl_a = a
  call set_coordinates(coord_bl)

  ! print*, get_potential([1d0,10d0,pi/2,0d0])

  include "print_art.f90"

  select case(view)
  case(1) !! xy plane
  !! Horizontal

    ifactor = 4
    imax = 640/ifactor
    imin = -640/ifactor

    do i = imin,imax
      write(*,'(a,a,f5.1,a,$)') char(13), "computing...", (imax+i)/(1.0*abs(imin)+abs(imax))*100.0, "%"

      sepgeod = 0.003125*ifactor
      pos = [15d0,sepgeod*real(i,dp),0d0]
      k = [-1d0,0d0,0d0]
      k0 = k


      !call ngeod_cart_rkf(pos,k,1d0,40d0,1d-5,2d0,1d-7,outputfile)
      ! call ngeod_cart_rkf(pos,k,1d0,50d0,1d-9,2d0,1d-11,outputfile)
      ! call ngeod_cart_rkf(pos,k,1d0,50d0,1d-8,5d0,1d-10,outputfile)
      call ngeod_cart_rkf(pos,k,1d0,50d0,dlmin,5d0,tol,outputfile)
      !write(logfile,*) i,0.25*i,k
      ! !  ngeod_cart_rkf(r  , k ,om,lmax,dlmin,dlmax,tol,outputfile)
      ! print*, acos(dot_product(k0,k)/sqrt(dot_product(k,k)*dot_product(k0,k0)))*180/pi
    end do

  case(2)
    !! Vertical

    ifactor = 64
    imax = 640/ifactor
    imin = -640/ifactor

    do i = imin,imax
      write(*,'(a,a,f5.1,a,$)') char(13), "computing...", (imax+i)/(1.0*abs(imin)+abs(imax))*100.0, "%"

      sepgeod = 0.003125*ifactor
      pos = [6d0,0.02d0,sepgeod*real(i,dp)]
      k = [-1d0,0d0,0d0]
      k0 = k
      ! call ngeod_cart_eul(pos,k,1d0,50d0,0.1d0,outputfile)


      ! call ngeod_cart_rkf(pos,k,1d0,40d0,1d-5,2d0,1d-7,outputfile)
      ! call ngeod_cart_rkf(pos,k,1d0,50d0,1d-8,2d0,1d-12,outputfile)

      call ngeod_cart_rkf(pos,k,1d0,20d0,dlmin,2d0,tol,outputfile)
      ! !  ngeod_cart_rkf(r  , k ,om,lmax,dlmin,dlmax,tol,outputfile)
      ! print*, acos(dot_product(k0,k)/sqrt(dot_product(k,k)*dot_product(k0,k0)))*180/pi
    end do
  end select

print*, ""
! print '(a)', "[Report]"
write(metafile,*) " dataid = ", dataid
write(metafile,*) " metricname = ", metricname
write(metafile,*) " m = ", m
write(metafile,*) " a = ", a
write(metafile,*) " q = ", q
write(metafile,*) " asp_ratio =", asp_ratio
write(metafile,*) " view = ", view
write(metafile,*) " dlmin = ", dlmin
write(metafile,*) " tol = ", tol
write(metafile,*) " imax = ", imax
write(metafile,*) " sepgeod = ", sepgeod
write(metafile,*) "8< - - - - - - - - - - - - - - -"


end program
