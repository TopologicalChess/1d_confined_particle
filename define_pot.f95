  !==================================================================
  ! Author: Victor Gustavo May Custodio.-----------------------------
  ! Date : 22/11/2020 -----------------------------------------------
module define_pot
  ! Define example potentials or defines new polynomial based potentials
  ! All potentials V(x) in an interval of length L must be normalized
  ! by the factor A=hbar^2/m; i.e.
  !     V(x) => V(x)/A.
  ! In this units the 1-D Time independent Schrödiger equation
  ! looks like -y''(x)+ V(x)y(x) = k^2 y(x)
  ! where k is the wave-number of eigenstates
  !   k_n^2 = 2*m*E_n/ hbar^2
  implicit none
  private
  public set_pot
  integer, parameter :: dp = kind(0.d0)
  contains
  subroutine set_pot(x,v,dx)
    !sets the interval of confinement and the potential.
    implicit none
    integer :: n, type = 0, i, example = 0
    character(len=100) :: potfile
    double precision :: v(:), x(:)
    double precision :: c, xmin, xmax, app, xmid
    double precision, intent(out) :: dx
    open(1, file = 'potential.dat')
    n = size(x)
    print*, 'Enter the interval in which the particle is confined (x_min, x_max)'
    read(*,*) xmin, xmax
    dx = (xmax-xmin)/float(n)
    xmid = (xmax + xmin)*0.5_dp
    do i = 1, n
      x(i) = xmin + float(i)*dx
    end do
    do while (type.gt.2 .or.  type.lt.1)
      print*, 'Press 1 to use an example potential'
      print*, 'Press 2 to input a .dat file'
      read(*,*) type
    end do
    if (type .eq. 1) then
      do while(example .gt. 3 .or. example .lt. 1)
        print*, 'Press 1 to use infinite potential well'
        print*, 'Press 2 to use confined quantum harmonic oscillator'
        print*, 'Press 3 to use approximated Dirac delta potential'
        read(*,*) example
      end do
      if (example .eq. 1) then
        v = 0.0
        do i = 1, n
          write(1,*) x(i), v(i)
        end do
      else if (example .eq. 2) then
        do i = 1, n
          v(i) = 0.25_dp*(x(i)-xmid)**2
          write(1,*) x(i), v(i)
        end do
      else if(example .eq. 3) then
        print*, 'We approximate the Dirac delta function by a sharp normal distribution,'
        print*, '  please enter a very small scale factor for this approximation (1.0E-10<sigma<1.0)'
        read(*,*) app
        do i = 1, n
          v(i) = 10*exp(-((x(i)-xmid)/app)**2)/(app*sqrt(4.0_dp*atan(1.0_dp)))
          write(1,*) x(i), v(i)
        end do
      end if
    else if (type .eq. 2) then
      print*, 'Enter file name (no more than 100 characters!)'
      read(*,*) potfile
      print*, potfile
      open(2, file = potfile)
      do i = 1, n
        read(2,*) x(i), v(i)
        write(1,*) x(i), v(i)
      end do
    end if
    close(1)
  end subroutine
end module
