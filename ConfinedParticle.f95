!==================================================================
! Author: Victor Gustavo May Custodio.-----------------------------
! Date : 22/11/2020 -----------------------------------------------
program ConfinedProgram
! This program read a .dat file for a well-scaled potential and computes
!  the eigenstates for a discretized
  !===========================================================
  ! ----------------INPUT-----(by the user)------------------
  ! n            - mesh (dimension of discretized Hilbert space).
  ! [xmin,xmax]  - interval of confinement (1D-box).
  ! v(n)         - potential.
  ! -----------------OUTPUT----------------------------------
  ! 4 .dat files:
  !     (1.) potential.dat   - stores x,v(x).
  !     (2.) wavenumbers.dat - stores wavenumbers k = 2mE/hbar^2.
  !     (3.) eigenstates.dat - stores the states psi_n for each k_n.
  !     (4.) plottingdata.dat - stores the potential and the first 4
  !                             states with their energies.
  ! 1 .plt file (already runned within the program) to plot the first
  !    four states and their energies.
  !   (*N-O-T-E*): This program assumes gnuplot 5.2 for plotting.
  !    -comment the last statement in this program to avoid errors if
  !     you do not have gnuplot
  !===========================================================
use define_pot
use eigenval
implicit none
integer, parameter :: dp = kind(0.d0)
double precision, allocatable :: x(:), v(:), H(:,:), k(:,:), Psi(:,:)
double precision :: dx, tol, L
integer :: n, i, j, m
integer, allocatable :: order(:)
tol = 1.0d-14
print*, 'Enter dimension of the discretized Hilbert Space'
read(*,*) n
allocate(v(n),x(n),H(n,n),k(n,2),Psi(n,n), order(n)) !initiates all arrays
call set_pot(x,v,dx)
! Set up the discretized Hamiltonian H(:,:)
do i = 1, n
  do j = 1, n
    if (i .eq. j) then
      H(i,j) = 2.0_dp/(dx*dx) + v(i)
    else if (i.eq.j-1 .or. i.eq.j+1) then
      H(i,j) = -1.0_dp/(dx**2)
    end if
  end do
end do
!Diagonalize the Hamiltonian and get the eigenvectors Psi
call Jacobi(H, Psi, tol, n)
!Get the eigenvectors from the diagonal
do i = 1, n
  k(i,1) = H(i,i)
  k(i,2) = float(i)
end do
!Order the eigenvalues in ascendent order
call sort(k,1,n)
!Write k
open(unit = 10, file = "wavenumbers.dat")
do i = 1, n
  write(10,*) k(i,1)
end do
close(10)
!Write Psi, but first, reorder them accordingly with the
!  new order of k
!  States are stored by row in the .dat file!! (not by columns)
open(unit = 20, file = "eigenstates.dat",  action="write", status="replace")
do i = 1, n
  do j = 1, n
    if(int(k(i,2)).eq.j) then
      write(20,*)  Psi(:,j)
      order(i) = j
    end if
  end do
end do
close(20)
!Write the plotting file
open(30, file='plottingdata.dat')
do i = 1,  n
write(30,*) x(i), v(i), k(1,1), Psi(i,order(1))  , k(2,1), Psi(i,order(2)), k(3,1), Psi(i,order(3)), k(4,1), Psi(i,order(4))
end do
close(30)
!Plots the first four eigenstates
 call system('gnuplot -p plots.plt')
end program
