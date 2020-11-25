!==================================================================
! Author: Victor Gustavo May Custodio.-----------------------------
! Date : 22/11/2020 -----------------------------------------------
module eigenval
  !  Toolbox for solving the eigenvalue-problems by the Jacobi method.
  !  It also includes the sorting subrouting to enlist all eigenvalues
  !    in ascending order.
  implicit none
  private
  integer, parameter :: dp = kind(0.d0)
  public Jacobi
  public sort
contains
  subroutine Jacobi(H,psi,error,n)
  !===========================================================
  ! ----------------INPUT-------------------------------------
  ! H(n,n) - Discretized Hamiltonian matrix.
  ! n      - mesh (dimension of discretized Hilbert space).
  ! error  - Abs tolerance [sum of (off-diagonal elements)^2].
  ! -----------------OUTPUT----------------------------------
  ! H(i,i)   - Eigenvalues.
  ! psi(i,j) - Eigenvectors by columns!
  !===========================================================
  implicit none
  integer i, j, k, n !counters and dimension
  double precision H(n,n),psi(n,n)
  double precision error, offd, bar
  double precision beta, coeff, c, s, cs, sc  !auxiliary local variables
  ! Initialize as the Identity matrix delta_(ij)
  psi = 0.0_dp
  do i=1,n
    psi(i,i) = 1.0_dp
  end do
  ! find the sum of all off-diagonal elementss suqared (b2)
  offd = 0.0_dp
  do i=1,n
    do j=1,n
      if (i.ne.j) then
        offd = offd + H(i,j)**2
      end if
    end do
  end do

  if (offd .le. error) return

  ! average for off-diagonal elements /2
  bar = 0.5_dp*offd/float(n*n)

  do while (offd.gt.error)
    do i=1,n-1
      do j=i+1,n
        if (H(j,i)**2 .le. bar) then
           cycle  ! skip 'basically zero' entries
        end if
        offd = offd - 2.0_dp*H(j,i)**2
        bar = 0.5_dp*offd/float(n*n)
  ! calculate coefficient c and s for Givens matrix
        beta = (H(j,j)-H(i,i))/(2.0_dp*H(j,i))
        coeff = 0.5_dp*beta/sqrt(1.0_dp+beta**2)
        s = sqrt(max(0.5_dp+coeff,0.0_dp))
        c = sqrt(max(0.5_dp-coeff,0.0_dp))
  ! recalculate rows i and j
        do k=1,n
          cs =  c*H(i,k)+s*H(j,k)
          sc = -s*H(i,k)+c*H(j,k)
          H(i,k) = cs
          H(j,k) = sc
        end do
  ! new matrix H_(k+1) from H_(k), and eigenvectors
        do k=1,n
          cs =  c*H(k,i)+s*H(k,j)
          sc = -s*H(k,i)+c*H(k,j)
          H(k,i) = cs
          H(k,j) = sc
          cs =  c*psi(k,i)+s*psi(k,j)
          sc = -s*psi(k,i)+c*psi(k,j)
          psi(k,i) = cs
          psi(k,j) = sc
        end do
      end do
    end do
  end do
  return
end subroutine

recursive subroutine sort(a, first, last)
  !===========================================================
  ! ----------------INPUT-------------------------------------
  ! a(first:last,2) - list to be sorted./ second index of array
  !                   keeps track of the original indices
  ! first         - first index of array to be sorted
  ! last          - last index of array to be sorted
  ! -----------------OUTPUT----------------------------------
  ! a(first:last) - sorted array
  !===========================================================
  implicit none
  double precision  a(:,:), x, t, s
  integer first, last
  integer i, j
  x = a((first+last)/2,1)
  i = first
  j = last
  do
     do while (a(i,1) < x)
        i=i+1
     end do
     do while (x < a(j,1))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i,1)
     s = a(i,2)
     a(i,1) = a(j,1)
     a(i,2) = a(j,2)
     a(j,1) = t
     a(j,2) = s
     i=i+1
     j=j-1
  end do
  if (first < i-1) call sort(a, first, i-1)
  if (j+1 < last)  call sort(a, j+1, last)
end subroutine
end module
