!
!  NURBS library for implementation of Isogeometric Analysis
!
  subroutine find_span(p,n,m,uvec,u,mid)
!
!  Determine the knot span idex
!  Algorithm A2.1 in The NURBS book.
!
!  Input: 
!         p - degree of the curve, 
!         n - number of control points,
!         m - the highest index of knot vector,
!         u - control parameter, 
!         uvec - a knot vector,
!  Output: The knot span index while (u < uvec(mid) || u >= uvec(mid+1));- that is 'i' if u <= u < u
!                                                                                          i        i+1
!
  implicit none 
  integer, intent(in) :: p
  integer, intent(in) :: n
  integer, intent(in) :: m
  double precision, dimension(0:m), intent(in) :: Uvec
  double precision, intent(in) :: u
  integer, intent(out) :: mid
! 
! Locals
!
  integer :: low, high

  if ( u >= Uvec(n+1)) then
    mid = n
    return
  endif
  if ( u<= Uvec(p)) then
    mid = p
    return
  endif

  low = p
  high = n+1
  mid = (low + high)/2

    do while (u < Uvec(mid) .or. u >= Uvec(mid+1))
      if ( u < Uvec(mid) ) then
        high = mid
      else
        low = mid
      endif
      mid = (low + high)/2
    enddo 

  return
  end

  subroutine basis_funs(p,n,m,uvec,u,nbasis)
!
!  Computes the nonvanishing basis functions.
!  Based on Eq. 2.5 and Cox-deBoor-Mansfield reccurence algorithm.
!  Algorithm A2.2 of The NURBS book, pp70. 
!
!  Input: span - specifies the span at which basis function to compute, 
!         u - the parametric value, 
!         p - degree of the curve, 
!         U - the knot vector
!  Output: N - the non-zero basis functions in the array N(0),...N(p) of length p+1.
!          This array has to be allocated in the calling function.
!
!  * Dopunjeno je tako da se unutar funkcije izracunava span
!

  implicit none
  integer, intent(in) :: p
  integer, intent(in) :: n
  integer, intent(in) :: m
  double precision, dimension(0:m), intent(in) :: uvec
  double precision, intent(in) :: u
  double precision, dimension(0:p), intent(out) :: nbasis
! 
! Locals
!
  integer :: i,j,r,span
  double precision :: temp, saved
  double precision, dimension(0:2*(p+1)-1), target :: left
  double precision, dimension(:), pointer :: right
  right => left(p+1:2*(p+1)-1)

  call find_span(p,n,m,uvec,u,span)

  i = span

  nbasis(0) = 1.0d0

  do j=1,p
    left(j) = u-Uvec(i+1-j)
    right(j) = Uvec(i+j)-u
    saved = 0.0d0
    do r=0,j-1 
      temp = nbasis(r)/(right(r+1)+left(j-r))
      nbasis(r) = saved+right(r+1) * temp
      saved = left(j-r) * temp
    end do
    nbasis(j) = saved
  end do

  return
  end

  subroutine basis_funs_and_derivatives(p,n,m,uvec,u,ders,d)
!
!  Compute the basis functions and their derivatives at 'u' of the NURBS curve.
!  Based on Eq. 2.10 and Algorithm A2.3 in The NURBS book, pp 72, 
!
!  The result is stored in the ders matrix, where ders is of 
!  size (d+1,p+1) and the derivative 
!  N'_i(u) = ders(1,i=span-p+j) where j = 0...p+1.
!
!
!  Input:
!  d - the degree of the derivation 
!  u - the parametric value
!  p - degree of the curve
!  span - the span for the basis functions
!  Uvec - the knot vector on which the Basis functions must be computed
! 
!  Output:
!  ders  A (d+1,p+1) matrix containing the basis functions and derivatives up to order d.
!

  implicit none
  integer, intent(in) :: p
  integer, intent(in) :: n
  integer, intent(in) :: m
  double precision, dimension(0:m), intent(in) :: uvec
  double precision, intent(in) :: u
  double precision, dimension(0:d,0:p), intent(out) :: ders
  integer, intent(in) :: d
! 
! Locals
!
  double precision, dimension(0:p) :: left
  double precision, dimension(0:p) :: right

  double precision, dimension(0:p,0:p) :: ndu
  double precision, dimension(0:1,0:p) :: a

  double precision :: saved,temp,dd
  integer :: j,r,span,s1,s2,rk,pk,j1,j2,k

  call find_span(p,n,m,uvec,u,span)

  ndu(0,0) = 1.0;
  do j=1,p
    left(j) = u-Uvec(span+1-j)
    right(j) = Uvec(span+j)-u
    saved = 0.0
    
    do r=0,j-1 
      ! Lower triangle
      ndu(j,r) = right(r+1)+left(j-r)
      temp = ndu(r,j-1)/ndu(j,r)
      ! Upper triangle
      ndu(r,j) = saved+right(r+1) * temp
      saved = left(j-r) * temp
    end do

    ndu(j,j) = saved
  end do

  ders(0,:) = ndu(:,p)  ! Load the basis functions

  ! Compute the derivatives Eq. 2.10 in The NURBS book

  do r=0,p
    s1 = 0 
    s2 = 1 ! alternate rows in array a
    a(0,0) = 1.0
    ! Compute the kth derivative
    do k=1,d
      dd = 0.0d0
      rk = r-k
      pk = p-k

      if(r>=k) then
        a(s2,0) = a(s1,0)/ndu(pk+1,rk)
        dd = a(s2,0)*ndu(rk,pk)
      endif

      if(rk>=-1) then
        j1 = 1
      else 
        j1 = -rk
      endif

      if(r-1 <= pk) then
        j2 = k-1 
      else 
        j2 = p-r
      endif

      do j=j1,j2
        a(s2,j) = (a(s1,j)-a(s1,j-1))/ndu(pk+1,rk+j)
        dd= dd + a(s2,j)*ndu(rk+j,pk)
      end do
      
      if(r<=pk) then
        a(s2,k) = -a(s1,k-1)/ndu(pk+1,r)
        dd = dd + a(s2,k)*ndu(r,pk)
      endif
      ders(k,r) = dd
      j = s1 
      s1 = s2 
      s2 = j ! Switch rows
    end do
  end do

  ! Multiply through by the correct factors
  r = p
  do k=1,d
    do j=0,p
      ders(k,j) = ders(k,j) * r
    end do
    r = r * p-k
  end do

! Strange scaling I have to do - check what bug:
  ders(2,:) = ders(2,:)*p/(p+1.0d0) 

  return
  end

  subroutine basis_funs_cpt(p,n,m,uvec,ncpt,ucpt,nb_cpt)
!
!  Computes the nonvanishing basis functions.
!  Based on Eq. 2.5 and Cox-deBoor-Mansfield reccurence algorithm.
!  Algorithm A2.2 of The NURBS book, pp70. 
!
!  Input: span - specifies the span at which basis function to compute, 
!         u - the parametric value, 
!         p - degree of the curve, 
!         U - the knot vector
!  Output: N - the non-zero basis functions in the array N(0),...N(p) of length p+1.
!          This array has to be allocated in the calling function.
!
!  * Dopunjeno je tako da se unutar funkcije izracunava span
!

  implicit none
  integer, intent(in) :: p
  integer, intent(in) :: n
  integer, intent(in) :: m
  double precision, dimension(0:m), intent(in) :: uvec
  integer, intent(in) :: ncpt
  double precision, dimension(ncpt), intent(in) :: ucpt
  double precision, dimension(ncpt,0:n), intent(out) :: nb_cpt
! 
! Locals
!
  integer :: i,j,r,span,ic
  double precision :: temp, saved
  double precision, dimension(0:p) :: nbasis
  double precision, dimension(0:2*(p+1)-1), target :: left
  double precision, dimension(:), pointer :: right
  right => left(p+1:2*(p+1)-1)

  nb_cpt = 0.0d0

  do ic=1,ncpt
  call find_span(p,n,m,uvec,ucpt(ic),span)

  i = span

  nbasis(0) = 1.0d0

  do j=1,p
    left(j) = ucpt(ic)-Uvec(i+1-j)
    right(j) = Uvec(i+j)-ucpt(ic)
    saved = 0.0d0
    do r=0,j-1 
      temp = nbasis(r)/(right(r+1)+left(j-r))
      nbasis(r) = saved+right(r+1) * temp
      saved = left(j-r) * temp
    end do
    nbasis(j) = saved
  end do

  nb_cpt(ic,i-p:i) = nbasis(:)
  enddo

  return
  end

  subroutine basis_funs_and_derivatives_cpt(p,n,m,uvec,ncpt,ucpt,ders_cpt,d)
!
!  Compute the basis functions and their derivatives at 'u' of the NURBS curve.
!  Based on Eq. 2.10 and Algorithm A2.3 in The NURBS book, pp 72, 
!
!  The result is stored in the ders matrix, where ders is of 
!  size (d+1,p+1) and the derivative 
!  N'_i(u) = ders(1,i=span-p+j) where j = 0...p+1.
!
!
!  Input:
!  d - the degree of the derivation 
!  u - the parametric value
!  p - degree of the curve
!  span - the span for the basis functions
!  Uvec - the knot vector on which the Basis functions must be computed
! 
!  Output:
!  ders  A (d+1,p+1) matrix containing the basis functions and derivatives up to order d.
!

  implicit none
  integer, intent(in) :: p
  integer, intent(in) :: n
  integer, intent(in) :: m
  double precision, dimension(0:m), intent(in) :: uvec
  integer, intent(in) :: ncpt
  double precision, dimension(ncpt), intent(in) :: ucpt
  double precision, dimension(ncpt,0:d,0:n), intent(out) :: ders_cpt
  integer, intent(in) :: d
! 
! Locals
!
  double precision, dimension(0:d,0:p) :: ders

  double precision, dimension(0:p) :: left
  double precision, dimension(0:p) :: right

  double precision, dimension(0:p,0:p) :: ndu
  double precision, dimension(0:1,0:p) :: a

  double precision :: saved,temp,dd
  integer :: j,r,span,s1,s2,rk,pk,j1,j2,k,ic

  ders_cpt = 0.0d0

  do ic=1,ncpt
  call find_span(p,n,m,uvec,ucpt(ic),span)

  ndu(0,0) = 1.0;
  do j=1,p
    left(j) = ucpt(ic)-Uvec(span+1-j)
    right(j) = Uvec(span+j)-ucpt(ic)
    saved = 0.0
    
    do r=0,j-1 
      ! Lower triangle
      ndu(j,r) = right(r+1)+left(j-r)
      temp = ndu(r,j-1)/ndu(j,r)
      ! Upper triangle
      ndu(r,j) = saved+right(r+1) * temp
      saved = left(j-r) * temp
    end do

    ndu(j,j) = saved
  end do

  ders(0,:) = ndu(:,p)  ! Load the basis functions

  ! Compute the derivatives Eq. 2.10 in The NURBS book

  do r=0,p
    s1 = 0 
    s2 = 1 ! alternate rows in array a
    a(0,0) = 1.0
    ! Compute the kth derivative
    do k=1,d
      dd = 0.0d0
      rk = r-k
      pk = p-k

      if(r>=k) then
        a(s2,0) = a(s1,0)/ndu(pk+1,rk)
        dd = a(s2,0)*ndu(rk,pk)
      endif

      if(rk>=-1) then
        j1 = 1
      else 
        j1 = -rk
      endif

      if(r-1 <= pk) then
        j2 = k-1 
      else 
        j2 = p-r
      endif

      do j=j1,j2
        a(s2,j) = (a(s1,j)-a(s1,j-1))/ndu(pk+1,rk+j)
        dd= dd + a(s2,j)*ndu(rk+j,pk)
      end do
      
      if(r<=pk) then
        a(s2,k) = -a(s1,k-1)/ndu(pk+1,r)
        dd = dd + a(s2,k)*ndu(r,pk)
      endif
      ders(k,r) = dd
      j = s1 
      s1 = s2 
      s2 = j ! Switch rows
    end do
  end do

  ! Multiply through by the correct factors
  r = p
  do k=1,d
    do j=0,p
      ders(k,j) = ders(k,j) * r
    end do
    r = r * p-k
  end do

! Strange scaling I have to do - check what bug:
  ders(2,:) = ders(2,:)*p/(p+1.0d0)

! Store everything in an array
  ders_cpt(ic,:,span-p:span) = ders(:,:)
  end do

  return
  end
