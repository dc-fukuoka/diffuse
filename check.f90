program main
  implicit none
  integer::imax,jmax,kmax

  namelist/param1/imax,jmax,kmax
  integer::i,j,k
  real(8),allocatable,dimension(:,:,:)::a,b
  real(8)::diff
  integer::len
  real(8),parameter::sigma=8.0d0
  real(8),parameter::tol=1.0d-12

  read(11,param1)
  allocate(a(imax,jmax,kmax),b(imax,jmax,kmax))
  open(45, file="data_out.orig", form="unformatted", access="stream")
  open(55, file="data_out",      form="unformatted", access="stream")

  read(45) a
  read(55) b

  write(6,*) "imax,jmax,kmax:",jmax,jmax,kmax
  write(6,*) "tol:",tol
  !$omp parallel do private(i,j,k,diff)
  do k=1,kmax
     do j=1,jmax
        do i=1,imax
           diff = abs((b(i,j,k)-a(i,j,k))/a(i,j,k))
           if (diff.ge.tol) then
              write(6,*) "i,j,k:",i,j,k
              write(6,*) "diff:", diff
              stop
           end if
        end do
     end do
  end do
  write(6,*) "check passed."
  close(45)
  close(55)
  deallocate(a,b)
  stop
end program main
