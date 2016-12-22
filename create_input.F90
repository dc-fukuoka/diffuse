program main
  implicit none
  integer::imax,jmax,kmax,iter_max
  namelist/param1/imax,jmax,kmax,iter_max
  integer::i,j,k
  real(8),allocatable,dimension(:,:,:)::a,b
  real(8)::diff
  integer::unit=10
  integer::len


  read(11,param1)
  open(unit=unit,file="data_in",form="unformatted",access="stream")

  allocate(a(0:imax+1,0:jmax+1,0:kmax+1),b(0:imax+1,0:jmax+1,0:kmax+1))
  do k=0,kmax+1
     do j=0,jmax+1
        do i=0,imax+1
           a(i,j,k) = k*100+j*10+i
           write(unit) a(i,j,k)
        end do
     end do
  end do
  rewind(unit)
  do k=0,kmax+1
     do j=0,jmax+1
        do i=0,imax+1
           read(unit) b(i,j,k)
           diff = b(i,j,k)-a(i,j,k)
           if (abs(diff).ge.1.0d-10) then
              write(6,*) "diff:", diff
              stop
           end if
!           write(6,*) "b(",i,j,k,"): ",b(i,j,k)
        end do
     end do
  end do
  close(unit)
  deallocate(a,b)
  stop
end program main
