program main
  implicit none
  integer::imax,jmax,kmax,iter_max
  namelist/param1/imax,jmax,kmax,iter_max
  integer::i,j,k
  real(8),allocatable,dimension(:,:,:)::a
  real(8)::diff
  integer::unit_read=10, unit_write=2222
  integer::len
  real(8),parameter::sigma=8.0d0
  

  read(11,param1)
  open(unit=unit_read,file="data_out",form="unformatted",access="stream")

  allocate(a(1:imax,1:jmax,1:kmax))
  write(6,*) "imax,jmax,kmax:",jmax,jmax,kmax
  do k=1,kmax
     do j=1,jmax
        do i=1,imax
           read(unit_read) a(i,j,k)
        end do
     end do
  end do
  
  k = kmax/2
  do j=1,jmax
     do i=1,imax
        write(unit_write,*) i,j,a(i,j,k)
        if (i.eq.imax) write(unit_write,*)
     end do
  end do

  write(9999) a
  close(unit_read)
  deallocate(a)
  stop
end program main
