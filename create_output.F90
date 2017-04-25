program main
  implicit none
  integer::imax,jmax,kmax
  integer::iter_max,tstep_max,freq_write
  namelist/param1/imax,jmax,kmax
  namelist/param3/iter_max,tstep_max,freq_write
  integer::i,j,k,l
  real(8),allocatable,dimension(:,:,:)::a
  real(8)::diff
  integer::unit_read=10, unit_write=2222
  integer::len
  real(8),parameter::sigma=8.0d0
  

  read(11,param1)
  read(11,param3)
  open(unit=unit_read,file="data_out",form="unformatted",access="stream")

  allocate(a(1:imax,1:jmax,1:kmax))
  write(6,*) "imax,jmax,kmax:",jmax,jmax,kmax
  write(6,*) "iter_max,tstep_max,freq_write:",iter_max,tstep_max,freq_write
  do l=1,freq_write
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
           if (i.eq.imax.and.j.eq.jmax) then
              write(unit_write,*)
              write(unit_write,*)
           end if
        end do
     end do
  end do

  write(9999) a
  close(unit_read)
  deallocate(a)
  stop
end program main
