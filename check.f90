program main
  implicit none
  integer::imax,jmax,kmax
  integer::iter_max,tstep_max,freq_write

  namelist/param1/imax,jmax,kmax
  namelist/param3/iter_max,tstep_max,freq_write
  integer::i,j,k,l
  real(8),allocatable,dimension(:,:,:,:)::a,b
  real(8)::diff
  integer::len
  real(8),parameter::sigma=8.0d0
  real(8),parameter::tol=1.0d-12

  read(11,param1)
  read(11,param3)
  allocate(a(imax,jmax,kmax,freq_write),b(imax,jmax,kmax,freq_write))
  open(45, file="data_out.orig", form="unformatted", access="stream")
  open(55, file="data_out",      form="unformatted", access="stream")

  read(45) a
  read(55) b

  write(6,*) "imax,jmax,kmax:",jmax,jmax,kmax
  write(6,*) "freq_write:",freq_write
  write(6,*) "tol:",tol
  !$omp parallel do private(i,j,k,l,diff)
  do l=1,freq_write
     do k=1,kmax
        do j=1,jmax
           do i=1,imax
              diff = abs((b(i,j,k,l)-a(i,j,k,l))/a(i,j,k,l))
              if (diff.ge.tol) then
                 write(6,*) "i,j,k,l:",i,j,k,l
                 write(6,*) "a:",a(i,j,k,l)
                 write(6,*) "b:",b(i,j,k,l)
                 write(6,*) "diff:", diff
                 stop
              end if
           end do
        end do
     end do
  end do
  write(6,*) "check passed."
  close(45)
  close(55)
  deallocate(a,b)
  stop
end program main
