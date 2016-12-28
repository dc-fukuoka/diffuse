program main
  implicit none
  integer::imax,jmax,kmax,iter_max,tstep_max
  namelist/param1/imax,jmax,kmax,iter_max,tstep_max
  integer::i,j,k
  real(8),allocatable,dimension(:,:,:)::a,b
  real(8)::diff
  integer::unit=10
  integer::len
  real(8),parameter::sigma=8.0d0
  

  read(11,param1)
  open(unit=unit,file="data_in",form="unformatted",access="stream")

  allocate(a(0:imax+1,0:jmax+1,0:kmax+1),b(0:imax+1,0:jmax+1,0:kmax+1))
  write(6,*) "imax,jmax,kmax:",jmax,jmax,kmax
  do k=0,kmax+1
     do j=0,jmax+1
        do i=0,imax+1
           !  a(i,j,k) = i*10000+j*10+k
           !  a(i,j,k) = sin(dble(i))+cos(dble(j))+sin(dble(k))
           a(i,j,k) = exp(-1.0d0*((i-imax/2)**2+(j-jmax/2)**2+(k-kmax/2)**2)/2/sigma/sigma)/sqrt(2*3.14159d0*sigma) ! exp(-r^2)
           write(unit) a(i,j,k)
        end do
     end do
  end do
#if 0
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
#endif
  k=kmax/2
  do j=1,jmax
     do i=1,imax
        write(100,*) i,j,a(i,j,k)
        if (i.eq.imax+1) write(100,*)
     end do
  end do
  close(unit)
  deallocate(a,b)
  stop
end program main
