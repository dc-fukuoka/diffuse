program main
  implicit none
  integer::imax,jmax,kmax
  
  namelist/param1/imax,jmax,kmax
  integer::i,j,k
  real(8),allocatable,dimension(:,:,:)::a,b
  real(8)::diff
  integer::unit=10
  integer::len
  real(8),parameter::sigma=8.0d0
  real(8),parameter::pi=4*atan(1.0d0)

  read(11,param1)
  open(unit=unit,file="data_in",form="unformatted",access="stream")

  allocate(a(0:imax+1,0:jmax+1,0:kmax+1),b(0:imax+1,0:jmax+1,0:kmax+1))
  write(6,*) "imax,jmax,kmax:",jmax,jmax,kmax
!$omp parallel do private(i,j,k)
  do k=0,kmax+1
     do j=0,jmax+1
        do i=0,imax+1
           !  a(i,j,k) = i*10000+j*10+k
           !  a(i,j,k) = sin(dble(i))+cos(dble(j))+sin(dble(k))
           a(i,j,k) = exp(-1.0d0*((i-imax/2)**2+(j-jmax/2)**2+(k-kmax/2)**2)/2/sigma/sigma)/sqrt(2*pi)/sigma ! exp(-r^2)
        end do
     end do
  end do
  write(unit) a

#ifdef _DEBUG
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
  k=kmax/2
  do j=1,jmax
     do i=1,imax
        write(100,*) i,j,a(i,j,k)
        if (i.eq.imax+1) write(100,*)
     end do
  end do
#endif

  close(unit)
  deallocate(a,b)
  stop
end program main
