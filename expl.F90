! ref: http://ushiro.jp/program/iter/cgsym.htm
module params
  implicit none
  real(8)::dt,dx,tol
  real(8)::diff_coef ! diffusion coefficient
  integer::imax, jmax, kmax ! global
  integer::imax_l, jmax_l, kmax_l
  integer::iter_max ! maximum iteration for CG method
  integer::tstep_max ! maximum # of time steps
  integer::freq_write ! frequency of writing the result
  real(8),parameter::sigma=8.0d0
  integer,parameter::unit_read=12,unit_write=13
end module params

program main
  use params
  implicit none
  namelist /param0/dt,dx,tol,diff_coef
  namelist /param1/imax,jmax,kmax
  namelist /param3/iter_max,tstep_max,freq_write
  real(8),allocatable,dimension(:,:,:)::a,anew
  real(8)::alpha,beta,coef1,coef2
  real(8)::r2,rnew2,pap,b2,eps
  real(8)::t=0.0d0
  real(8)::dclock,t0,time
  real(8)::sevenstencil
  integer::i,j,k,iter,tstep

  read(11,param0)
  read(11,param1)
  read(11,param3)

  allocate(a(0:imax+1,0:jmax+1,0:kmax+1))
  allocate(anew(0:imax+1,0:jmax+1,0:kmax+1))
  
  write(6,"(a,4(1pe14.5))") "dt,dx,tol,diff_coef:",dt,dx,tol,diff_coef
  write(6,"(a,4i5)") "imax,jmax,kmax:",imax,jmax,kmax
  write(6,"(a,2i8)") "iter_max,tstep_max:",iter_max,tstep_max
  write(6,"(a,i4)" ) "freq_write:",freq_write
  write(6,"(a,1pe14.5)") "3*diff_coef*dt/dx/dx:",3*diff_coef*dt/dx/dx
  b2 = 0.0d0
  a  = 0.0d0
  anew = 0.0d0

  coef1 = dt/dx/dx*diff_coef
  coef2  = 1.0d0-6.0d0*dt/dx/dx*diff_coef

  open(unit=unit_read,file="data_in",form="unformatted",access="stream")

  ! read initial data
  read(unit_read) a
  close(unit_read)

  k=kmax/2
  do j=1,jmax
     do i=1,imax
        write(599,*) i,j,a(i,j,k)
        if (i.eq.imax+1) write(599,*)
     end do
  end do

  t0 = dclock()
  open(unit=unit_write,file="data_out",form="unformatted",access="stream")
  do tstep=1,tstep_max
     if (mod(tstep,tstep_max/10).eq.0) write(6,*) "tstep:",tstep
     ! explicit method
     !$omp parallel private(i,j,k)
     !$omp do
     do k=1,kmax
        do j=1,jmax
           do i=1,imax
              anew(i,j,k) = coef1*(a(i+1,j,  k  )+a(i-1,j,  k  )  &
                                  +a(i,  j+1,k  )+a(i,  j-1,k  )  &
                                  +a(i,  j,  k+1)+a(i,  j,  k-1)) &
                            +coef2*a(i,  j,  k  )
           end do
        end do
     end do
     !$omp end do
     !$omp do
     do k=1,kmax
        do j=1,jmax
           do i=1,imax
              a(i,j,k) = anew(i,j,k)
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel
     if (mod(tstep,tstep_max/freq_write).eq.0) then
        write(unit_write) a(1:imax,1:jmax,1:kmax)
     end if
  end do ! tstep

  close(unit_write)
  time = dclock()-t0
  write(6,*) "time[s]:",time

  k=kmax/2
  do j=1,jmax
     do i=1,imax
        write(600,*) i,j,a(i,j,k)
        if (i.eq.imax+1) write(600,*)
     end do
  end do

  deallocate(a,anew)
  
  stop
end program main
