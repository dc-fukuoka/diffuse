! ref: http://ushiro.jp/program/iter/cgsym.htm
! Crank-Nicolson method
! time direction second order accuracy
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
  real(8),allocatable,dimension(:,:,:)::a,anew,r,rnew,p,pnew,x,xnew,ap
  real(8)::alpha,beta,coef1,coef2
  real(8)::r2,rnew2,pap,b2,eps
  real(8)::t=0.0d0
  real(8)::dclock,t0,time
  integer::i,j,k,iter,tstep
#ifdef _CN
  real(8),allocatable,dimension(:,:,:)::ba
  real(8)::coef3,coef4
#endif
  
  read(11,param0)
  read(11,param1)
  read(11,param3)

  allocate(a(0:imax+1,0:jmax+1,0:kmax+1))
  allocate(anew(0:imax+1,0:jmax+1,0:kmax+1))
  allocate(r(0:imax+1,0:jmax+1,0:kmax+1))
  allocate(rnew(0:imax+1,0:jmax+1,0:kmax+1))
  allocate(p(0:imax+1,0:jmax+1,0:kmax+1))
  allocate(pnew(0:imax+1,0:jmax+1,0:kmax+1))
  allocate(x(0:imax+1,0:jmax+1,0:kmax+1))
  allocate(xnew(0:imax+1,0:jmax+1,0:kmax+1))
  allocate(ap(0:imax+1,0:jmax+1,0:kmax+1))
#ifdef _CN
  allocate(ba(0:imax+1,0:jmax+1,0:kmax+1))
#endif
  
  write(6,"(a,4(1pe14.5))") "dt,dx,tol,diff_coef:",dt,dx,tol,diff_coef
  write(6,"(a,4i5)") "imax,jmax,kmax:",imax,jmax,kmax
  write(6,"(a,2i8)") "iter_max,tstep_max:",iter_max,tstep_max
  write(6,"(a,i4)" ) "freq_write:",freq_write
  b2 = 0.0d0

#ifdef _CN
  ! 1/2 factor is needed for Crank-Nicolson method
  coef1 = -0.5d0*dt/dx/dx*diff_coef
  coef2 = 1.0d0+6.0d0/2*dt/dx/dx*diff_coef
  coef3 = -1.0d0*coef1
  coef4 = 1.0d0-6.0d0/2*dt/dx/dx*diff_coef
  write(6,*) "using Crank-Nicolson method"
#else
  coef1 = -1.0d0*dt/dx/dx*diff_coef
  coef2 = 1.0d0+6.0d0*dt/dx/dx*diff_coef
#endif

  open(unit=unit_read,file="data_in",form="unformatted",access="stream")

  ! read initial data
  read(unit_read) a
  close(unit_read)

#ifdef _DEBUG
  k=kmax/2
  do j=1,jmax
     do i=1,imax
        write(599,*) i,j,a(i,j,k)
        if (i.eq.imax+1) write(599,*)
     end do
  end do
#endif

  t0 = dclock()
  open(unit=unit_write,file="data_out",form="unformatted",access="stream")
  do tstep=1,tstep_max
     if (mod(tstep,tstep_max/10).eq.0) write(6,*) "tstep:",tstep
#ifdef _CN
     ! A*anew = B*a
     ! calculate B*a
     do k=1,kmax
        do j=1,jmax
           do i=1,imax
              ! left diagonal precondition
              ! Aanew = Ba
              ! P = diag(Anew)
              ! P^{-1}Aanew = P^{-1}Ba
              ! diag(P^{-1})= coef2
              ! diag(B) = coef4
              ba(i,j,k) = coef3/coef2*(a(i+1,j,  k  )+a(i-1,j,  k  )  &
                                      +a(i,  j+1,k  )+a(i,  j-1,k  )  &
                                      +a(i,  j,  k+1)+a(i,  j,  k-1)) &
                          +coef4/coef2*a(i,  j,  k  )
           end do
        end do
     end do
#endif
     b2 = 0.0d0
     !$omp parallel private(i,j,k)
     !$omp do
     do k=0,kmax+1
        do j=0,jmax+1
           do i=0,imax+1
#ifdef _CN
              x(i,j,k) = ba(i,j,k) ! initial guess
#else
              ! left diagonal precondition
              ! P^{-1}Aanew = P^{-1}a
              ! rhs = P^{-1}a
              x(i,j,k) = a(i,j,k)/coef2 ! initial guess
#endif
           end do
        end do
     end do
     !$omp end do

     ! g = i+(imax+2)*j+(imax+2)*(jmax+2)*k
     !             i+1       i-1       j+1              j-1              k+1                       k-1
     ! coef1*(anew(g+1)+anew(g-1)+anew(g+(imax+2))+anew(g-(imax+2))+anew(g+(imax+2)*(jmax+2))+anew(g-(imax+2)*(jmax+2)))+coef2*anew(g) = a(g)
     !
     ! Crank-Nicolson method:
     ! coef1*(anew(g+1)+anew(g-1)+anew(g+(imax+2))+anew(g-(imax+2))+anew(g+(imax+2)*(jmax+2))+anew(g-(imax+2)*(jmax+2)))+coef2*anew(g) = ba(g)
     ! ba(g) = coef3*(a(g+1)+a(g-1)+a(g+(imax+2))+a(g-(imax+2))+a(g+(imax+2)*(jmax+2))+anew(g-(imax+2)*(jmax+2)))+coef4*a(g)
     !

     !$omp do reduction(+:b2)
     do k=1,kmax
        do j=1,jmax
           do i=1,imax
#ifdef _CN
              ! left diagonal precondition
              r(i,j,k) = ba(i,j,k)-(coef1/coef2*(x(i+1,j,  k  )+x(i-1,j,  k  )  &
                                                +x(i,  j+1,k  )+x(i,  j-1,k  )  &
                                                +x(i,  j,  k+1)+x(i,  j,  k-1)) &
                                                +x(i,  j,  k  ))
              b2       = b2 + ba(i,j,k)*ba(i,j,k)
#else
              ! left diagonal precondition
              ! Aanew = a
              ! P^{-1}Aanew = P^{-1}a
              r(i,j,k) = a(i,j,k)/coef2-(coef1/coef2*(x(i+1,j,  k  )+x(i-1,j,  k  )  &
                                                     +x(i,  j+1,k  )+x(i,  j-1,k  )  &
                                                     +x(i,  j,  k+1)+x(i,  j,  k-1)) &
                                                     +x(i,  j,  k  ))
              b2       = b2 + a(i,j,k)*a(i,j,k)/coef2**2
#endif
              p(i,j,k) = r(i,j,k)
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel
     ! exclude t loop
     do iter=1,iter_max
        rnew2 = 0.0d0
        r2    = 0.0d0
        pap   = 0.0d0
        alpha = 0.0d0
        beta  = 0.0d0
        eps   = 0.0d0
        !     ap    = 0.0d0
        !     xnew  = 0.0d0
        !     rnew  = 0.0d0
        !$omp parallel private(i,j,k)
        !$omp do
        do k=1,kmax
           do j=1,jmax
              do i=1,imax
                 ! left diagonal precondition
                 ap(i,j,k) = coef1/coef2*(p(i+1,j,  k  )+p(i-1,j,  k  )  &
                                         +p(i,  j+1,k  )+p(i,  j-1,k  )  &
                                         +p(i,  j,  k+1)+p(i,  j,  k-1)) &
                                         +p(i,  j,  k  )
              end do
           end do
        end do
        !$omp end do
        !$omp do reduction(+:r2,pap)
        do k=1,kmax
           do j=1,jmax
              do i=1,imax
                 r2  = r2  + r(i,j,k)*r(i,j,k)
                 pap = pap + p(i,j,k)*ap(i,j,k)
              end do
           end do
        end do
        !$omp end do
        !     write(6,*) "iter,r2,pap:",iter,r2,pap
        !$omp single
        alpha = r2/pap
        !$omp end single
        !$omp do reduction(+:rnew2)
        do k=1,kmax
           do j=1,jmax
              do i=1,imax
                 xnew(i,j,k) = x(i,j,k)+alpha*p(i,j,k)
                 rnew(i,j,k) = r(i,j,k)-alpha*ap(i,j,k)
                 rnew2       = rnew2+rnew(i,j,k)*rnew(i,j,k)
              end do
           end do
        end do
        !$omp end do
        !$omp end parallel
        eps = sqrt(rnew2)/sqrt(b2)
        if (eps.le.tol) then
           !        write(6,*) "iter,residual,tol:",iter,eps,tol
           !$omp parallel do
           do k=1,kmax
              do j=1,jmax
                 do i=1,imax
                    anew(i,j,k) = xnew(i,j,k)
                 end do
              end do
           end do
           exit
        else if (iter.eq.iter_max.and.eps.ge.tol) then
           write(6,*) "did not converge. residual,tol:",eps,tol
           stop
        end if
        beta = rnew2/r2
        if (isnan(alpha).or.isnan(beta)) then
           write(6,*) "alpha,beta:",alpha,beta
           stop
        end if

        ! if (mod(iter,iter_max/(imax+2)/(jmax+2)).eq.0) then
        !    write(6,*) "iter,eps:",iter,eps
        ! end if
        !$omp parallel private(i,j,k)
        !$omp do
        do k=1,kmax
           do j=1,jmax
              do i=1,imax
                 pnew(i,j,k) = rnew(i,j,k)+beta*p(i,j,k)
              end do
           end do
        end do
        !$omp end do
        !$omp do
        do k=1,kmax
           do j=1,jmax
              do i=1,imax
                 x(i,j,k) = xnew(i,j,k)
                 r(i,j,k) = rnew(i,j,k)
                 p(i,j,k) = pnew(i,j,k)
              end do
           end do
        end do
        !$omp end do
        !$omp end parallel
     end do ! iter

     !$omp parallel do
     do k=1,kmax
        do j=1,jmax
           do i=1,imax
              a(i,j,k) = anew(i,j,k)
           end do
        end do
     end do

     if (mod(tstep,tstep_max/freq_write).eq.0) then
        write(unit_write) a(1:imax,1:jmax,1:kmax)
     end if
  end do ! tstep

  close(unit_write)
  time = dclock()-t0
  write(6,*) "time[s]:",time

#ifdef _DEBUG
  k=kmax/2
  do j=1,jmax
     do i=1,imax
        write(600,*) i,j,a(i,j,k)
        if (i.eq.imax+1) write(600,*)
     end do
  end do
#endif

  deallocate(a,anew,r,rnew,p,pnew,x,xnew,ap)
#ifdef _CN
  deallocate(ba)
#endif
  
  stop
end program main
