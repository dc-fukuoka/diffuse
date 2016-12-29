! ref: http://ushiro.jp/program/iter/cgsym.htm
module params
  implicit none
  integer,parameter::imax=64,jmax=64,kmax=64
  integer,parameter::max_iter=(imax+2)*(jmax+2)*(kmax+2),max_tstep=100
  !  integer,parameter::max_iter=1024*1014*128
  real(8),parameter::dt=1.0d-6, dx=1.0d-2
  real(8),parameter::tol=1.0d-32
  real(8),parameter::diff_coef=1.0d2
  real(8),parameter::sigma=8.0d0
end module params

program main
  use params
  implicit none
  real(8),dimension(0:imax+1,0:jmax+1,0:kmax+1)::a,anew,r,rnew,p,pnew,x,xnew,ax,ap
  real(8)::alpha,beta,coef1,coef2
  real(8)::r2,rnew2,pap,b2,eps
  real(8)::t=0.0d0
  real(8)::dclock,t0,time
  real(8)::sevenstencil
  integer::i,j,k,iter,tstep

  write(6,*) "dt/dx/dx:",dt/dx/dx
  b2 = 0.0d0

  coef1 = -1.0d0*dt/dx/dx*diff_coef
  coef2  = 1.0d0+6.0d0*dt/dx/dx*diff_coef
  !  write(6,*) "coef1,coef2:",coef1,coef2

  ! initial
  !$omp parallel do
  do k=0,kmax+1
     do j=0,jmax+1
        do i=0,imax+1
           !           a(i,j,k)    = cos(dble(i))+sin(dble(j))+cos(dble(k))
           a(i,j,k)    = exp(-1.0d0*((i-imax/2)**2+(j-jmax/2)**2+(k-kmax/2)**2)/2/sigma/sigma)/sqrt(2*3.14159d0*sigma) ! exp(-r^2)
        end do
     end do
  end do

  k=kmax/2
  do j=1,jmax
     do i=1,imax
        write(599,*) i,j,a(i,j,k)
        if (i.eq.imax+1) write(599,*)
     end do
  end do

  t0 = dclock()
  do tstep=1,max_tstep
     if (mod(tstep,max_tstep/10).eq.0) write(6,*) "tstep:",tstep
     b2 = 0.0d0
     !$omp parallel private(i,j,k)
     !$omp do
     do k=0,kmax+1
        do j=0,jmax+1
           do i=0,imax+1
              x(i,j,k) = a(i,j,k) ! initial guess
           end do
        end do
     end do
     !$omp end do

     ! g = i+(imax+2)*j+(imax+2)*(jmax+2)*k
     !             i+1       i-1       j+1              j-1              k+1                       k-1
     ! coef1*(anew(g+1)+anew(g-1)+anew(g+(imax+2))+anew(g-(imax+2))+anew(g+(imax+2)*(jmax+2))+anew(g-(imax+2)*(jmax+2)))+coef2*anew(g) = a(g)
     !

     ! boundary value of ax is zero, how to handle it?
     !$omp do
     do k=1,kmax
        do j=1,jmax
           do i=1,imax
              ax(i,j,k) = coef1*(x(i+1,j,  k  )+x(i-1,j,  k  )  &
                                +x(i,  j+1,k  )+x(i,  j-1,k  )  &
                                +x(i,  j,  k+1)+x(i,  j,  k-1)) &
                          +coef2*x(i,  j,  k  )
           end do
        end do
     end do
     !$omp end do
     !$omp do reduction(+:b2)
     do k=1,kmax
        do j=1,jmax
           do i=1,imax
              r(i,j,k) = a(i,j,k)-ax(i,j,k)
              p(i,j,k) = r(i,j,k)
              b2       = b2 + a(i,j,k)*a(i,j,k)
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel
     ! exclude t loop
     do iter=1,max_iter
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
                 ap(i,j,k) = coef1*(p(i+1,j,  k  )+p(i-1,j,  k  )  &
                                   +p(i,  j+1,k  )+p(i,  j-1,k  )  &
                                   +p(i,  j,  k+1)+p(i,  j,  k-1)) &
                             +coef2*p(i,  j,  k  )
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
        else if (iter.eq.max_iter.and.eps.ge.tol) then
           write(6,*) "did not converge. residual,tol:",eps,tol
           stop
        end if
        beta = rnew2/r2
        if (isnan(alpha).or.isnan(beta)) then
           write(6,*) "alpha,beta:",alpha,beta
           stop
        end if

        ! if (mod(iter,max_iter/(imax+2)/(jmax+2)).eq.0) then
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
  end do ! tstep

  write(777) a
  time = dclock()-t0
  write(6,*) "time[s]:",time

  k=kmax/2
  do j=1,jmax
     do i=1,imax
        write(600,*) i,j,a(i,j,k)
        if (i.eq.imax+1) write(600,*)
     end do
  end do

  stop
end program main
