! cannot find a good way to write multiple statements in fortran macro function...
#define myassert(ret) if (ret.ne.0) write(6,'(2a,i0,a,i0)') __FILE__, ": ", __LINE__, ": ret: ", ret

module params
  use mpi
  implicit none
  integer,parameter::ndims=3
  character(len=32)::data_in = "data_in"
  real(8)::dt,dx,tol
  integer::imax, jmax, kmax ! global
  integer::imax_l, jmax_l, kmax_l
  integer::idiv, jdiv, kdiv ! # of division for each dimensions
  integer::iter_max ! maximum iteration for CG method
  integer::iam,np,ierr
end module params

! to be able to pass allocatable variables, need to use module...
module mysubs
contains
  subroutine myinit()
    use params
    implicit none

    call mpi_init(ierr)
#ifdef _DEBUG
    myassert(ierr) ! test
#endif ! _DEBUG
    call mpi_comm_rank(mpi_comm_world,iam,ierr)
    call mpi_comm_size(mpi_comm_world,np ,ierr)
    return
  end subroutine myinit

  subroutine free_type(type)
    use params
    implicit none
    integer,intent(inout)::type
    if (type.ne.0) call mpi_type_free(type,ierr)
    return
  end subroutine free_type

  subroutine free_data(array)
    use params
    implicit none
    real(8),allocatable,dimension(:,:,:),intent(inout)::array
    if(allocated(array)) deallocate(array)
    return
  end subroutine free_data

  subroutine myfini(a_l,anew_l,ifiletype)
    use params
    implicit none
    real(8),allocatable,dimension(:,:,:),intent(inout)::a_l,anew_l
    integer,intent(inout)::ifiletype

    call free_data(a_l)
    call free_data(anew_l)
    call free_type(ifiletype)
    call mpi_finalize(ierr)
    return 
  end subroutine myfini

  subroutine read_inputs(unit)
    use params
    implicit none
    integer,intent(in)::unit
    namelist /param0/dt,dx,tol
    namelist /param1/imax,jmax,kmax,iter_max
    namelist /param2/idiv,jdiv,kdiv
    read(unit, param0)
    read(unit, param1)
    read(unit, param2)
    if (mod(imax,idiv).ne.0.or.mod(jmax,jdiv).ne.0.or.mod(kmax,kdiv).ne.0) then
       if (iam.eq.0) write(6,*) "imax/idiv, jmax/jdiv, kmax/kdiv must be able to be devided."
       call mpi_finalize(ierr)
    else if (np.ne.idiv*jdiv*kdiv) then
       if (iam.eq.0) write(6,*) "np must equal to idiv*jdiv*kdiv"
       call mpi_finalize(ierr)
    end if
    if (iam.eq.0) then
       write(6,"(a,3(1pe14.5))") "dt,dx,tol:",dt,dx,tol
       write(6,"(a,4i5)") "imax,jmax,kmax,iter_max:",imax,jmax,kmax,iter_max
       write(6,"(a,3i3)") "idiv,jdiv,kdiv:",idiv,jdiv,kdiv
    end if
    imax_l = imax/idiv
    jmax_l = jmax/jdiv
    kmax_l = kmax/kdiv
    return
  end subroutine read_inputs

  subroutine allocate_arrays(a_l,anew_l)
    use params
    implicit none
    real(8),allocatable,dimension(:,:,:),intent(out)::a_l,anew_l
    integer::i,j,k

    allocate(a_l(0:imax_l+1,0:jmax_l+1,0:kmax_l+1))
    allocate(anew_l(0:imax_l+1,0:jmax_l+1,0:kmax_l+1))
!$omp parallel do
    do k=0,kmax_l+1
       do j=0,jmax_l+1
          do i=0,imax_l+1
             a_l(i,j,k) = 0.0d0
             anew_l(i,j,k) = 0.0d0
          end do
       end do
    end do
    return
  end subroutine allocate_arrays

  subroutine create_datatypes(ifiletype, &
       src_i,dest_i,src_j,dest_j,src_k,dest_k, &
       if_update_i,if_update_j,if_update_k)
    use params
    implicit none
    integer,intent(out)::ifiletype
    integer,intent(out)::src_i,dest_i,src_j,dest_j,src_k,dest_k
    logical,dimension(2),intent(inout)::if_update_i,if_update_j,if_update_k ! 1: + direction, 2: - direction
    integer::comm_cart,iam_cart
    integer::idivs(ndims),coords(ndims)
    logical::prds(ndims)
    integer::i,j,k
    integer(kind=mpi_offset_kind)::idisp=0
    integer::sizes(ndims),subsizes(ndims),starts(ndims),oder,subdatatype

    prds(:) = .false.
    idivs(1) = idiv
    idivs(2) = jdiv
    idivs(3) = kdiv

    call mpi_cart_create(mpi_comm_world,ndims,idivs,prds,.false.,comm_cart,ierr)
    call mpi_comm_rank(comm_cart,iam_cart,ierr)
    call mpi_cart_coords(comm_cart,iam_cart,ndims,coords,ierr)

    ! for read the input data
    sizes(1)    = imax+2
    sizes(2)    = jmax+2
    sizes(3)    = kmax+2
    subsizes(1) = imax_l+2
    subsizes(2) = jmax_l+2
    subsizes(3) = kmax_l+2
    ! coords(*) by MPI is C order, start from 0, not 1
    starts(1)   = coords(3)*imax_l
    starts(2)   = coords(2)*jmax_l
    starts(3)   = coords(1)*kmax_l
    call mpi_type_create_subarray(ndims,sizes,subsizes,starts,mpi_order_fortran,mpi_real8,ifiletype,ierr)
    call mpi_type_commit(ifiletype,ierr)

    ! do not define datatypes for halo exchange, use temporary buffers
    
    call mpi_cart_shift(comm_cart,2,1,src_i,dest_i,ierr)
    call mpi_cart_shift(comm_cart,1,1,src_j,dest_j,ierr)
    call mpi_cart_shift(comm_cart,0,1,src_k,dest_k,ierr)

    ! there is no need to send/recv on the boundaries
    ! this was useful for OpenACC/MIC offload.
    if (coords(3).eq.0)      if_update_i(1) = .false.
    if (coords(3).eq.idiv-1) if_update_i(2) = .false.
    if (coords(2).eq.0)      if_update_j(1) = .false.
    if (coords(2).eq.jdiv-1) if_update_j(2) = .false.
    if (coords(1).eq.0)      if_update_k(1) = .false.
    if (coords(1).eq.kdiv-1) if_update_k(2) = .false.

    return
  end subroutine create_datatypes

  subroutine read_initial_data(a_l,a_g,ifiletype,iam,np)
    use params
    implicit none
    real(8),dimension(0:imax_l+1,0:jmax_l+1,0:kmax_l+1),intent(out)::a_l
    real(8),dimension(:,:,:),allocatable,intent(out)::a_g
    integer,intent(in)::ifiletype,iam,np
    integer::fh,info
    integer(kind=mpi_offset_kind)::idisp=0
    integer::istat(mpi_status_size)
    integer::i,j,k

    call mpi_info_create(info,ierr)
    call mpi_info_set(info,"striping_factor","8",ierr)
    call mpi_file_open(mpi_comm_world,data_in,mpi_mode_rdonly,info,fh,ierr)
    call mpi_file_set_view(fh,idisp,mpi_real8,ifiletype,"native",info,ierr)
    call mpi_file_read_all(fh,a_l(0,0,0),(imax_l+2)*(jmax_l+2)*(kmax_l+2),mpi_real8,istat,ierr)
    call mpi_file_close(fh,ierr)

#if 0
    ! read global
    if (iam.eq.0) then
       allocate(a_g(0:imax+1,0:jmax+1,0:kmax+1))
       open(999,file="data_in",form="unformatted",access="stream")
       do k=0,kmax+1
          do j=0,jmax+1
             do i=0,imax+1
                read(999) a_g(i,j,k)
                write(1000,"(a,3i2,a,f6.2)") "a_g(",i,j,k,"): ",a_g(i,j,k)
             end do
          end do
       end do
       close(999)
       deallocate(a_g)
    end if
    ! read locals
    do k=0,kmax_l+1
       do j=0,kmax_l+1
          do i=0,imax_l+1
             write(200+iam,"(a,3i2,a,f6.2)") "a_l(",i,j,k,"): ",a_l(i,j,k)
          end do
       end do
    end do
    ! ok
#endif

    return
  end subroutine read_initial_data

  ! incomplete
  subroutine exchange_halo(a_in, &
       src_i,dest_i,src_j,dest_j,src_k,dest_k, &
       if_update_i,if_update_j,if_update_k)
    use params
    implicit none
    real(8),dimension(0:imax+1,0:jmax+1,0:kmax+1),intent(inout)::a_in
    integer,intent(in)::src_i,dest_i,src_j,dest_j,src_k,dest_k
    logical,dimension(2),intent(in)::if_update_i,if_update_j,if_update_k
    real(8),dimension(:,:,:),allocatable::buf_i,buf_j,buf_k
    integer::ireqs(ndims*4),istats(mpi_status_size,ndims*4)
#ifdef _DEBUG
    real(8):: debugval = -1.0d0
#endif

    ! temporary buffer for halo exchange
    allocate(buf_i(jmax_l,kmax_l,4)) ! j-k plane
    allocate(buf_j(imax_l,kmax_l,4)) ! i-k plane
    allocate(buf_k(imax_l,jmax_l,4)) ! i-j plane
    ! need to be a face, not a line
    ! i direction
    buf_i(1:jmax_l,1:kmax_l,1) = a_in(1,     1:jmax_l,1:kmax_l)
    buf_i(1:jmax_l,1:kmax_l,2) = debugval
    buf_i(1:jmax_l,1:kmax_l,3) = a_in(imax_l,1:jmax_l,1:kmax_l)
    buf_i(1:jmax_l,1:kmax_l,4) = debugval
    call mpi_isend(buf_i(1,1,1),jmax_l*kmax_l,mpi_real8,dest_i,0,mpi_comm_world,ireqs(1), ierr)
    call mpi_irecv(buf_i(1,1,2),jmax_l*kmax_l,mpi_real8,src_i, 0,mpi_comm_world,ireqs(2), ierr)
    call mpi_isend(buf_i(1,1,3),jmax_l*kmax_l,mpi_real8,dest_i,1,mpi_comm_world,ireqs(3), ierr)
    call mpi_irecv(buf_i(1,1,4),jmax_l*kmax_l,mpi_real8,src_i, 1,mpi_comm_world,ireqs(4), ierr)

    ! j direction
    buf_j(1:imax_l,1:kmax_l,1) = a_in(1:imax_l,1,     1:kmax_l)
    buf_j(1:imax_l,1:kmax_l,2) = debugval
    buf_j(1:imax_l,1:kmax_l,3) = a_in(1:imax_l,jmax_l,1:kmax_l)
    buf_j(1:imax_l,1:kmax_l,4) = debugval
    call mpi_isend(buf_j(1,1,1),imax_l*kmax_l,mpi_real8,dest_j,2,mpi_comm_world,ireqs(5), ierr)
    call mpi_irecv(buf_j(1,1,2),imax_l*kmax_l,mpi_real8,src_j, 2,mpi_comm_world,ireqs(6), ierr)
    call mpi_isend(buf_j(1,1,3),imax_l*kmax_l,mpi_real8,dest_j,3,mpi_comm_world,ireqs(7), ierr)
    call mpi_irecv(buf_j(1,1,4),imax_l*kmax_l,mpi_real8,src_j, 3,mpi_comm_world,ireqs(8), ierr)

    ! k direction
    buf_k(1:imax_l,1:jmax_l,1) = a_in(1:imax_l,1:jmax_l,1     )
    buf_k(1:imax_l,1:jmax_l,2) = debugval
    buf_k(1:imax_l,1:jmax_l,3) = a_in(1:imax_l,1:jmax_l,kmax_l)
    buf_k(1:imax_l,1:jmax_l,4) = debugval
    call mpi_isend(buf_k(1,1,1),imax_l*jmax_l,mpi_real8,dest_k,4,mpi_comm_world,ireqs(9), ierr)
    call mpi_irecv(buf_k(1,1,2),imax_l*jmax_l,mpi_real8,src_k, 4,mpi_comm_world,ireqs(10),ierr)
    call mpi_isend(buf_k(1,1,3),imax_l*jmax_l,mpi_real8,dest_k,5,mpi_comm_world,ireqs(11),ierr)
    call mpi_irecv(buf_k(1,1,4),imax_l*jmax_l,mpi_real8,src_k, 5,mpi_comm_world,ireqs(12),ierr)
    call mpi_waitall(ndims*4,ireqs,istats,ierr)

    ! i direction
    if (if_update_i(1)) a_in(imax_l+1,1:jmax_l,1:kmax_l) = buf_i(1:jmax_l,1:kmax_l,2)
    if (if_update_i(2)) a_in(0,       1:jmax_l,1:kmax_l) = buf_i(1:jmax_l,1:kmax_l,4)
    ! j drection
    if (if_update_j(1)) a_in(1:imax_l,jmax_l+1,1:kmax_l) = buf_j(1:imax_l,1:kmax_l,2)
    if (if_update_j(2)) a_in(1:imax_l,0,       1:kmax_l) = buf_j(1:imax_l,1:kmax_l,4)
    ! k direction
    if (if_update_k(1)) a_in(1:imax_l,1:jmax_l,kmax_l+1) = buf_k(1:imax_l,1:jmax_l,2)
    if (if_update_k(2)) a_in(1:imax_l,1:jmax_l,0       ) = buf_k(1:imax_l,1:jmax_l,4)
    
    deallocate(buf_i,buf_j,buf_k)
    return
  end subroutine exchange_halo

  subroutine diffuse(a_l,anew_l,ifiletype, &
       src_i,dest_i,src_j,dest_j,src_k,dest_k, &
       if_update_i,if_update_j,if_update_k)
    use params
    implicit none
    real(8),dimension(0:imax_l+1,0:jmax_l+1,0:kmax_l+1),intent(inout)::a_l,anew_l
    real(8),dimension(:,:,:),allocatable::r_l,rnew_l,p_l,pnew_l,x_l,xnew_l
    integer,intent(in)::ifiletype
    integer,intent(in)::src_i,dest_i,src_j,dest_j,src_k,dest_k
    logical,dimension(2),intent(in)::if_update_i,if_update_j,if_update_k
    real(8)::coef1,coef2,alpha,beta
    real(8)::r2,pap,norm2,r2_l,pap_l,norm2_l
    integer::i,j,k,tstep,iter

    coef1 = -1.0d0*dt/dx/dx
    coef2 = 6.0d0*dt/dx/dx+1

    allocate(r_l(0:imax+1,0:jmax+1,0:kmax+1),rnew_l(0:imax+1,0:jmax+1,0:kmax+1))
    allocate(p_l(0:imax+1,0:jmax+1,0:kmax+1),pnew_l(0:imax+1,0:jmax+1,0:kmax+1))
    allocate(x_l(0:imax+1,0:jmax+1,0:kmax+1),xnew_l(0:imax+1,0:jmax+1,0:kmax+1))

    ! g = i+(imax+2)*j+(imax+2)*(jmax+2)*k (fortran order)
    ! coef1*(anew(g+1)+anew(g-1)+anew(g+imax+2)+anew(g-(imax+2))+anew(g+(imax+2)*(jmax+2))+anew(g-(imax+2)*(jmax+2)))+coef2anew(g) = a(g)
    
    ! CG method

    do k=0,kmax_l+1
       do j=0,jmax_l+1
          do i=0,imax_l+1
             x_l(i,j,k) = a_l(i,j,k) ! initial guess
          end do
       end do
    end do

    do k=0,kmax_l+1
       do j=0,jmax_l+1
          do i=0,imax_l+1
             r_l(i,j,k) = a_l(i,j,k) - (coef1*(x_l(i+1,j,k)+x_l(i-1,j,k)+x_l(i,j+1,k)+x_l(i,j-1,k)+x_l(i,j,k+1)+x_l(i,j,k-1))+coef2*x_l(i,j,k))
             p_l(i,j,k) = r_l(i,j,k)
          end do
       end do
    end do

    do tstep=1,1 ! time step, test
       do iter=1,iter_max ! CG method iteration
          norm2   = 0.0d0
          r2      = 0.0d0
          pap     = 0.0d0
          norm2_l = 0.0d0
          r2_l    = 0.0d0
          pap_l   = 0.0d0
          do k=0,kmax_l+1
             do j=0,jmax_l+1
                do i=0,imax_l+1
                   r2_l = r2_l + r_l(i,j,k)*r_l(i,j,k)
                   pap_l = pap_l + p_l(i,j,k)*(coef1*(p_l(i+1,j,k)+p_l(i-1,j,k)+p_l(i,j+1,k)+p_l(i,j,k+1)+p_l(i,j,k-1))+coef2*p_l(i,j,k))
                end do
             end do
          end do
          
          call mpi_allreduce(r2_l, r2, 1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
          call mpi_allreduce(pap_l,pap,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
          
          alpha = r2/pap
          do k=0,kmax_l+1
             do j=0,jmax_l+1
                do i=0,imax_l+1
                   xnew_l(i,j,k) = x_l(i,j,k)+alpha*p_l(i,j,k)
                   rnew_l(i,j,k) = r_l(i,j,k)-alpha*(coef1*(p_l(i+1,j,k)+p_l(i-1,j,k)+x_l(i,j+1,k)+p_l(i,j,k+1)+p_l(i,j,k-1))+coef2*p_l(i,j,k))
                   norm2_l       = norm2_l + rnew_l(i,j,k)*rnew_l(i,j,k)
                end do
             end do
          end do
          
          call mpi_allreduce(norm2_l,norm2,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
          
#ifdef _DEBUG
          if (iam.eq.0) then
             if (mod(iter,imax_l*jmax_l).eq.0) write(6,"(a,i5,2(1pe14.5))") "iter,residual,tol:",iter,sqrt(norm2),tol
          end if
#endif
          if (sqrt(norm2).le.tol) then
             write(6,"(a,i5,2(1pe14.5))") "iter,residual,tol:",iter,sqrt(norm2),tol
             do k=0,kmax_l+1
                do j=0,jmax_l+1
                   do i=0,imax_l+1
                      anew_l(i,j,k) = xnew_l(i,j,k)
                   end do
                end do
             end do
             exit
          else if (iter.eq.iter_max.and.sqrt(norm2).le.tol) then
             if (iam.eq.0) then
                write(6,"(a,2(1pe14.5))") "the result was not converged, residual, tolerance:",sqrt(norm2),tol
                call mpi_finalize(ierr)
                stop
             end if
          end if
          
          beta = norm2/r2
          
          do k=0,kmax+1
             do j=0,jmax+1
                do i=0,imax+1
                   pnew_l(i,j,k) = rnew_l(i,j,k)+beta*p_l(i,j,k)
                end do
             end do
          end do
          
          do k=0,kmax+1
             do j=0,jmax+1
                do i=0,imax+1
                   x_l(i,j,k) = xnew_l(i,j,k)
                   r_l(i,j,k) = rnew_l(i,j,k)
                   p_l(i,j,k) = pnew_l(i,j,k)
                end do
             end do
          end do

          ! segfault here..
          call exchange_halo(a_l,src_i,dest_i,src_j,dest_j,src_k,dest_k, &
               if_update_i,if_update_j,if_update_k)
          call exchange_halo(x_l,src_i,dest_i,src_j,dest_j,src_k,dest_k, &
               if_update_i,if_update_j,if_update_k)
          call exchange_halo(p_l,src_i,dest_i,src_j,dest_j,src_k,dest_k, &
               if_update_i,if_update_j,if_update_k)
          
          ! df: debug
          do j=1,jmax
             write(300+iam,"(16f8.2)") a_l(0,j,1:kmax_l)
          end do
          write(300+iam,*) "---"

       end do ! iter
       
       do k=0,kmax+1
          do j=0,jmax+1
             do i=0,imax+1
                a_l(i,j,k) = anew_l(i,j,k)
             end do
          end do
       end do
       
    end do ! tstep

    deallocate(r_l,rnew_l,p_l,pnew_l,x_l,xnew_l)
    return
  end subroutine diffuse
end module mysubs

program main
  use mpi
  use params
  use mysubs
  implicit none
  real(8),allocatable,dimension(:,:,:)::a_g ! debug
  real(8),allocatable,dimension(:,:,:)::a_l, anew_l
  integer::ifiletype
  integer::src_i,dest_i,src_j,dest_j,src_k,dest_k
  logical,dimension(2)::if_update_i=.true.,if_update_j=.true.,if_update_k=.true.

  call myinit()
  call read_inputs(11)
  call allocate_arrays(a_l,anew_l)
  call create_datatypes(ifiletype, &
       src_i,dest_i,src_j,dest_j,src_k,dest_k, &
       if_update_i,if_update_j,if_update_k)
  call read_initial_data(a_l,a_g,ifiletype,iam,np)
  call diffuse(a_l,anew_l,ifiletype, &
       src_i,dest_i,src_j,dest_j,src_k,dest_k, &
       if_update_i,if_update_j,if_update_k)
  call myfini(a_l,anew_l,ifiletype)
  stop
end program main
