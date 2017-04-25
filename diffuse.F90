! cannot find a good way to write multiple statements in fortran macro function...
#define myassert(ret) if (ret.ne.0) write(6,'(2a,i0,a,i0)') __FILE__, ": ", __LINE__, ": ret: ", ret

module params
  use mpi
  implicit none
  integer,parameter::ndims=3
  character(len=32)::data_in = "data_in"
  real(8)::dt,dx,tol
  real(8)::diff_coef ! diffusion coefficient
  integer::imax, jmax, kmax ! global
  integer::imax_l, jmax_l, kmax_l
  integer::idiv, jdiv, kdiv ! # of division for each dimensions
  integer::iter_max ! maximum iteration for CG method
  integer::tstep_max ! maximum # of time steps
  integer::freq_write ! frequency of writing the result
  integer::iam,np,ierr
end module params

! to be able to return allocated arrays, need to use module...
module mysubs
  contains
  subroutine myinit(fh)
    use params
    implicit none
    integer,intent(out)::fh
    integer::info
    character(8)::data_out="data_out"
    
    ! call mpi_info_set(info,"striping_factor","8",ierr)
    info = mpi_info_null
    
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,iam,ierr)
    call mpi_comm_size(mpi_comm_world,np ,ierr)
    call mpi_file_open(mpi_comm_world,data_out,mpi_mode_wronly+mpi_mode_create,info,fh,ierr)
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

  subroutine myfini(a_l,anew_l,ifiletype_read,ifiletype_write,fh,comm_cart)
    use params
    implicit none
    real(8),allocatable,dimension(:,:,:),intent(inout)::a_l,anew_l
    integer,intent(inout)::ifiletype_read,fh,ifiletype_write,comm_cart
    character(8)::data_out="data_out"

    call free_data(a_l)
    call free_data(anew_l)
    call free_type(ifiletype_read)
    call free_type(ifiletype_write)
    call mpi_comm_free(comm_cart,ierr)
    call mpi_file_close(fh,ierr)
    call mpi_finalize(ierr)
    return
  end subroutine myfini

  subroutine read_inputs(unit)
    use params
    implicit none
    integer,intent(in)::unit
    namelist /param0/dt,dx,tol,diff_coef
    namelist /param1/imax,jmax,kmax
    namelist /param2/idiv,jdiv,kdiv
    namelist /param3/iter_max,tstep_max,freq_write
    read(unit, param0)
    read(unit, param1)
    read(unit, param2)
    read(unit, param3)
    if (mod(imax,idiv).ne.0.or.mod(jmax,jdiv).ne.0.or.mod(kmax,kdiv).ne.0) then
       if (iam.eq.0) write(6,*) "imax/idiv, jmax/jdiv, kmax/kdiv must be able to be devided."
       call mpi_finalize(ierr)
       stop
    else if (np.ne.idiv*jdiv*kdiv) then
       if (iam.eq.0) write(6,*) "np must equal to idiv*jdiv*kdiv =",idiv*jdiv*kdiv
       call mpi_finalize(ierr)
       stop
    end if
    if (iam.eq.0) then
       write(6,"(a,4(1pe14.5))") "dt,dx,tol,diff_coef:",dt,dx,tol,diff_coef
       write(6,"(a,4i5)") "imax,jmax,kmax:",imax,jmax,kmax
       write(6,"(a,2i8)") "iter_max,tstep_max:",iter_max,tstep_max
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
    do k=0,kmax_l+1
       do j=0,jmax_l+1
          do i=0,imax_l+1
             a_l(i,j,k)    = 0.0d0
             anew_l(i,j,k) = 0.0d0
          end do
       end do
    end do
    return
  end subroutine allocate_arrays

  subroutine create_datatypes(ifiletype_read,ifiletype_write, &
       src_i,dest_i,src_j,dest_j,src_k,dest_k, &
       if_update_i,if_update_j,if_update_k,comm_cart)
    use params
    implicit none
    integer,intent(out)::ifiletype_read,ifiletype_write
    integer,intent(out)::src_i,dest_i,src_j,dest_j,src_k,dest_k
    logical,dimension(2),intent(inout)::if_update_i,if_update_j,if_update_k ! 1: + direction, 2: - direction
    integer,intent(out)::comm_cart
    integer::iam_cart
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
    ! includes halo
    sizes(1)    = imax+2
    sizes(2)    = jmax+2
    sizes(3)    = kmax+2
    subsizes(1) = imax_l+2
    subsizes(2) = jmax_l+2
    subsizes(3) = kmax_l+2
    ! coords(*) by MPI_Cart_coords() is C order, (0,0,0)->(0,0,1), ...
    starts(1)   = coords(1)*imax_l
    starts(2)   = coords(2)*jmax_l
    starts(3)   = coords(3)*kmax_l
    call mpi_type_create_subarray(ndims,sizes,subsizes,starts,mpi_order_fortran,mpi_real8,ifiletype_read,ierr)
    call mpi_type_commit(ifiletype_read,ierr)

    ! for write the output data
    ! excludes halo
    sizes(1)    = imax
    sizes(2)    = jmax
    sizes(3)    = kmax
    subsizes(1) = imax_l
    subsizes(2) = jmax_l
    subsizes(3) = kmax_l
    starts(1)   = coords(1)*imax_l
    starts(2)   = coords(2)*jmax_l
    starts(3)   = coords(3)*kmax_l
    call mpi_type_create_subarray(ndims,sizes,subsizes,starts,mpi_order_fortran,mpi_real8,ifiletype_write,ierr)
    call mpi_type_commit(ifiletype_write,ierr)

    ! do not define derived datatypes for halo exchange, use temporary buffers

    call mpi_cart_shift(comm_cart,0,1,src_i,dest_i,ierr)
    call mpi_cart_shift(comm_cart,1,1,src_j,dest_j,ierr)
    call mpi_cart_shift(comm_cart,2,1,src_k,dest_k,ierr)

    ! there is no need to exchange on the boundaries for non periodic boundary
    if_update_i(:) = .true.
    if_update_j(:) = .true.
    if_update_k(:) = .true.
    if (coords(1).eq.0)      if_update_i(2) = .false. ! - boundary on i direction
    if (coords(1).eq.idiv-1) if_update_i(1) = .false. ! + boundary on i direction
    if (coords(2).eq.0)      if_update_j(2) = .false. ! - boundary on j direction
    if (coords(2).eq.jdiv-1) if_update_j(1) = .false. ! + boundary on j direction
    if (coords(3).eq.0)      if_update_k(2) = .false. ! - boundary on k direction
    if (coords(3).eq.kdiv-1) if_update_k(1) = .false. ! + boundary on k direction

    return
  end subroutine create_datatypes

  subroutine read_initial_data(a_l,ifiletype_read)
    use params
    implicit none
    real(8),dimension(0:imax_l+1,0:jmax_l+1,0:kmax_l+1),intent(out)::a_l
    integer,intent(in)::ifiletype_read
    integer::fh
    integer(kind=mpi_offset_kind)::idisp=0
    integer::istat(mpi_status_size)
    integer::i,j,k

    call mpi_file_open(mpi_comm_world,data_in,mpi_mode_rdonly,mpi_info_null,fh,ierr)
    call mpi_file_set_view(fh,idisp,mpi_real8,ifiletype_read,"native",mpi_info_null,ierr)
    call mpi_file_read_all(fh,a_l(0,0,0),(imax_l+2)*(jmax_l+2)*(kmax_l+2),mpi_real8,istat,ierr)
    call mpi_file_close(fh,ierr)

    return
  end subroutine read_initial_data

  subroutine exchange_halo(a_in, &
       buf_i, buf_j, buf_k, &
       src_i,dest_i,src_j,dest_j,src_k,dest_k, &
       if_update_i,if_update_j,if_update_k,comm_cart)
    use params
    implicit none
    real(8),dimension(0:imax_l+1,0:jmax_l+1,0:kmax_l+1),intent(inout)::a_in
    real(8),dimension(1:jmax_l,1:kmax_l,4),intent(inout)::buf_i
    real(8),dimension(1:imax_l,1:kmax_l,4),intent(inout)::buf_j
    real(8),dimension(1:imax_l,1:jmax_l,4),intent(inout)::buf_k
    integer,intent(in)::src_i,dest_i,src_j,dest_j,src_k,dest_k
    logical,dimension(2),intent(in)::if_update_i,if_update_j,if_update_k
    integer,intent(in)::comm_cart
    integer::ireqs(ndims*4),istats(mpi_status_size,ndims*4)
    integer::i,j,k

    ! need to be a face, not a line
    ! i direction
    buf_i(1:jmax_l,1:kmax_l,1) = a_in(1,     1:jmax_l,1:kmax_l) ! send to -i, src_i
    buf_i(1:jmax_l,1:kmax_l,2) = 0.0d0
    buf_i(1:jmax_l,1:kmax_l,3) = a_in(imax_l,1:jmax_l,1:kmax_l) ! send to +i, dest_i
    buf_i(1:jmax_l,1:kmax_l,4) = 0.0d0
    call mpi_isend(buf_i(1,1,1),jmax_l*kmax_l,mpi_real8,src_i, 0,comm_cart,ireqs(1), ierr)
    call mpi_irecv(buf_i(1,1,2),jmax_l*kmax_l,mpi_real8,dest_i,0,comm_cart,ireqs(2), ierr)
    call mpi_isend(buf_i(1,1,3),jmax_l*kmax_l,mpi_real8,dest_i,1,comm_cart,ireqs(3), ierr)
    call mpi_irecv(buf_i(1,1,4),jmax_l*kmax_l,mpi_real8,src_i, 1,comm_cart,ireqs(4), ierr)

    ! j direction
    buf_j(1:imax_l,1:kmax_l,1) = a_in(1:imax_l,1,     1:kmax_l) ! send to -j, src_j
    buf_j(1:imax_l,1:kmax_l,2) = 0.0d0
    buf_j(1:imax_l,1:kmax_l,3) = a_in(1:imax_l,jmax_l,1:kmax_l) ! send to +j, dest_j
    buf_j(1:imax_l,1:kmax_l,4) = 0.0d0
    call mpi_isend(buf_j(1,1,1),imax_l*kmax_l,mpi_real8,src_j, 2,comm_cart,ireqs(5), ierr)
    call mpi_irecv(buf_j(1,1,2),imax_l*kmax_l,mpi_real8,dest_j,2,comm_cart,ireqs(6), ierr)
    call mpi_isend(buf_j(1,1,3),imax_l*kmax_l,mpi_real8,dest_j,3,comm_cart,ireqs(7), ierr)
    call mpi_irecv(buf_j(1,1,4),imax_l*kmax_l,mpi_real8,src_j, 3,comm_cart,ireqs(8), ierr)

    ! k direction
    buf_k(1:imax_l,1:jmax_l,1) = a_in(1:imax_l,1:jmax_l,1     ) ! send to -k, src_k
    buf_k(1:imax_l,1:jmax_l,2) = 0.0d0
    buf_k(1:imax_l,1:jmax_l,3) = a_in(1:imax_l,1:jmax_l,kmax_l) ! send to +k, dest_k
    buf_k(1:imax_l,1:jmax_l,4) = 0.0d0
    call mpi_isend(buf_k(1,1,1),imax_l*jmax_l,mpi_real8,src_k, 4,comm_cart,ireqs(9), ierr)
    call mpi_irecv(buf_k(1,1,2),imax_l*jmax_l,mpi_real8,dest_k,4,comm_cart,ireqs(10),ierr)
    call mpi_isend(buf_k(1,1,3),imax_l*jmax_l,mpi_real8,dest_k,5,comm_cart,ireqs(11),ierr)
    call mpi_irecv(buf_k(1,1,4),imax_l*jmax_l,mpi_real8,src_k, 5,comm_cart,ireqs(12),ierr)
    call mpi_waitall(ndims*4,ireqs,istats,ierr)

    ! if (coords(1).eq.0)      if_update_i(2) = .false. ! - boundary on i direction
    ! if (coords(1).eq.idiv-1) if_update_i(1) = .false. ! + boundary on i direction
    ! if (coords(2).eq.0)      if_update_j(2) = .false. ! - boundary on j direction
    ! if (coords(2).eq.jdiv-1) if_update_j(1) = .false. ! + boundary on j direction
    ! if (coords(3).eq.0)      if_update_k(2) = .false. ! - boundary on k direction
    ! if (coords(3).eq.kdiv-1) if_update_k(1) = .false. ! + boundary on k direction

    ! i direction
    if (if_update_i(1)) a_in(imax_l+1,1:jmax_l,1:kmax_l) = buf_i(1:jmax_l,1:kmax_l,2) ! receive from -i direction i=1      -> imax_l+1
    if (if_update_i(2)) a_in(0,       1:jmax_l,1:kmax_l) = buf_i(1:jmax_l,1:kmax_l,4) ! receive from +i direction i=imax_l -> 0
    ! j drection
    if (if_update_j(1)) a_in(1:imax_l,jmax_l+1,1:kmax_l) = buf_j(1:imax_l,1:kmax_l,2) ! receive from -j direction j=1      -> jmax_l+1
    if (if_update_j(2)) a_in(1:imax_l,0,       1:kmax_l) = buf_j(1:imax_l,1:kmax_l,4) ! receive from +j direction j=jmax_l -> 0
    ! k direction
    if (if_update_k(1)) a_in(1:imax_l,1:jmax_l,kmax_l+1) = buf_k(1:imax_l,1:jmax_l,2) ! receive from -k direction k=1      -> kmax_l+1
    if (if_update_k(2)) a_in(1:imax_l,1:jmax_l,0       ) = buf_k(1:imax_l,1:jmax_l,4) ! receive from +k direction k=kmax_l -> 0

    return
  end subroutine exchange_halo

  subroutine diffuse(a_l,anew_l,ifiletype_read, &
       src_i,dest_i,src_j,dest_j,src_k,dest_k, &
       if_update_i,if_update_j,if_update_k,comm_cart, &
       ifiletype_write,fh)
    use params
    implicit none
    real(8),dimension(0:imax_l+1,0:jmax_l+1,0:kmax_l+1),intent(inout)::a_l,anew_l
    real(8),dimension(:,:,:),allocatable::r_l,rnew_l,p_l,pnew_l,x_l,xnew_l
    integer,intent(in)::ifiletype_read
    integer,intent(in)::src_i,dest_i,src_j,dest_j,src_k,dest_k
    integer,intent(in)::comm_cart
    integer,intent(in)::ifiletype_write
    logical,dimension(2),intent(in)::if_update_i,if_update_j,if_update_k
    real(8),dimension(:,:,:),allocatable::buf_i,buf_j,buf_k
    real(8)::coef1,coef2,alpha,beta,eps
    real(8)::r2,pap,rnew2,r2_l,pap_l,rnew2_l,b2,b2_l
    integer::i,j,k,tstep,iter
    integer(kind=mpi_offset_kind)::count_write=0
    integer,intent(in)::fh


    ! temporary buffer for halo exchange
    allocate(buf_i(jmax_l,kmax_l,4)) ! j-k plane
    allocate(buf_j(imax_l,kmax_l,4)) ! i-k plane
    allocate(buf_k(imax_l,jmax_l,4)) ! i-j plane
    allocate(r_l(0:imax_l+1,0:jmax_l+1,0:kmax_l+1),rnew_l(0:imax_l+1,0:jmax_l+1,0:kmax_l+1))
    allocate(p_l(0:imax_l+1,0:jmax_l+1,0:kmax_l+1),pnew_l(0:imax_l+1,0:jmax_l+1,0:kmax_l+1))
    allocate(x_l(0:imax_l+1,0:jmax_l+1,0:kmax_l+1),xnew_l(0:imax_l+1,0:jmax_l+1,0:kmax_l+1))

    coef1 = -1.0d0*dt/dx/dx*diff_coef
    coef2 =  1.0d0+6.0d0*dt/dx/dx*diff_coef

    do k=0,kmax_l+1
       do j=0,jmax_l+1
          do i=0,imax_l+1
             r_l(i,j,k)    = 0.0d0
             rnew_l(i,j,k) = 0.0d0
             p_l(i,j,k)    = 0.0d0
             pnew_l(i,j,k) = 0.0d0
             x_l(i,j,k)    = 0.0d0
             xnew_l(i,j,k) = 0.0d0
          end do
       end do
    end do

    ! g = i+(imax+2)*j+(imax+2)*(jmax+2)*k (fortran order)
    ! coef1*(anew(g+1)+anew(g-1)+anew(g+imax+2)+anew(g-(imax+2))+anew(g+(imax+2)*(jmax+2))+anew(g-(imax+2)*(jmax+2)))+coef2*anew(g) = a(g)
    ! 7-stencil symmetric matrix, CG method can be used

    do tstep=1,tstep_max ! time step
#if 1
       if (iam.eq.0) then
          if (mod(tstep,tstep_max/10).eq.0) then
             write(6,*) "tstep:",tstep
          end if
       end if
#endif

       call exchange_halo(a_l,buf_i,buf_j,buf_k, &
            src_i,dest_i,src_j,dest_j,src_k,dest_k, &
            if_update_i,if_update_j,if_update_k,comm_cart)

       do k=0,kmax_l+1
          do j=0,jmax_l+1
             do i=0,imax_l+1
                x_l(i,j,k) = a_l(i,j,k) ! initial guess
             end do
          end do
       end do

       b2   = 0.0d0
       b2_l = 0.0d0

       call exchange_halo(x_l,buf_i,buf_j,buf_k, &
            src_i,dest_i,src_j,dest_j,src_k,dest_k, &
            if_update_i,if_update_j,if_update_k,comm_cart)

       do k=1,kmax_l
          do j=1,jmax_l
             do i=1,imax_l
                r_l(i,j,k) = a_l(i,j,k)-(coef1*(x_l(i+1,j,  k  )+x_l(i-1,j,  k  )  &
                                               +x_l(i,  j+1,k  )+x_l(i,  j-1,k  )  &
                                               +x_l(i,  j,  k+1)+x_l(i,  j,  k-1)) &
                                         +coef2*x_l(i,  j,  k  ))
                p_l(i,j,k) = r_l(i,j,k)
                b2_l       = b2_l+a_l(i,j,k)*a_l(i,j,k)
             end do
          end do
       end do

       call mpi_allreduce(b2_l, b2, 1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

       do iter=1,iter_max ! CG method iteration
          rnew2   = 0.0d0
          r2      = 0.0d0
          pap     = 0.0d0
          r2_l    = 0.0d0
          rnew2_l = 0.0d0
          pap_l   = 0.0d0
          eps     = 0.0d0
          alpha   = 0.0d0
          beta    = 0.0d0

          ! xnew_l(:,:,:)  = 0.0d0
          ! rnew_l(:,:,:)  = 0.0d0
          ! pnew_l(:,:,:)  = 0.0d0

          call exchange_halo(p_l,buf_i,buf_j,buf_k, &
               src_i,dest_i,src_j,dest_j,src_k,dest_k, &
               if_update_i,if_update_j,if_update_k,comm_cart)

          do k=1,kmax_l
             do j=1,jmax_l
                do i=1,imax_l
                   r2_l  = r2_l+r_l(i,j,k)*r_l(i,j,k)
                   pap_l = pap_l+p_l(i,j,k)*(coef1*(p_l(i+1,j,  k  )+p_l(i-1,j,  k  )  &
                                                   +p_l(i,  j+1,k  )+p_l(i,  j-1,k  )  &
                                                   +p_l(i,  j,  k+1)+p_l(i,  j,  k-1)) &
                                             +coef2*p_l(i,  j,  k  ))
                end do
             end do
          end do

          call mpi_allreduce(r2_l, r2, 1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
          call mpi_allreduce(pap_l,pap,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
#if 0
          if (iam.eq.0) then
             if (iter.eq.1) then
                write(6,*) "coef1,coef2:",coef1,coef2
                write(6,*) "r2_l,r2:",r2_l,r2
                write(6,*) "pap_l,pap;",pap_l,pap
             end if
          end if
#endif
          alpha = r2/pap
          do k=1,kmax_l
             do j=1,jmax_l
                do i=1,imax_l
                   xnew_l(i,j,k) = x_l(i,j,k)+alpha*p_l(i,j,k)
                   rnew_l(i,j,k) = r_l(i,j,k)-alpha*(coef1*(p_l(i+1,j,  k  )+p_l(i-1,j,  k  )  &
                                                           +p_l(i,  j+1,k  )+p_l(i,  j-1,k  )  &
                                                           +p_l(i,  j,  k+1)+p_l(i,  j,  k-1)) &
                                                     +coef2*p_l(i,  j,  k  ))
                   rnew2_l       = rnew2_l+rnew_l(i,j,k)*rnew_l(i,j,k)
                end do
             end do
          end do

          call mpi_allreduce(rnew2_l,rnew2,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

          eps = sqrt(rnew2)/sqrt(b2)
#if 0
          if (iam.eq.0) then
             if (mod(iter,imax_l*jmax_l).eq.0) then
                write(6,*) "rnew2_l,rnew2:",rnew2_l,rnew2
                write(6,"(a,i5,2(1pe14.5))") "iter,residual,tol:",iter,eps,tol
             end if
          end if
#endif
          if (eps.le.tol) then
#if 0
             if (iam.eq.0) then
                write(6,*) "the result converged."
                write(6,"(a,i5,2(1pe14.5))") "iter,residual,tol:",iter,eps,tol
             end if
#endif
             do k=1,kmax_l
                do j=1,jmax_l
                   do i=1,imax_l
                      anew_l(i,j,k) = xnew_l(i,j,k)
                   end do
                end do
             end do
             exit
          else if (iter.eq.iter_max.and.eps.ge.tol) then
             if (iam.eq.0) then
                write(6,"(a,2i7,2(1pe14.5))") "the result did not converge, tstep,iter,residual, tolerance:",tstep,iter,eps,tol
             end if
             call mpi_finalize(ierr)
             stop
          end if

          beta = rnew2/r2

          do k=1,kmax_l
             do j=1,jmax_l
                do i=1,imax_l
                   pnew_l(i,j,k) = rnew_l(i,j,k)+beta*p_l(i,j,k)
                end do
             end do
          end do

          do k=1,kmax_l
             do j=1,jmax_l
                do i=1,imax_l
                   x_l(i,j,k) = xnew_l(i,j,k)
                   r_l(i,j,k) = rnew_l(i,j,k)
                   p_l(i,j,k) = pnew_l(i,j,k)
                end do
             end do
          end do

       end do ! iter

       ! converged
       do k=1,kmax_l
          do j=1,jmax_l
             do i=1,imax_l
                a_l(i,j,k) = anew_l(i,j,k)
             end do
          end do
       end do

       if (mod(tstep,tstep_max/freq_write).eq.0) then
          call write_output(a_l,ifiletype_write,fh,count_write)
          count_write=count_write+1
       end if

    end do ! tstep

    deallocate(r_l,rnew_l,p_l,pnew_l,x_l,xnew_l)
    deallocate(buf_i,buf_j,buf_k)

    return
  end subroutine diffuse

  subroutine write_output(a_l,ifiletype_write,fh,ntimes)
    use params
    implicit none
    real(8),dimension(0:imax_l+1,0:jmax_l+1,0:kmax_l+1),intent(in)::a_l
    real(8),dimension(imax_l,jmax_l,kmax_l)::tmp
    integer,intent(in)::ifiletype_write,fh
    integer::info
    integer(kind=mpi_offset_kind)::idisp=0, ntimes
    integer::istat(mpi_status_size)
    integer::i,j,k

    idisp = ntimes*sizeof(a_l(1,1,1))*imax*jmax*kmax

    tmp(1:imax_l,1:jmax_l,1:kmax_l) = a_l(1:imax_l,1:jmax_l,1:kmax_l)

    ! C-like writing, exclude datasize information
    info = mpi_info_null
    call mpi_file_set_view(fh,idisp,mpi_real8,ifiletype_write,"native",info,ierr)
    call mpi_file_write_all(fh,tmp(1,1,1),imax_l*jmax_l*kmax_l,mpi_real8,istat,ierr)
  end subroutine write_output
end module mysubs

program main
  use mpi
  use params
  use mysubs
  implicit none
  real(8),allocatable,dimension(:,:,:)::a_l,anew_l
  integer::fh,ifiletype_read,ifiletype_write
  integer::src_i,dest_i,src_j,dest_j,src_k,dest_k
  logical,dimension(2)::if_update_i=.true.,if_update_j=.true.,if_update_k=.true.
  integer::comm_cart
  real(8)::time,t0

  call myinit(fh)
  call read_inputs(11)
  call allocate_arrays(a_l,anew_l)
  call create_datatypes(ifiletype_read,ifiletype_write, &
       src_i,dest_i,src_j,dest_j,src_k,dest_k, &
       if_update_i,if_update_j,if_update_k,comm_cart)
  call read_initial_data(a_l,ifiletype_read)
  call mpi_barrier(mpi_comm_world,ierr)
  t0 = mpi_wtime()
  call diffuse(a_l,anew_l,ifiletype_read, &
       src_i,dest_i,src_j,dest_j,src_k,dest_k, &
       if_update_i,if_update_j,if_update_k,comm_cart,ifiletype_write,fh)
  call mpi_barrier(mpi_comm_world,ierr)
  time = mpi_wtime() - t0
  if (iam.eq.0) write(6,*) "time[s]:",time
  call myfini(a_l,anew_l,ifiletype_read,ifiletype_write,fh,comm_cart)
  stop
end program main
