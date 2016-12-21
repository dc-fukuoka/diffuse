! cannot find a good way to write multiple statements in fortran macro function...
#define myassert(ret) if (ret.ne.0) write(6,'(2a,i0,a,i0)') __FILE__, ": ", __LINE__, ": ret: ", ret

module params
  use mpi
  implicit none
  integer,parameter::ndims=3
  character(len=32)::data_in = "data_in"
  real(8)::dt,dx
  integer::imax, jmax, kmax ! global
  integer::imax_l, jmax_l, kmax_l
  integer::idiv, jdiv, kdiv ! # of division for each dimensions
  integer::ierr
end module params

! to be able to pass allocatable variables, need to use module...
module mysubs
contains
  subroutine myinit(iam,np)
    use params
    implicit none
    integer,intent(out)::iam,np

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
    real(8),allocatable,dimension(:,:,:)::array
    if(allocated(array)) deallocate(array)
    return
  end subroutine free_data

  subroutine myfini(a_l,anew_l,type_j,type_k,ifiletype)
    use params
    implicit none
    real(8),allocatable,dimension(:,:,:),intent(inout)::a_l,anew_l
    integer,intent(inout)::type_j,type_k,ifiletype

    call free_data(a_l)
    call free_data(anew_l)
    call free_type(type_j)
    call free_type(type_k)
    call free_type(ifiletype)
    call mpi_finalize(ierr)
    return 
  end subroutine myfini

  subroutine read_inputs(unit,iam,np)
    use params
    implicit none
    integer,intent(in)::unit,iam,np
    namelist /param0/dt,dx,imax,jmax,kmax,idiv,jdiv,kdiv
    read(unit, param0)
    if (mod(imax,idiv).ne.0.or.mod(jmax,jdiv).ne.0.or.mod(kmax,kdiv).ne.0) then
       if (iam.eq.0) write(6,*) "imax/idiv, jmax/jdiv, kmax/kdiv must be able to be devided."
       call mpi_finalize(ierr)
    else if (np.ne.idiv*jdiv*kdiv) then
       if (iam.eq.0) write(6,*) "np must equal to idiv*jdiv*kdiv"
       call mpi_finalize(ierr)
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

  subroutine create_datatypes(type_j,type_k,ifiletype, &
       src_i,dest_i,src_j,dest_j,src_k,dest_k, &
       if_update_i,if_update_j,if_update_k)
    use params
    implicit none
    integer,intent(out)::type_j,type_k,ifiletype
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

    call mpi_type_vector(jmax_l,1,imax_l+2,             mpi_real8,type_j,ierr)
    call mpi_type_vector(kmax_l,1,(imax_l+2)*(jmax_l+2),mpi_real8,type_k,ierr)
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
    call mpi_type_commit(type_j,ierr)
    call mpi_type_commit(type_k,ierr)
    call mpi_type_commit(ifiletype,ierr)

    call mpi_cart_shift(comm_cart,0,1,src_i,dest_i,ierr)
    call mpi_cart_shift(comm_cart,1,1,src_j,dest_j,ierr)
    call mpi_cart_shift(comm_cart,2,1,src_k,dest_k,ierr)

    ! there is no need to send/recv on the boundaries
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

#ifdef _DEBUG
    ! read global
    if (iam.eq.0) then
       write(6,*) "imax,jmax,kmax:",imax,jmax,kmax
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
             write(100+iam,"(a,3i2,a,f6.2)") "a_l(",i,j,k,"): ",a_l(i,j,k)
          end do
       end do
    end do
    ! ok
#endif ! _DEBUG

    return
  end subroutine read_initial_data

  subroutine diffuse(a_l,anew_l,type_j,type_k,ifiletype, &
       src_i,dest_i,src_j,dest_j,src_k,dest_k, &
       if_update_i,if_update_j,if_update_k)
    use params
    implicit none
    real(8),dimension(0:imax_l+1,0:jmax_l+1,0:kmax_l+1),intent(inout)::a_l,anew_l
    integer,intent(in)::type_j,type_k,ifiletype
    integer,intent(in)::src_i,dest_i,src_j,dest_j,src_k,dest_k
    logical,dimension(2),intent(in)::if_update_i,if_update_j,if_update_k

    real(8)::alpha
    real(8)::beta

    alpha= -1.0d0*dt/dx/dx
    beta = 6*dt/dx/dx+1
    ! g = i+(imax+2)*j+(imax+2)*(jmax+2)*k (fortran order)
    ! alpa*a^(n+1)(g+1)+alpa*a^(n+1)(g-1)+alpa*a^(n+1)(g+imax+2)+alpa*a^(n+1)(g-(imax+2))+alpa*a^(n+1)(g+(imax+2)*(jmax+2))+alpa*a^(n+1)(g-(imax+2)*(jmax+2))+beta*alpa*a^(n+1)(g)
    !  = a^(n)(g)

    ! petsc here
    
    return
  end subroutine diffuse
end module mysubs

program main
  use mpi
  use params
  use mysubs
  implicit none
  integer::iam,np
  real(8),allocatable,dimension(:,:,:)::a_g ! debug
  real(8),allocatable,dimension(:,:,:)::a_l, anew_l
  integer::type_j,type_k,ifiletype
  integer::src_i,dest_i,src_j,dest_j,src_k,dest_k
  logical,dimension(2)::if_update_i=.true.,if_update_j=.true.,if_update_k=.true.

  call myinit(iam,np)
  call read_inputs(11,iam,np)
  call allocate_arrays(a_l,anew_l)
  call create_datatypes(type_j,type_k,ifiletype, &
       src_i,dest_i,src_j,dest_j,src_k,dest_k, &
       if_update_i,if_update_j,if_update_k)
  call read_initial_data(a_l,a_g,ifiletype,iam,np)
  call diffuse(a_l,anew_l,type_j,type_k,ifiletype, &
       src_i,dest_i,src_j,dest_j,src_k,dest_k, &
       if_update_i,if_update_j,if_update_k)
  call myfini(a_l,anew_l,type_j,type_k,ifiletype)
  stop
end program main
