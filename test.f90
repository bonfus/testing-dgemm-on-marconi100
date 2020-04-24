SUBROUTINE count_cores()
    USE OMP_LIB
    implicit none
    include 'mpif.h'

    integer io, nme, ierror

    call MPI_COMM_RANK(MPI_COMM_WORLD, io, ierror)

    if(io==0)    WRITE(*,*) "Number of cores: ", OMP_GET_MAX_THREADS()
    RETURN
END SUBROUTINE

SUBROUTINE dimentionMatrix(dimen, io_key)
implicit none

    integer :: num1, num2
    integer(kind=4)::dimen, io_key
    CHARACTER(100) :: arg1char, arg2char

    IF(COMMAND_ARGUMENT_COUNT().NE.2)THEN
        WRITE(*,*)'Error. Parameters are required!. Ej. 1000 1'
        STOP
    ENDIF

    CALL GET_COMMAND_ARGUMENT(1,arg1char)
    CALL GET_COMMAND_ARGUMENT(2,arg2char)

    READ(arg1char,*) num1
    READ(arg2char,*) num2

    dimen = num1
    io_key = num2

!    WRITE(*,*) "Matrix Dimension: ", dimen, " x ", dimen

END SUBROUTINE

SUBROUTINE outputMatrix(m, n_matrix, name_file)
implicit none
    real(kind=8),allocatable::m(:,:)
    character(len = 15) :: name_file
    real::aux
    integer :: n_matrix, i, j, n
    
    n = n_matrix
    allocate(m(n,n))

    open(1, file = name_file, status='new')
    call random_seed()
    !$OMP PARALLEL DO
    do i = 1, n
        do j = 1, n
            call random_number(m(i, j))
            write(1,*)  m(i,j)
        end do
    end do
    !$OMP END PARALLEL DO
    close(1)
END SUBROUTINE

SUBROUTINE ioMatrix(m, n_matrix, name_file)
    implicit none
        real(kind=8),allocatable::m(:,:)
        character(len = 15) :: name_file
        integer :: n_matrix, i, j, n
        
        n = n_matrix
        allocate(m(n,n))
    
        open(1, file = name_file, status='new')
        call random_seed()
        !$OMP PARALLEL DO
        do i = 1, n
            do j = 1, n
                call random_number(m(i, j))
                write(1,*)  m(i,j)
            end do
        end do
        !$OMP END PARALLEL DO
        close(1)

        ! opening file for reading
        open (2, file = name_file, status = 'old')
        !$OMP PARALLEL DO
        do i = 1,n
            do j = 1, n
                read(2,*) m(i,j)
            end do 
        end do 
        !$OMP END PARALLEL DO
        close(2)
END SUBROUTINE

    
PROGRAM matrixMultiplication
implicit none   

include 'mpif.h'

    integer(kind=4)::n=10, key=0 ! 10 by default
    real(kind=8),allocatable::a(:,:),b(:,:),c(:,:),d(:,:),e(:,:)
    real(kind=8)::alpha,beta
    integer(kind=4)::i,j,k,lda,ldb,lde
    real(kind=4)::start, finish, var_triple, var_matmul, var_dgemm
    real(kind=4)::var_triple_sum, var_matmul_sum, var_dgemm_sum

    integer io, nme, ierror

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nme, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, io, ierror)

    if(io==0)    WRITE(*,*)'****Parallel Matrix Multiplication****'
    call count_cores()

    !   validate input
    call dimentionMatrix(n, key)
    !   allocate memory
    allocate(a(n,n),b(n,n),c(n,n),d(n,n),e(n,n))
    alpha=1.0;beta=1.0
    lda = n; ldb = n; lde = n

    IF(key.EQ.1)THEN
        
        call cpu_time(start)

        call ioMatrix(a, n, 'data1.dat')
        call ioMatrix(b, n, 'data2.dat')

        call cpu_time(finish)
        
        if(io==0)        write(unit=6,fmt=100) "Time to create matrix: ",finish-start," s."

    ELSE
        
        call cpu_time(start)

        call random_seed()
        do j=1, n
            do i=1, n
                call random_number(a(i,j))
                call random_number(b(i,j))
            enddo
        enddo
        
        call cpu_time(finish)
        if(io==0)        write(unit=6,fmt=100) "Time to create matrix: ",finish-start," s."

    ENDIF
    
    !---------------------------------------------------------------------
    ! a) triple do-loop.
    call cpu_time(start)
    c=0.0D0
    if (n < 1500) then
    !$OMP PARALLEL DO
    do j=1, n
        do i=1, n
            do k=1, n
                c(i,j)=c(i,j)+a(i,k)*b(k,j)
            enddo
        enddo
    enddo
    !$OMP END PARALLEL DO
    endif
    call cpu_time(finish)
    if(io==0)      write(unit=6,fmt=100) "Triple do-loop with threads: ",finish-start," s."
    var_triple = (finish-start)
    
    !---------------------------------------------------------------------
    ! b) matmul(a,b).
    call cpu_time(start)
    if (n < 1500) then
    d=0.0D0
    d=matmul(a,b)
    endif
    call cpu_time(finish)
    if(io==0)      write(unit=6,fmt=100) "Matmul(a,b): ",finish-start," s."
    var_matmul = (finish-start)

    !---------------------------------------------------------------------
    ! c) dgemm - INTEL MKL.
    call cpu_time(start)
    e=0.0D0
    call dgemm("N","N",n,n,n,alpha,a,lda,b,ldb,beta,e,lde)
    call cpu_time(finish)
    if(io==0)      write(unit=6,fmt=100) "DGEMM: ",finish-start," s."
    var_dgemm = (finish-start)

    !---------------------------------------------------------------------
    call MPI_Reduce(var_triple, var_triple_sum, 1, MPI_FLOAT, MPI_SUM, 0,        MPI_COMM_WORLD, ierror);
    call MPI_Reduce(var_matmul, var_matmul_sum, 1, MPI_FLOAT, MPI_SUM, 0,        MPI_COMM_WORLD, ierror);
    call MPI_Reduce(var_dgemm, var_dgemm_sum, 1,   MPI_FLOAT, MPI_SUM, 0,        MPI_COMM_WORLD, ierror);
    
    if(io==0) write(unit=6,fmt=100) "DGEMM AVG: ",var_dgemm_sum/nme," s."

    if(io==0)      open (1, file = 'data_times.dat', action='write', position='append') !check the file
    if(io==0)      write(1, *) nint(n * 1000.0) * 1E-3,var_triple,var_matmul,var_dgemm
    if(io==0)      close(1)

    deallocate(a,b,c,d,e)

    call MPI_FINALIZE(ierror)

    100 format(A,F8.3,A)

END PROGRAM
