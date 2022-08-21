program sdriver

! Sample Size


! Timing Variables
real (kind=8) :: cpuStart, cpuEnd, cputime
real (kind=8) :: wallStart, wallEnd, walltime
external cputime
external walltime
!& Other Vars

#ifdef DOT
real (kind=8), DIMENSION(:), allocatable :: vA, vB
real (kind=8) :: dot
#elif VVM
real (kind=8), DIMENSION(:), allocatable :: vA, vB
real (kind=8), DIMENSION(:,:), allocatable :: C
#elif MVV
real (kind=8), DIMENSION(:), allocatable :: vA, vB
real (kind=8), DIMENSION(:,:), allocatable :: A
#elif MMM
real (kind=8), DIMENSION(:,:), allocatable :: A, B, C
#elif DLS
real (kind=8), DIMENSION(:,:), allocatable :: A
real (kind=8), DIMENSION(:), allocatable :: vA, vB
#endif

real (kind=8) :: flops, flops2, ans
character (len=8) :: carg1, carg2
integer numThreads, N

! Number Of Threads We Want to Use
! ****** Select 1 for serial ******

call get_command_argument(1,carg1)
call get_command_argument(2,carg2)
read(carg1, '(i8)') N
read(carg2, '(i8)') numThreads


#ifdef DOT
       
        ! Allocate memory to vectors
        allocate( vA(N), stat=ierr)
        allocate( vB(N), stat=ierr)   
        
        ! Fill Vectors
        do i=1, N
                vA(i) = dble(i)
                vB(i) = 1.0/dble(i)
        enddo
        
        ! Get Dot Product of Vectors and Time
        wallStart = walltime()
        cpuStart = cputime()

        call omp_set_num_threads(1)
!        ans = ddot(N,vA,0,vB,0)
        ans = dot(numThreads, N, vA, vB)

        cpuEnd = cputime()
        wallEnd = walltime()

        flops = (2.0D0*N/(cpuEnd-cpuStart)/1e6)
        flops2 = (2.0D0*N/(wallEnd-wallStart)/1e6)        
        
        ! Clean Up
        if (allocated(vA)) deallocate(vA)
        if (allocated(vB)) deallocate(vB)

        ! Show Results
        print *, " Dot Product : ", ans 
        print *, " Time        : ", (cpuEnd-cpuStart)
        print *. " Wall Time   : ", (wallEnd-wallStart)
        print *, " Mega Flops  : ", flops
        print *, " Wall Flops  : ", flops2
        
#elif VVM
        
        ! Allocate memory
        allocate( vA(N), stat=ierr)
        allocate( vB(N), stat=ierr)
        allocate( C(N,N), stat=ierr)

        ! Fill Vectors
        do i=1, N
                vA(i) = dble(i)
                vB(i) = 1/dble(i)
        enddo
        
        ! Run and Time VVM
        wallStart = walltime()
        cpuStart = cputime()
        call vvm( numThreads, N, vA, vB, C) 
        cpuEnd = cputime()
        wallEnd = walltime()

        flops = N*N/(cpuEnd-cpuStart)/1e6
        flops2 = N*N/(wallEnd-wallStart)/1e6
        ! Prints Resultant Matrix
!        print *, "Resultant Matrix: "
!        do i=1, N
!                do j=1, N
!                        write(*,"(ES8.1,$)") C(i,j)
!                enddo
!                write(*,*) ''
!        enddo

        print *, " *** VVM ***  "
        print *, " N Size     : ", N
        print *, " Time       : ", (cpuEnd-cpuStart)
        print *, " Wall Time  : ", (wallEnd-wallStart)
        print *, " Mega Flops : ", flops
        print *, " Wall Flops : ", flops2

        ! Clean Up
        if (allocated(vA)) deallocate(vA)
        if (allocated(vB)) deallocate(vB)
        if (allocated(C)) deallocate(C)

#elif MVV

        ! Allocate Memory
        allocate( A(N,N), stat=ierr)
        allocate( vA(N), stat=ierr)
        allocate( vB(N), stat=ierr)

        ! Fill matrix and vector
        do i=1, N
                vA(i) = i
                do j=1, N
                        A(i,j) = 1
                enddo
        enddo

        ! Initialize vB        
        vB = 0.0
     
        ! Run and Time MVV
        wallStart = walltime()
        cpuStart = cputime()
        call mvv(numThreads, N, A, vA, vB)
        cpuEnd = cputime()
        wallEnd = walltime()
        
        ! Calculate Megaflops
        flops = 3.0*N*N/(cpuEnd-cpuStart)/1e6
        flops2 = 3.0*N*N/(wallEnd-wallStart)/1e6


        ! Prints Resultant Vector
!        print *, "Resultant Vector"
!        do i=1, N
!                write(*,"(ES8.1,$)") vB(i)
!        enddo

        print *, " *** MVV ***   "
        print *, " N           : ", N
        print *, " Time        : ", (cpuEnd-cpuStart)
        print *, " Wall Time   : ", (wallEnd-wallStart)
        print *, " Cpu Flops   : ", flops
        print *, " Wall mFlops : ", flops2


#elif MMM

        ! Allocate Memory
        allocate( A(N,N), stat=ierr)
        allocate( B(N,N), stat=ierr)
        allocate( C(N,N), stat=ierr)

        do i=1, N
                do j=1, N
                        A(i,j) = i
                        B(i,j) = j
                enddo
        enddo

        cpuStart = cputime()
        call mmm(numThreads, N, A, B, C)
        cpuEnd = cputime()

        ! Calculate Megaflops
        flops = 2.0*dble(N)**3.0/ (cpuEnd-cpuStart) / 1.0e6

        ! Print Results
        print *, " *** MMM ***  "
        print *, " N          : ", N
        print *, " Time       : ", (cpuEnd-cpuStart)
        print *, " Mega Flops : ", flops
        

#elif DLS
        
        allocate( A(N,N), stat=ierr)
        allocate( vA(N), stat=ierr)
        allocate( vB(N), stat=ierr)

        do i=1, N
                vA(i) = 0
                do j=1, N
                        A(i,j) = j+1
                enddo
        enddo

        cpuStart = cputime()
        call dls(numThreads, N, A, vA, vB)
        cpuEnd = cputime()

        ! Calculate Megaflops
        flops = (2.0/3.0)*N*N*N/(cpuEnd-cpuStart)/100000.0

        ! Print Results
        print *, " *** DLS ***  "
        print *, " N          : ", N
        print *, " Time       : ", (cpuEnd-cpuStart)
        print *, " Mega Flops : ", flops

#endif

end program sdriver
















