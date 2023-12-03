program tm3d
    implicit none
    real(8), parameter :: PI=4.0*ATAN(1.0)
    integer, parameter :: prsc = 4
    integer, parameter :: range1 = 4
    integer, parameter :: range2 = 24
    integer, parameter :: numloops = 1  
    complex(8), dimension(4,4), parameter :: gamma1 = reshape([0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0],[4,4])
    complex(8), dimension(4,4), parameter :: gamma2 = reshape([cmplx(0,0),cmplx(0,0),cmplx(0,0),&
    cmplx(0,-1),cmplx(0,0), cmplx(0,0),cmplx(0,1),cmplx(0,0),cmplx(0,0),cmplx(0,-1),cmplx(0,0),&
    cmplx(0,0),cmplx(0,1),cmplx(0,0),cmplx(0,0),cmplx(0,0)],[4,4])
    complex(8), dimension(4,4), parameter :: gamma3 = reshape([0,0,1,0,0,0,0,-1,1,0,0,0,0,-1,0,0],[4,4])
    complex(8), dimension(4,4), parameter :: gamma4 = reshape([cmplx(0,0),cmplx(0,0),cmplx(0,-1),cmplx(0,0),&
    cmplx(0,0),cmplx(0,0),cmplx(0,0),cmplx(0,-1),cmplx(0,1),cmplx(0,0),cmplx(0,0),cmplx(0,0),cmplx(0,0),&
    cmplx(0,1),cmplx(0,0),cmplx(0,0)],[4,4])    
    complex(8), dimension(2,2) :: pauli1 = cmplx(reshape([[0,1,1,0]],[2,2]))
    complex(8), dimension(2,2):: pauli2 = reshape([cmplx(0,0),cmplx(0,-1),cmplx(0,1),cmplx(0,0)],[2,2])
    complex(8), dimension(2,2):: pauli3 = reshape([1,0,0,-1],[2,2])
    complex(8), dimension(2,2):: pauli4 = cmplx(0,-1)*reshape([1,0,0,1],[2,2])
    real(8), dimension(2*range1*prsc+1,range2*prsc+1,numloops) :: cherns
    real(8), dimension(2*range1*prsc+1,range2*prsc+1) :: cCherns
    real :: r
    integer :: k, loopnum, l, t1,t2

    call system_clock(t1,r)
    do loopnum = 1,numloops
        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) SHARED(cherns, loopnum)
        do k = 1,2*range1*prsc+1
            do l = 1,range2*prsc+1
                cherns(k,l,loopnum) = ly(4,75,real(k-range1*prsc-1)/real(prsc),real(l-1)/real(prsc),0.0)
            end do
        end do
        !$OMP END PARALLEL DO
    end do
    


!Writing results to be processed
    cCherns = averageOverIndex(cherns)
    open(1, file = 'data1.dat', status = 'old')
        do k = 1,2*range1*prsc+1
            do l = 1,range2*prsc+1
                write(1,*) cCherns(k,l)
            end do
        end do
    close(1)
    !Timer
    call system_clock(t2)
    print*, 'time taken: ',real((t2-t1))/real(r), 'seconds'

contains

    pure function delta(i,j) result(res)
        integer, intent(in) :: i,j
        integer :: res
        if(i == j) then
            res = 1
        else 
            res = 0
        end if 

        end function delta

    function transferm(ls,is,m,d1,d2) result(transfermatrix)
        integer, intent(in) :: ls, is
        real, intent(in) :: m, d1, d2
        complex(8), dimension(8*(ls**2),8*(ls**2)) :: transfermatrix0, transfermatrix
        complex(8), dimension(2*(ls**2),2*(ls**2)) :: h1,h2,hopping1,hopping2,hopping1i,hopping2i
        complex(4), dimension(2*(ls**2),2*(ls**2)) :: c1,c2,e1,e2
        complex(8), dimension(ls**2,ls**2) :: s, s1, t0, m0, t0i
        real(8), dimension(ls**2) :: randomM0, randomT0
        real(8), dimension(is,(ls**2)) :: dvt, dvm
        integer :: i,j

        randomT0 = 0.0
        randomM0 = 0.0

        s = shift(ls**2,1)
        s1 = shift(ls**2,-1)

        c1 = 0.5 * kron(pauli1,identityM(ls**2)) + &
        0.5 * kron(pauli2,kron(s+s1,identityM(ls))) + &
        0.5 * kron(pauli3,kron(identityM(ls),s+s1)) + &
        kron(pauli4, cmplx(0,-0.5) *  identityM(ls**2) + kron(s-s1,identityM(ls)) + &
        kron(identityM(ls),s-s1))
        
        c2 = 0.5 * kron(pauli1,identityM(ls**2))+ &
        0.5 * kron(pauli2,kron(transpose(s+s1),identityM(ls)))+ &
        0.5 * kron(pauli3,kron(identityM(ls),transpose(s+s1))) + &
        kron(pauli4, cmplx(0,-0.5) *  identityM(ls**2) + kron(transpose(s-s1),identityM(ls)) + &
        kron(identityM(ls),transpose(s-s1)))

        e1 = invM(c1)
        e2 = invM(c2)

        call random_number(dvt)
        call random_number(dvm)

        transfermatrix = identityM(8*(ls**2))

        do i = 1,is

            randomT0 = 1 + d2*(dvt(i,:)-0.5)
            randomM0 = m + d1*(dvm(i,:)-0.5)
            
            t0 = 0
            do j = 1,ls**2
                t0(j,j) = randomT0(j)
            end do 

            m0 = 0
            do j = 1,ls**2
                m0(j,j) = randomM0(j)
            end do 

            t0i = 0
            do j = 1,ls**2
                t0i(j,j) = randomT0(j)**(-1.0)
            end do 

            hopping1 = matmul(kron(identityM(2),t0),c1)
            hopping2 = matmul(kron(identityM(2),t0),c2)
            hopping1i = matmul(e1,kron(identityM(2),t0i))
            hopping2i = matmul(e2,kron(identityM(2),t0i))
            h1 = hopping1 + kron(pauli4,m0)

            transfermatrix0 = 0.0

            transfermatrix0(1:2*(ls**2),1:2*(ls**2)) = -matmul(conjg(transpose(hopping1i)),conjg(transpose(h1)))
            transfermatrix0(2*(ls**2)+1:4*(ls**2),2*(ls**2)+1:4*(ls**2)) = -matmul(conjg(transpose(hopping2i)),h1)

            transfermatrix0(1:2*(ls**2),1:2*(ls**2)) = -conjg(transpose(matmul(conjg(transpose(h1)),hopping1i)))
            transfermatrix0(2*(ls**2)+1:4*(ls**2),2*(ls**2)+1:4*(ls**2)) = -matmul(conjg(transpose(hopping2i)),h1)

            transfermatrix0(1:2*(ls**2),4*(ls**2)+1:6*(ls**2)) = -matmul(conjg(transpose(hopping1i)),hopping2)
            transfermatrix0(2*(ls**2)+1:4*(ls**2),6*(ls**2)+1:8*(ls**2)) = -matmul(conjg(transpose(hopping2i)),hopping1)

            transfermatrix0(4*(ls**2)+1:8*(ls**2),1:4*(ls**2)) = identityM(4*(ls**2))
            transfermatrix = matmul(transfermatrix,transfermatrix0)

        end do 


        end function transferm

    function identityM(n) result(eye)
        integer, intent(in) :: n
        integer :: i
        complex(8), dimension(n,n) :: eye

        eye = 0
        do i = 1,n
            eye(i,i) = 1
        end do

        end function identityM

    function hslice(n,dvt,dvm) result(H)
        integer, intent(in) :: n
        complex(8), dimension(n,n) :: S, invS
        complex(8), dimension(4 * (n ** 2),4 * (n ** 2)) :: H
        complex(8), dimension(n**2,n**2) :: randomT, randomM
        real(8), dimension(:) :: dvt, dvm
        integer :: i
        
        randomT = 0
        randomM = 0
        do i = 1,n**2
            randomT(i,i) = dvt(i)
            randomM(i,i) = dvm(i)
        end do 

        S = shift(n,1)
        invS = shift(n,-1)

        H = 0.5 * kron(gamma1,matmul(randomT,identityM(n**2))) + &
        0.5 * kron(gamma2,matmul(randomT,kron(S+invS,identityM(n)))) + &
        0.5 * kron(gamma3,matmul(randomT,kron(identityM(n),S+invS))) + &
        kron(gamma4, cmplx(0,-0.5) *  matmul(randomT,identityM(n**2) + kron(S-invS,identityM(n)) + &
        kron(identityM(n),S-invS))) + kron(gamma4, randomM)

        end function hslice

    function kron(A, B) result(C)
        complex(8), dimension(:,:), intent(in)  :: A, B
        complex(8), dimension(:,:), allocatable :: C
        integer :: i, j, m, n, p, q
    
        m = size(A, 1)
        n = size(A, 2)
        p = size(B, 1)
        q = size(B, 2)
    
        allocate(C(m*p, n*q))

        C = 0.0
    
        do i = 1, m
            do j = 1, n
                C((i-1)*p+1:i*p, (j-1)*q+1:j*q) = A(i, j) * B
            end do
        end do

        end function kron
        
    function shift(n,q) result(shiftM)
        integer, intent(in) :: n,q
        complex(8), dimension(n,n) :: shiftM
        integer :: i,j

        shiftM = 0

        do i = 1,n
            do j = 1,n
                if (i == mod(j+q,n)) then 
                    shiftM(i,j) = 1
                else
                    shiftM(i,j) = 0
                end if
            end do
        end do

        end function shift

    function averageOverIndex(T) result(collapsedT)
        real(8), intent(in), dimension(:,:,:) :: T
        real(8), dimension(size(T,1),size(T,2)) :: collapsedT
        integer :: i
    
        collapsedT = 0
        do i = 1, size(T,3)
            collapsedT = collapsedT + T(:,:,i)
        end do
        collapsedT = collapsedT / real(size(T,3))
        end function averageOverIndex

    function ly(ls,is,m,d1,d2) result(lyapunov)
        integer, intent(in) :: ls, is
        real, intent(in) :: m, d1, d2
        real(8) :: lyapunov
        real(8), dimension(8 * (ls ** 2)) :: eigvals
        complex(8), dimension(8 * (ls ** 2),8 * (ls ** 2)) :: tm
        complex(8), dimension(:), allocatable ::  work
        real(8), dimension(3 * 8 * (ls ** 2) - 2) :: rwork 
        integer :: info, lwork

        tm = transferm(ls,is,m,d1,d2)
        

        lyapunov = dlog(norm2(abs(tm)))/real(2*ls)



        end function ly

    function invT(diag,sign) result(Ti)
        real(8), dimension(:), intent(in) :: diag
        integer, intent(in) :: sign
        real(8), dimension(size(diag,1),size(diag,1)) :: Ti
        real(8), dimension(4*size(diag,1),4*size(diag,1)) :: hopping
        integer :: i, j, n
        
        n = size(diag,1)
        Ti = 0

        if(sign == 1) then
            do i = 1,n-1,2
                Ti(i,n) = ((-1)**(delta(mod(i,4),3)))
            end do 

            do i = 1,n-1,2
                Ti(n,i) = ((-1)**(delta(mod(i,4),3)))
            end do 

            do i = 1,n-1,2
                do j = 2,i-1,2
                    Ti(j,i) = ((-1)**(delta(mod(i+j,4),3)))
                    Ti(i,j) = ((-1)**(delta(mod(i+j,4),3)))
                end do 
            end do
        else
            do i = 1,n-1,2
                Ti(i,n) = mod(i,2)*(-1)
            end do 

            do i = 1,n-1,2
                Ti(n,i) = mod(i,2)*((-1)**(mod(i+j+1,2)))
            end do 

            do i = 1,n-1,2
                do j = 2,i-1,2
                    Ti(j,i) = ((-1)**(mod(i+j+1,2)))
                    Ti(i,j) = ((-1)**(mod(i+j,2)))
                end do 
            end do
        end if
        end function invT
        
    function invM(A) result(Ainv)
        complex(4), dimension(:,:), intent(in) :: A
        complex(4), dimension(size(A,1),size(A,2)) :: Ainv
        
        complex(4), dimension(size(A,1)) :: work  ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: n, info
        
        ! External procedures defined in LAPACK
        external DGETRF
        external DGETRI
        
        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)
        
        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call CGETRF(n, n, Ainv, n, ipiv, info)
        
        if (info /= 0) then
            stop 'Matrix is numerically singular!'
        end if
        
        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call CGETRI(n, Ainv, n, ipiv, work, n, info)
        
        if (info /= 0) then
            stop 'Matrix inversion failed!'
        end if
        end function invM

    

end program tm3d
