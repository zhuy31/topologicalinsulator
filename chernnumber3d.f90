program threedim
    implicit none
    real, parameter :: PI=4.0*ATAN(1.0)
    integer, parameter :: prsc = 4
    integer, parameter :: range1 = 4
    integer, parameter :: range2 = 6
    integer, parameter :: numloops = 4
    integer, parameter :: latticesites = 8
    complex(4), dimension(4,4), parameter :: gamma1 = reshape([0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0],[4,4])
    complex(4), dimension(4,4), parameter :: gamma2 = reshape([cmplx(0,0),cmplx(0,0),cmplx(0,0),&
    cmplx(0,-1),cmplx(0,0), cmplx(0,0),cmplx(0,1),cmplx(0,0),cmplx(0,0),cmplx(0,-1),cmplx(0,0),&
    cmplx(0,0),cmplx(0,1),cmplx(0,0),cmplx(0,0),cmplx(0,0)],[4,4])
    complex(4), dimension(4,4), parameter :: gamma3 = reshape([0,0,1,0,0,0,0,-1,1,0,0,0,0,-1,0,0],[4,4])
    complex(4), dimension(4,4), parameter :: gamma4 = reshape([cmplx(0,0),cmplx(0,0),cmplx(0,-1),cmplx(0,0),&
    cmplx(0,0),cmplx(0,0),cmplx(0,0),cmplx(0,-1),cmplx(0,1),cmplx(0,0),cmplx(0,0),cmplx(0,0),cmplx(0,0),&
    cmplx(0,1),cmplx(0,0),cmplx(0,0)],[4,4])
    real(4), dimension(2*range1*prsc+1,range2*prsc+1,numloops) :: cherns
    real(4), dimension(2*range1*prsc+1,range2*prsc+1) :: cCherns
    real :: r
    integer :: k, loopnum, l, t1,t2

    call system_clock(t1,r)

   
    do loopnum = 1,numloops
        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) SHARED(cherns, loopnum)
        do k = 1,2*range1*prsc+1
            do l = 1,range2*prsc+1
                cherns(k,l,loopnum) = chernNumber(latticesites,real(k-range1*prsc-1)/real(prsc),real(l-1)/real(prsc),real(l-1)/real(2*prsc))
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

    elemental function sgn(x) result(y)
        real, intent(in) :: x
        real :: y
        if (x > 0) then
            y = 1
        else if (x == 0) then
            y = 0
        else
            y = -1
        end if

    end function sgn



    function shift(n,q) result(shiftM)
    integer, intent(in) :: n,q
    complex(4), dimension(n,n) :: shiftM
    integer :: i,j

        shiftM = 0

        do i = 1,n
            do j = 1,n
                if (mod(i,n) == mod(j+q,n)) then 
                    shiftM(i,j) = 1
                else
                    shiftM(i,j) = 0
                end if
            end do
        end do

        end function shift

    function kron(A, B) result(C)
        complex(4), dimension(:,:), intent(in)  :: A, B
        complex(4), dimension(:,:), allocatable :: C
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
    
    function identityM(n) result(eye)
        integer, intent(in) :: n
        integer :: i
        complex(4), dimension(n,n) :: eye

        eye = 0
        do i = 1,n
            eye(i,i) = 1
        end do

        end function identityM

    function outer(u,v) result(resultant)
        complex(4), dimension(:), intent(in) :: u,v
        complex(4), dimension(size(u),size(u)) :: resultant
        integer :: i

        do i = 1, size(v)
            resultant(:,i) = conjg(v(i)) * u
        end do

        end function outer

    function randomDiag(m,d,n) result(rDiag)
        real, intent(in) :: d, m
        integer, intent(in) :: n
        complex(4), dimension(n,n) :: rDiag
        real(4), dimension(n) :: diag, diag1
        integer :: i
        rDiag = 0
        call random_number(diag)
        call random_number(diag1)
        do i = 1,n
            diag(i) = (-2*log(diag(i))) ** (0.5) * cos(2 * PI * diag1(i))
        end do 
        do i = 1,n
            rDiag(i,i) = m + d * (diag(i) + 0.5 * diag(I) ** 2 + (1.0/6.0) * diag(I) ** 3)
        end do 

        end function randomDiag

    function hamiltonian(n,m,d1,d2) result(H)
        integer, intent(in) :: n
        complex(4), dimension(n,n) :: S, invS
        real, intent(in) :: m, d1, d2
        complex(4), dimension(4 * (n ** 3),4 * (n ** 3)) :: H
        complex(4), dimension(n**3,n**3) :: randomT, randomM


        randomT = randomDiag(1.0,d2,n**3)
        randomM = randomDiag(m,d1,n**3)
        S = shift(n,1)
        invS = shift(n,-1)

        H = 0.5 * kron(gamma1,matmul(randomT,kron(S+invS,identityM(n**2)))) + &
        0.5 * kron(gamma2,matmul(randomT,kron(kron(identityM(n),S+invS),identityM(n)))) + &
        0.5 * kron(gamma3,matmul(randomT,kron(identityM(n**2),S+invS))) + &
        kron(gamma4, cmplx(0,-0.5) *  matmul(randomT,kron(S-invS,identityM(n**2)) + kron(kron(identityM(n),S-invS),identityM(n)) + &
        kron(identityM(n**2),S-invS))) + kron(gamma4, randomM)

    end function hamiltonian

    function mSgn(M,n) result(sgnM)
        integer, intent(in) :: n 
        complex(4), dimension(4 * (n ** 3),4 * (n ** 3)), intent(in) :: M
        real, dimension(4 * (n ** 3)) :: eigvals
        complex(4), dimension(4 * (n ** 3),4 * (n ** 3)) :: temp, sgnM
        real(4), dimension(4 * (n ** 3),4 * (n ** 3)):: es
        complex(4), dimension(:), allocatable ::  work
        real, dimension(3 * 4 * (n ** 3) - 2) :: rwork 
        integer :: info, i, lwork

 
        sgnM = 0.0
        temp = M
        rwork = 0.0
        lwork = -1
        allocate(work(1))
        call cheev('V', 'U', 4 * (n ** 3), temp, 4 * (n ** 3), eigvals, work, lwork, rwork, info)
        lwork = ceiling(real(work(1)))
        deallocate(work)
    
        allocate(work(lwork))
        call cheev('V', 'U', 4 * (n ** 3), temp, 4 * (n ** 3), eigvals, work, lwork, rwork, info)
        deallocate(work)

        es = 0
        do i = 1, 4 * (n ** 3)
            es(i,i) = eigvals(i)
        end do 

        es = sgn(es)

        sgnM = matmul(temp, matmul(es, conjg(transpose(temp))))

    end function mSgn


    function partial(M,index) result(pU)
        integer, intent(in) :: index
        complex(4), dimension(:,:), intent(in) :: M
        complex(4), dimension(size(M,1),size(M,1)) :: pU, derivative
        integer :: i, j, n, x, y
        pU = 0
        derivative = 0
        n = int((real(size(M,1))/2.0 )**(1.0/3.0))
    
        if (index == 1) then 
            do i = 1,size(M,1)
                do j = 1,size(M,1)
                    x = mod(mod(mod(j, n ** 3), n ** 2), n)
                    y = mod(mod(mod(i, n ** 3), n ** 2), n)
                    derivative(i,j) =  (mod(abs(x-y+(n/2)),n) - (n/2))
                end do
            end do
        else if (index == 2) then
            do i = 1,size(M,1)
                do j = 1,size(M,1)
                    x = mod(mod(j, n ** 3), n ** 2) / n
                    y = mod(mod(i, n ** 3), n ** 2) / n 
                    derivative(i,j) = (mod(abs(x-y+(n/2)),n) - (n/2))
                end do
            end do
        else 
            do i = 1,size(M,1)
                do j = 1,size(M,1)
                    x = mod(j, n ** 3) / (n ** 2)
                    y = mod(i, n ** 3) / (n ** 2)
                    derivative(i,j) =  (mod(abs(x-y+(n/2)),n) - (n/2))
                end do
            end do    
        end if
    
    
        pU = M * derivative

        end function partial
    
    function tr(M) result(trace)
        complex(4), dimension(:,:), intent(in) :: M
        integer :: i,n
        complex(4) :: trace
    
        n = size(M,1)
        trace = 0
    
        do i=1,n
            trace = trace + M(i,i)
        end do
    
        end function tr
    
    function chernNumber(n,m,d1,d2) result(ch)
        integer, intent(in) :: n
        complex(4), dimension(4 * (n ** 3), 4 * (n ** 3)) :: Ham
        complex(4), dimension(2 * (n ** 3), 2 * (n ** 3)) :: U, d1U, d2U, d3U, integrand
        real(4), intent(in) :: m, d1, d2
        real(4) :: ch 
        
        Ham = mSgn(hamiltonian(n,m,d1,d2),n)
        U = Ham(1:2 * (n ** 3),2 * (n ** 3)+1:4 * (n ** 3))
    
        d1U = partial(U,1)
        d2U = partial(U,2)
        d3U = partial(U,3)
    
        integrand = matmul(conjg(transpose(U)),matmul(d1U,matmul(d2U,d3U)) + matmul(d3U,matmul(d1U,d2U)) + &
        matmul(d2U,matmul(d3U,d1U)) - matmul(d1U,matmul(d3U,d2U)) - &
        matmul(d2U,matmul(d1U,d3U)) - matmul(d3U,matmul(d2U,d1U)))
        
        ch = real(real(PI/3.0) * (tr(integrand)/( n ** 3))*cmplx(0,1))

        end function chernNumber

    function averageOverIndex(T) result(collapsedT)
        real, intent(in), dimension(:,:,:) :: T
        real, dimension(size(T,1),size(T,2)) :: collapsedT
        integer :: i
    
        collapsedT = 0
        do i = 1, size(T,3)
            collapsedT = collapsedT + T(:,:,i)
        end do
        collapsedT = collapsedT / real(size(T,3))
    end function averageOverIndex


        
end program threedim
