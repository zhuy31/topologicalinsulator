program onedim
    implicit none
    real, parameter :: PI=4.0*ATAN(1.0)
    integer, parameter :: prsc = 4
    integer, parameter :: range1 = 4
    integer, parameter :: range2 = 32
    integer, parameter :: numloops = 4
    integer, parameter :: latticesites = 100
    complex(4), dimension(2,2) :: pauli1 = cmplx(reshape([[0,1,1,0]],[2,2]))
    complex(4), dimension(2,2):: pauli2 = reshape([cmplx(0,0),cmplx(0,-1),cmplx(0,1),cmplx(0,0)],[2,2])
    real(4), dimension(2*range1*prsc+1,range2*prsc+1,numloops) :: cherns
    real(4), dimension(2*range1*prsc+1,range2*prsc+1) :: cCherns
    real :: r
    integer :: k, loopnum, l, t1,t2

    call system_clock(t1,r)
    

    !$OMP PARALLEL DO COLLAPSE(3) NUM_THREADS(12) DEFAULT(NONE) SHARED(cherns)
        do loopnum = 1,numloops
            do k = 1,2*range1*prsc+1
                do l = 1,range2*prsc+1
                    cherns(k,l,loopnum) = chernNumber(latticesites,&
                    real(k-1)/real(prsc)-real(range1),real(l-1)/real(prsc),0.5*real(l-1)/real(prsc))
                end do
            end do
        end do
    !$OMP END PARALLEL DO
    

    !Writing results to be processed
    cCherns = averageOverIndex(cherns)
    open(1, file = 'cherndata1d.dat', status = 'old')
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



function randomHermitian(n) result(rh)
    integer, intent(in) :: n
    real, dimension(n,n) :: m1,m2
    complex(4), dimension(n,n) :: rh
    complex(4), dimension(n,n) :: cx

    call random_number(m1)
    call random_number(m2)
    cx = cmplx(m1,m2)
    rh = 0.5 * (cx + transpose(conjg(cx)))

    end function randomHermitian

function normalize(v) result(unitV)
    complex(4), dimension(:), intent(in) :: v
    complex(4), dimension(size(v)) :: unitV
    complex(4) :: norm
    integer :: i
    norm = 0
    do i = 1,size(v)
        norm = norm + abs(v(i))**2
    end do
    unitV = v * (1/norm)

    end function normalize

function shift(n) result(shiftM)
    integer, intent(in) :: n
    complex(4), dimension(n,n) :: shiftM
    integer :: i,j

    do i = 1,n
        do j = 1,n
            if (mod(i,n) == mod(j+1,n)) then 
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

function randomDiag(mass,d,n) result(rDiag)
    real, intent(in) :: d, mass
    integer, intent(in) :: n
    complex(4), dimension(n,n) :: rDiag
    real(4), dimension(n) :: diag, diag1
    integer :: i
    rDiag = 0
    call random_number(diag1)
    call random_number(diag)
    do i = 1,n
        rDiag(i,i) = (mass) + d*10*log(diag(i)+0.01)
    end do 

    end function randomDiag

function hamiltonian(n,mass,d1,d2) result(H)
    integer, intent(in) :: n
    complex(4), dimension(n,n) :: S
    real, intent(in) :: mass, d1, d2
    complex(4), dimension(2 * n,2 * n) :: H
    complex(4), dimension(n,n) :: randomT, randomM

    randomT = randomDiag(1.0,d2,n)
    randomM = randomDiag(mass,d1,n)

    S = shift(n)
    H = 0.5 * cmplx(0,-1) * kron(pauli1, matmul(randomT,S - transpose(S))) + kron(pauli2, &
    randomM + 0.5 * matmul(randomT,(S + transpose(S))))

    end function hamiltonian

    function spectralInv(M,n) result(sgnM)
        integer, intent(in) :: n 
        complex(4), dimension(2*n,2*n), intent(in) :: M
        real, dimension(2*n) :: eigvals
        complex(4), dimension(2*n,2*n) :: temp, sgnM
        real(4), dimension(2*n,2*n) :: es
        complex(4), dimension(:), allocatable ::  work
        real, dimension(6*n - 2) :: rwork
        integer :: info,  i, lwork
        

        sgnM = 0
        temp = M
    
        lwork = -1
        allocate(work(1))
        call cheev('V', 'U', 2*n, temp, 2*n, eigvals, work, lwork, rwork, info)
        lwork = ceiling(real(work(1)))
        deallocate(work)
    
    
        allocate(work(lwork))
        work = 0
        call cheev('V', 'U', 2*n, temp, 2*n, eigvals, work, lwork, rwork, info)
        deallocate(work)
        
        es = 0

        do i = 1, 2*n
            es(i,i) = eigvals(i)
        end do 

        es = sgn(es)
        sgnM = matmul(temp, matmul(es, conjg(transpose(temp))))

    end function spectralInv


    function matrixDerivative(U,n) result(pU)
        integer, intent(in) :: n
        complex(4), dimension(n,n), intent(in) :: U
        complex(4), dimension(n,n) :: pU
        integer :: i, j
    
        do i = 1,n
            do j = 1,n
                pU(i,j) = cmplx(U(i,j) * (mod(abs(j-i+(n/2)),n) - (n/2)))
            end do
        end do
        
        end function matrixDerivative
    
    function tr(M,n) result(trace)
        integer, intent(in) :: n
        complex(4), dimension(n,n), intent(in) :: M
        integer :: i
        complex(4) :: trace
    
        trace = 0
        do i=1,n
            trace = trace + M(i,i)
        end do
    
        end function tr
    
    function chernNumber(n,mass,d1,d2) result(ch)
        integer, intent(in) :: n
        complex(4), dimension(2*n,2*n) :: H
        real, intent(in) :: mass, d1, d2
        complex(4), dimension(n,n) :: U, dU
        real(4) :: ch
        
        H = spectralInv(hamiltonian(n,mass,d1,d2),n)
        U = H(1:n,n+1:2*n)
        dU = matrixDerivative(U,n)
        ch = real(tr(matmul(conjg(transpose(U)),dU),n))/(real(n))
        
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
        
end program onedim
