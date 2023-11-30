program onedim

    use random

    implicit none

    real, parameter :: PI=4.0*ATAN(1.0)
    integer, parameter :: prsc = 64
    integer, parameter :: range1 = 2
    integer, parameter :: range2 = 6
    integer, parameter :: numloops = 8
    integer, parameter :: latticesites = 200
    complex(4), dimension(2,2) :: pauli1 = cmplx(reshape([[0,1,1,0]],[2,2]))
    complex(4), dimension(2,2):: pauli2 = reshape([cmplx(0,0),cmplx(0,-1),cmplx(0,1),cmplx(0,0)],[2,2])
    real(4), dimension(2*range1*prsc+1,range2*prsc+1,numloops) :: cherns
    real(4), dimension(2*range1*prsc+1,range2*prsc+1) :: cCherns
    real :: r
    integer :: k, loopnum, l, t1,t2

    call system_clock(t1,r)

    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(NONE) SHARED(cherns)
    do loopnum = 1,numloops
        do k = 1,2*range1*prsc+1
            do l = 1,range2*prsc+1
                cherns(k,l,loopnum) = real(transferMatrix(latticesites,real(k-range1*prsc-1)/real(prsc),&
                real(l-1)/real(prsc),real(l-1)/real(2*prsc)))
            end do
        end do
    end do
    !$OMP END PARALLEL DO


!Writing results to be processed
    cCherns = averageOverIndex(cherns)
    open(1, file = 'data1.dat', status = 'old')
        do k = 1,2*range1*prsc+1
            do l = 1,range2*prsc+1
                write(1,*) cCherns(k,l)
            end do
        end do
    close(1)

    call system_clock(t2)
    print*, 'time taken: ',real((t2-t1))/real(r), 'seconds'

contains

function transferMatrix(n,m,d1,d2) result(ly)
    integer, intent(in) :: n
    real, intent(in) :: d1,d2,m
    complex(4), dimension(2,2) :: transfer
    complex(4) ::  entry1, entry2
    real(4), dimension(n) :: disordervector
    real :: tn, mn, tni, mni, ly
    integer :: i


    call random_number(disordervector)

    do i = 1,n
        disordervector(i) = disordervector(i)-0.5
    end do 

    tn = 1 + d2*(disordervector(1))
    mn = m + d1*(disordervector(1))


    transfer = reshape([cmplx(0,-1)*(mn**2)/(tn * mn), cmplx(-1.0,0.0), cmplx(1.0,0.0), cmplx(0.0,0.0)],[2,2])

    do i = 2,n
        tn = 1 + d2*(disordervector(i))
        mn = m + d1*(disordervector(i))
        tni = 1 + d2*(disordervector(i-1))
        mni = m + d1*(disordervector(i-1))
        entry1 = cmplx(0,-1)*(mn**2 - tni**2)/(tn * mn)
        entry2 = -(tni * mni)/(tn * mn)
        transfer = matmul(transfer, reshape([entry1, entry2, cmplx(1.0,0.0), cmplx(0.0,0.0)],[2,2]))
    end do 
    
    ly = log(norm2(abs(transfer)))/real(n)

    if(ly /= ly) then
        ly = 0
    else if(abs(ly)<0.01) then 
        ly = 0
    else 
        ly = 1/ly
    end if

    end function transferMatrix

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


function hamiltonian(n,m,d1,d2,omega) result(H)
    integer, intent(in) :: n
    complex(4), dimension(n,n) :: S
    real, intent(in) :: m, d1, d2
    real(4), dimension(n), intent(in) :: omega
    complex(4), dimension(2 * n,2 * n) :: H
    complex(4), dimension(n,n) :: randomT, randomM
    integer :: i 

    randomT = 0.0
    randomM = 0.0
    do i = 1,n
        randomM(i,i) = m + d1 * (omega(i)-0.5)
        randomT(i,i) = 1 + d2 * (omega(i)-0.5)
    end do 


    S = shift(n)
    H = 0.5 * cmplx(0,-1) * kron(pauli1, matmul(randomT,S - transpose(S))) + kron(pauli2, &
    randomM + 0.5 * matmul(randomT,(S + transpose(S))))

    end function hamiltonian

    function evs(M,n) result(eigvals)
        integer, intent(in) :: n 
        complex(4), dimension(2*n,2*n), intent(in) :: M
        real, dimension(2*n) :: eigvals
        complex(4), dimension(2*n,2*n) :: temp
        complex(4), dimension(:), allocatable ::  work
        real, dimension(6*n - 2) :: rwork
        integer :: info, lwork
        
        temp = M

        lwork = -1
        allocate(work(1))
        call cheev('N', 'U', 2*n, temp, 2*n, eigvals, work, lwork, rwork, info)
        lwork = ceiling(real(work(1)))
        deallocate(work)
    
    
        allocate(work(lwork))
        work = 0
        call cheev('N', 'U', 2*n, temp, 2*n, eigvals, work, lwork, rwork, info)
        deallocate(work)
        
    end function evs

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
