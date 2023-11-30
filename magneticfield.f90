program magneticfield
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
    real(4) :: r
    integer :: k, loopnum, l, t1,t2, iter

    
    call system_clock(t1,r)

    do iter = 1,4
        do loopnum = 1,numloops
            !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) SHARED(cherns, loopnum,iter)
            do k = 1,2*range1*prsc+1
                do l = 1,range2*prsc+1
                    cherns(k,l,loopnum) = real(chernNumber(latticesites,real(k-range1*prsc-1)/real(prsc),&
                    real(l-1)/real(prsc),real(l-1)/real(2*prsc),0.75,iter))
                end do
            end do
            !$OMP END PARALLEL DO
        end do

        !Writing results to be processed
        cCherns = averageOverIndex(cherns)
        open(1, file = 'data1.dat', status = 'old',position='append')
            do k = 1,2*range1*prsc+1
                do l = 1,range2*prsc+1
                    write(1,*) cCherns(k,l)
                end do
            end do
        close(1)
    end do 

    !Timer
    call system_clock(t2)
    print*, 'time taken: ',real((t2-t1))/real(r), 'seconds'

contains
    pure function skd(x,y,q,n) result(z)
        integer, intent(in) :: x,y,q,n
        complex :: z
        
        if (mod(x,n) == mod(mod(y+q,n)+n,n)) then
            z = cmplx(1,0)
        else 
            z = cmplx(0,0)
        end if 

        end function skd

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

    function randomDiag(m,d,n,index) result(rDiag)
        real, intent(in) :: d, m
        integer, intent(in) :: n, index
        complex(4), dimension(n,n) :: rDiag
        real(4), dimension(n) :: diag, diag1
        integer :: i
        rDiag = 0

        call random_number(diag)
        call random_number(diag1)

        do i = 1,n
            diag(i) = (-2*log(diag(i))) ** (0.5) * cos(2 * PI * diag1(i))
        end do 

        if(index == 1) then
            do i = 1,n
                rDiag(i,i) = m + d * (diag(i))
            end do 
        else if(index == 2) then
            do i = 1,n
                rDiag(i,i) = m + d * (diag(i) + 0.5*(diag(i)**2))
            end do 
        else if(index == 3) then
            do i = 1,n
                rDiag(i,i) = m + d * (diag(i) + 0.5*(diag(i)**2) + (1.0/6.0) * (diag(i)**3))
            end do 
        else if(index == 4) then
            do i = 1,n
                rDiag(i,i) = m + d * exp(diag(i))
            end do 
        else
            rDiag = 0.0
        end if 
    

        end function randomDiag

    function hamiltonian(n,m,d1,d2,B,index) result(H)
        integer, intent(in) :: n,index
        real, intent(in) :: B
        complex(4), dimension(n,n) :: S, invS
        real, intent(in) :: m, d1, d2
        complex(4), dimension(4 * (n ** 3),4 * (n ** 3)) :: H
        complex(4), dimension(n**3,n**3) :: randomT, randomM


        randomT = randomDiag(1.0,d2,n**3,index)
        randomM = randomDiag(m,d1,n**3,index)
        S = shift(n,1)
        invS = shift(n,-1)

        H = kron(gamma1,matmul(randomT,mi1(n,1,1,B))) + &
        kron(gamma2,matmul(randomT,mi1(n,2,1,B))) + &
        kron(gamma3,matmul(randomT,mi1(n,3,1,B))) + &
        kron(gamma4,matmul(randomT,mi1(n,1,-1,B)+mi1(n,2,-1,B)+mi1(n,3,-1,B))) + kron(gamma4, randomM)
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
                    derivative(i,j) =  (mod(abs(x-y+(n/2)),n) - (n/2)) * cmplx(1,0)
                end do
            end do
        else if (index == 2) then
            do i = 1,size(M,1)
                do j = 1,size(M,1)
                    x = mod(mod(j, n ** 3), n ** 2) / n
                    y = mod(mod(i, n ** 3), n ** 2) / n 
                    derivative(i,j) = (mod(abs(x-y+(n/2)),n) - (n/2)) * cmplx(1,0)
                end do
            end do
        else 
            do i = 1,size(M,1)
                do j = 1,size(M,1)
                    x = mod(j, n ** 3) / (n ** 2)
                    y = mod(i, n ** 3) / (n ** 2)
                    derivative(i,j) =  (mod(abs(x-y+(n/2)),n) - (n/2)) * cmplx(1,0)
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
    
    function chernNumber(n,m,d1,d2,B,index) result(ch)
        integer, intent(in) :: n, index
        real, intent(in) :: B
        complex(4), dimension(4 * (n ** 3), 4 * (n ** 3)) :: Ham
        complex(4), dimension(2 * (n ** 3), 2 * (n ** 3)) :: U, d1U, d2U, d3U, integrand
        real(4), intent(in) :: m, d1, d2
        complex(4) :: ch 
        
        Ham = mSgn(hamiltonian(n,m,d1,d2,B,index),n)
        U = Ham(1:2 * (n ** 3),2 * (n ** 3)+1:4 * (n ** 3))
    
        d1U = partial(U,1)
        d2U = partial(U,2)
        d3U = partial(U,3)
    
        integrand = matmul(conjg(transpose(U)),matmul(d1U,matmul(d2U,d3U)) + matmul(d3U,matmul(d1U,d2U)) + &
        matmul(d2U,matmul(d3U,d1U)) - matmul(d1U,matmul(d3U,d2U)) - &
        matmul(d2U,matmul(d1U,d3U)) - matmul(d3U,matmul(d2U,d1U)))
        
        ch = real(PI/3.0) * (tr(integrand)/( n ** 3)) * cmplx(0,1)

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
        
    function mI1(n,index,sign,B) result(matrix)
        integer, intent(in) :: n, index, sign
        complex(4), dimension(n ** 3, n ** 3) :: matrix
        real(4), intent(in) :: B
        integer :: i, j
        matrix = 0

        if(index == 1) then
            if(sign == 1) then
                do i = 0,(n**3)-1
                    do j = 0,(n**3)-1
                        matrix(i+1,j+1) = (skd(i/(n**2),j/(n**2),1,n)*exp(cmplx(0,1)*B*(PI/n)*((mod(i,n**2)/n)+(mod(i,n))))+&
                        skd(i/(n**2),j/(n**2),-1,n)*exp(cmplx(0,-1)*B*(PI/n)*((n/2 - mod(i,n**2)/n)+(n/2 - mod(i,n)))))*&
                        skd((mod(i,n**2)/n),(mod(j,n**2)/n),0,n)*skd(mod(i,n),mod(j,n),0,n)*0.5
                    end do 
                end do 
            else
                do i = 0,(n**3)-1
                    do j = 0,(n**3)-1
                        matrix(i+1,j+1) = (skd(i/(n**2),j/(n**2),1,n)*exp(cmplx(0,1)*B*(PI/n)*((mod(i,n**2)/n)+(mod(i,n))))-&
                        skd(i/(n**2),j/(n**2),-1,n)*exp(cmplx(0,-1)*B*(PI/n)*((n/2 - mod(i,n**2)/n)+(n/2 - mod(i,n)))))*&
                        skd((mod(i,n**2)/n),(mod(j,n**2)/n),0,n)*skd(mod(i,n),mod(j,n),0,n)*cmplx(0,-0.5)
                    end do 
                end do 
            end if
        else if(index == 2) then
            if(sign == 1) then
                do i = 0,(n**3)-1
                    do j = 0,(n**3)-1
                        matrix(i+1,j+1) = (skd(mod(i,n**2)/n,mod(j,n**2)/n,1,n)*exp(cmplx(0,1)*B*(PI/n)*(mod(i,n) - (i / (n**2))))+&
                        skd(mod(i,n**2)/n,mod(j,n**2)/n,-1,n)*exp(cmplx(0,-1)*B*(PI/n)*(n/2 - mod(i,n) - (n/2 - i / (n**2)))))*&
                        skd(mod(i,n),mod(j,n),0,n)*skd(i/(n**2),j/(n**2),0,n)*0.5
                    end do 
                end do 
            else
                do i = 0,(n**3)-1
                    do j = 0,(n**3)-1
                        matrix(i+1,j+1) = (skd(mod(i,n**2)/n,mod(j,n**2)/n,1,n)*exp(cmplx(0,1)*B*(PI/n)*(mod(i,n) - (i / (n**2))))-&
                        skd(mod(i,n**2)/n,mod(j,n**2)/n,-1,n)*exp(cmplx(0,-1)*B*(PI/n)*(n/2 - mod(i,n) - (n/2 - i / (n**2)))    ))*&
                        skd(mod(i,n),mod(j,n),0,n)*skd(i/(n**2),j/(n**2),0,n)*cmplx(0,-0.5)
                    end do 
                end do 
            end if
        else 
            if(sign == 1) then
                do i = 0,(n**3)-1
                    do j = 0,(n**3)-1
                        matrix(i+1,j+1) = (skd(mod(i,n),mod(j,n),1,n)*exp(cmplx(0,1)*B*(PI/n)*(-(i / (n**2))-(mod(i,n**2)/n)))+&
                        skd(mod(i,n),mod(j,n),-1,n)*exp(cmplx(0,-1)*B*(PI/n)*(-(n/2 - i / (n**2))-(n/2 - mod(i,n**2)/n))))*&
                        skd((i/(n**2)),(j/(n**2)),0,n)*skd(mod(i,n**2)/n,mod(j,n**2)/n,0,n)*0.5
                    end do 
                end do 
            else
                do i = 0,(n**3)-1
                    do j = 0,(n**3)-1
                        matrix(i+1,j+1) = (skd(mod(i,n),mod(j,n),1,n)*exp(cmplx(0,1)*B*(PI/n)*(-(i / (n**2))-(mod(i,n**2)/n)))-&
                        skd(mod(i,n),mod(j,n),-1,n)*exp(cmplx(0,-1)*B*(PI/n)*(-(n/2 - i / (n**2))-(2 - mod(i,n**2)/n))))*&
                        skd((i/(n**2)),(j /(n**2)),0,n)*skd(mod(i,n**2)/n,mod(j,n**2)/n,0,n)*cmplx(0,-0.5)
                    end do 
                end do 
            end if
        end if 
        end function mI1

end program magneticfield
