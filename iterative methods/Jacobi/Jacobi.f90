program Jacobi

implicit none

real ( kind = 8 ), allocatable :: A(:,:), D(:,:),inv_D(:,:),U(:,:),L(:,:),F(:,:)
real ( kind = 8 ), allocatable :: b(:)
real ( kind = 8 ), allocatable :: x(:),err(:),xold(:)
real :: bi,ai,epsi,max,sum1
integer ::  n,i,j,k,itr

write(*,*)'matrix dimension n = '
read(*,*)n

allocate(A(n,n),b(n),x(n),xold(n),err(n),D(n,n),inv_D(n,n),U(n,n),L(n,n),F(n,n))

! assign values to matrix
!input matrix b 
do i  =1,n
   write(*,"(A10,I2)")'value of b',i
   read(*,*)bi
   b(i) = bi
end do
! formating b into a matrix output
write(*,*)"b = "
do i=1,size(b,1)
 write(*,'(20G12.4)') b(i)
end do

!input matrix A
do i=1,n
      do j = 1, n
         write(*,"(A10,I2,I2)")'value of a',i,j
          read(*,*)ai
         A(i, j) = ai
      end do
end do
! formating A into a matrix output
write(*,*)"A = "
do i=1,size(A,1)
 write(*,'(20G12.4)') A(i,:)
end do

! initial guess x
x=0.0

write(*,*)'number of iteration'
read(*,*)itr
write(*,*)'Tolerance'
read(*,*)epsi

! Matrix manipulation
! Diagonal of Matrix A
do i=1,n
	do j=1,n
		if(i==j) then
			D(i,j)=A(i,j)
		else
			D(i,j)=0
		endif
	enddo 
enddo

! Lower Triangle Matrix L
do i=1,n
	do j=1,n
		if(i==j) then
                        L(i,j)=0
                else
                        L(i,j)=A(i,j)
                endif
	enddo 
enddo
do i=1,n-1
	do j=i+1,n
		L(i,j)=0 
	enddo 
enddo


! Upper Triangle Matrix U
do i=1,n
        do j=1,n
                if(i==j) then
                        U(i,j)=0  ! diagonal is 0
                else
                        U(i,j)=A(i,j) ! the rest is equal to A
                endif
        enddo
enddo
do i=1,n
	do j=i+1,n
                U(j,i)=0
	enddo
enddo

!inverse of D
do i=1,n
	do j=1,n
                if(i==j) then
                        inv_D(i,j)=1/D(i,j)
                else
                        inv_D(i,j)=0
                endif
        enddo
enddo

!solving Ax=b

k=0
do i=1,itr
	F=-(L+U)
	x= matmul(matmul(F,x)+b,inv_D) ! matmul for matrix multiplication
        err = matmul(A,x)-b
        max = sum(err)/size(err)! get the average value of an array
        k=k+1
	if (abs(max)<epsi)then
		exit
	endif
end do
write(*,"(A3,I3,A26)")'for ',k,' iteration the solution is '
! formating the output of x into a column
do i=1,size(x,1)
	write(*,'(A1,I1,A1,20G12.4)')'x',i,"=", x(i)
end do

!method 2 without the matrices
do k=1,itr
        xold = x
        do i=1,n
                sum1 = 0
                do j=1,n
                        if (i/=j) then
                                sum1 = sum1 + A(i,j)*xold(j)
                        endif
                enddo
                x(i) = (1/A(i,i))*(b(i)-sum1)
        enddo
        err = matmul(A,x)-b
        max = sum(err)/size(err)! get the minimum value of an array
        f=f+1
        if (abs(max)<epsi)then
                exit
        endif
enddo
 
End Program Jacobi
