Program gSeidel

implicit none
real ( kind = 8 ), allocatable :: A(:,:)
real ( kind = 8 ), allocatable :: b(:)
real ( kind = 8 ), allocatable :: x(:),err(:),xold(:)
real :: bi,ai,epsi,max,sum1,sum2
integer ::  n,i,j,k,itr,f

write(*,*)'matrix dimension n = '
read(*,*)n

allocate(A(n,n),b(n),x(n),err(n),xold(n))

! assign values to matrix
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

write(*,*)'number of iteration'
read(*,*)itr

write(*,*)'Tolerance'
read(*,*)epsi

! initial guess x
x=0.0
! solving Ax = b
f=0
do k=1,itr
        xold = x
        do i=1,n
                sum1 = 0
                do j=1,i-1
                        sum1 = sum1 + A(i,j)*x(j)
                end do
                sum2 = 0
                do j=i+1,n
                        sum2 = sum2 + A(i,j)*xold(j)
                enddo
                x(i) = (1/A(i,i))*(b(i)-sum1-sum2)
        enddo
        err = matmul(A,x)-b
        max = sum(err)/size(err)! get the average value of an array
        f=f+1
        if (abs(max)<epsi)then
                exit
        endif
enddo

write(*,"(A3,I3,A26)")'for ',f,' iteration the solution is '
! formating the output of x into a column
do i=1,size(x,1)
        write(*,'(A1,I1,A1,20G12.4)')'x',i,"=", x(i)
end do

End Program gSeidel
