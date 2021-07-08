program newton

        real:: a ,f ,df, epsi, x

        write(*,*) "find the cube root of:"
        read(*,*) a
        write(*,*) "initial x0?:"
        read(*,*) x
        write(*,*) "tolerance ?:"
        read(*,*) epsi

        i=0
        do
              !change your func
                f = x**3 - a  !function 
                df = 3*x**2  ! derivative
                
                xn = x - f / df ! the next value 
                x = xn
                i = i + 1
                write(*,*)"step", i, "x = ", xn, "f(x) =", f, "f'(x)=", df

                !stop the iteration when it reaches abs(f)< epsi
                if (abs(f)< epsi) then
                        exit
                end if
        end do

        write(*,200)'the cube root of',a,'is',xn
        200  FORMAT(A16,F6.2,A3,F9.6) ! formating the output 

end program newton 