!   Método de Lax difusivo 
! 
!   Solução das equações de Saint-Venant utilizando o método de Lax difusivo 
!   para avançar a solução no tempo.
!   
!   Os testes suportados por esse programa são: Escoamentos subcrítico e 
!   transcrítico.

program DB_Lax
    implicit none
    
    real                :: tic, toc
    double precision    :: t, tmax
    integer :: k, ntmax
    
    common /temp/ t, tmax, ntmax
    
    call cpu_time(tic)
    
!   Leitura do arquivo .ini que contém as condições iniciais 
    call leitura

!   Inicia as variaveis que serão utilizadas ao longo do programa x, zb, U
    call start 

!   Loop temporal para avançar a solução no tempo
    do  k = 1,ntmax
        
!       Atualização do passo de tempo
        call CFLcon(k)
        
!       Aplicação do método de Lax difusivo para avançar a solução
        call LaxD
               
!       Aplicação das condições de fronteira (fronteira simetrica)
        call boundary
      
        if (t.GE.tmax) then
!           Escrita das informações ao final do tempo maximo de execução        
            call output
            
!           Executa o plot da solução no tempo final
            call execute_command_line('gnuplot -p plot_db_Lax.plt')
        
            call cpu_time(toc)
            
            write(*,*) '--------------------------------'           
            write(*,*) 'Numero de iterações = ', k
            write(*,*) 'Tempo de execução   = ', toc-tic, 'segundos'
            write(*,*) '--------------------------------'
            exit
        end if        
    end do

end program DB_Lax
!--------------------------------------------------------------------
!Leitura dos dados de entrada
!--------------------------------------------------------------------
subroutine leitura
    implicit none
    
    double precision     :: alpha, cr, g, dx, dt, xmin, xmax
    double precision     :: b, n, hu, hd, vu, vd, t, tmax, dam 
    integer :: m, ntmax
    
    common /cons/ n, cr, b, alpha
    common /domi/ xmin, xmax, dam
    common /cini/ hu, hd, vu, vd    
    common /temp/ t, tmax, ntmax
    common /malh/ dx, dt, m
    common /grav/ g
    
    open(unit = 1, file = 'db_Lax.ini', status = 'unknown')
    
    read(1,*) xmin      ! Ponto inicial do canal (Teste subcrítico)
    read(1,*) xmax      ! Ponto final do canal   
    read(1,*) dam       ! Posição da barragem 
    read(1,*) b         ! Largura do canal
    read(1,*) hu        ! Altura da lamina d'água à esquerda da barragem   
    read(1,*) hd        ! Altura da lamina d'água à direita da barragem       
    read(1,*) vu        ! Velocidade à esquerda da barragem
    read(1,*) vd        ! Velocidade à direita da barragem
    read(1,*) tmax      ! Tempo maximo de execução
    read(1,*) n         ! Coeficiente de rugosidade de Manning
    read(1,*) cr        ! Numero de Courant
    read(1,*) alpha     ! Coeficiente de difusão do método de Lax difusivo
    read(1,*) ntmax     ! Numero maximo de iterações no tempo
    read(1,*) m         ! Numero de pontos/celulas no dominio
    read(1,*) g         ! Aceleração da gravidade
    
    close(1)
    
    write(*,*) '--------------------------------'
    write(*,*) 'Dados iniciais:'
    write(*,*) '--------------------------------'
    write(*,*) 'xmin  = ', xmin
    write(*,*) 'xmax  = ', xmax
    write(*,*) 'dam   = ', dam
    write(*,*) 'b     = ', b
    write(*,*) 'hu    = ', hu
    write(*,*) 'hd    = ', hd
    write(*,*) 'vu    = ', vu
    write(*,*) 'vd    = ', vd
    write(*,*) 'tmax  = ', tmax
    write(*,*) 'Manng = ', n
    write(*,*) 'CFL   = ', cr
    write(*,*) 'alpha = ', alpha
    write(*,*) 'ntmax = ', ntmax
    write(*,*) 'm     = ', m
    write(*,*) 'g     = ', g
    write(*,*) '--------------------------------'
    
end subroutine leitura
!--------------------------------------------------------------------
!Inicialização do vetor de variaveis x, zb, U
!--------------------------------------------------------------------
subroutine start
    implicit none
    
    double precision     :: g, dx, dt, xmin, xmax
    double precision     :: hu, hd, vu, vd, dam 
    integer :: i, m
    integer, parameter :: mx = 5000
    double precision, dimension (2,-1:mx+1) :: U
    double precision, dimension (-1:mx+1) :: x, zb
    
    common /domi/ xmin, xmax, dam
    common /cini/ hu, hd, vu, vd
    common /malh/ dx, dt, m
    common /vari/ U, x, zb
    common /grav/ g
    
    dx = (xmax-xmin)/(dble(m-1)) 

    do i = 1,m
        x(i)  = xmin+(i-1)*dx
        zb(i) = 0
        
        if (x(i).LE.dam) then
            U(1,i) = hu
            U(2,i) = hu*vu
        else
            U(1,i) = hd
            U(2,i) = hd*vd
        end if
    end do

end subroutine start
!--------------------------------------------------------------------
!Atualização do passo de tempo
!--------------------------------------------------------------------
subroutine CFLcon(k)
    implicit none

    double precision     :: alpha, cr, g, dx, dt, n, b, t, tmax 
    integer :: m, k, ntmax
    integer, parameter :: mx = 5000
    double precision, dimension (2,-1:mx+1) :: U
    double precision, dimension (-1:mx+1) :: x, zb
    
    common /cons/ n, cr, b, alpha
    common /temp/ t, tmax, ntmax
    common /malh/ dx, dt, m
    common /vari/ U, x, zb
    common /grav/ g

    dt = (1-alpha)*cr*dx*minval(1/(abs(U(2,1:m)/U(1,1:m))+sqrt(g*U(1,1:m))))

    if (k.LE.5) then
        dt = 0.2*dt
    end if

    if ((t+dt).GT.tmax) then
        dt = tmax-t
    end if
    
    t = t+dt

end subroutine CFLcon
!--------------------------------------------------------------------
!Método de Lax difusivo
!--------------------------------------------------------------------
subroutine LaxD
    implicit none
    
    double precision     :: alpha, cr, g, dx, dt, b, n 
    integer :: i, m
    integer, parameter :: mx = 5000
    double precision, dimension (2) :: FL, FR, S 
    double precision, dimension (2,-1:mx+1) :: U, Uz
    double precision, dimension (-1:mx+1) :: x, zb
    
    common /cons/ n, cr, b, alpha
    common /malh/ dx, dt, m
    common /vari/ U, x, zb
    common /grav/ g
    
    Uz = U
    
    do i = 2,m-1
        call flux(FL,Uz(:,i-1))
        call flux(FR,Uz(:,i+1))
        
        call source(S,Uz(:,i))
        
        U(:,i) = alpha*Uz(:,i)+(1-alpha)*0.5*(Uz(:,i+1)+Uz(:,i-1)) &
               -0.5*(dt/dx)*(FR-FL)+dt*S
    end do

end subroutine LaxD
!-------------------------------------------------------------------- 
!Avaliação do vetor de fluxo no ponto U_i no tempo t 
!--------------------------------------------------------------------
subroutine flux(Fk,Uk)
    implicit none
    
    double precision :: g, hk, vk
    double precision, dimension (2) :: Fk, Uk
    
    common /grav/ g
    hk = Uk(1)
    vk = Uk(2)/hk
    
    Fk(1) = hk*vk
    Fk(2) = hk*vk**2+0.5*g*hk**2
    
end subroutine flux
!--------------------------------------------------------------------
!Termo fonte - Equação de Manning + Derivada de zb
!--------------------------------------------------------------------
subroutine source(S,Uk)
    implicit none
    
    double precision :: alpha, cr, g, n, b
    double precision :: hk, vk, Rh, sf
    double precision, dimension (2) :: S, Uk

    common /cons/ n, cr, b, alpha
    common /grav/ g

    hk = Uk(1)
    vk = Uk(2)/hk

    Rh = (b*hk)/(b+2*hk)
    
    sf = (abs(vk)*vk*n**2)/(Rh**(4./3.))
    
    S(1) = 0
    S(2) = -g*hk*sf
    
end subroutine source
!--------------------------------------------------------------------
!Condições de fronteira
!--------------------------------------------------------------------
subroutine boundary
    implicit none
    
    double precision  :: dx, dt
    integer :: m
    integer, parameter :: mx = 5000
    double precision, dimension (2,-1:mx+1) :: U
    double precision, dimension (-1:mx+1) :: x, zb
    
    common /malh/ dx, dt, m
    common /vari/ U, x, zb
    
    U(1,1) = U(1,2)
    U(2,1) = 0

    U(1,m) = U(1,m-1)
    U(2,m) = 0


end subroutine boundary
!--------------------------------------------------------------------
!Geração do arquivo de saída
!--------------------------------------------------------------------
subroutine output
    implicit none
    
    double precision  :: g, dx, dt
    integer :: i, m
    integer, parameter :: mx = 5000
    double precision, dimension (2,-1:mx+1) :: U
    double precision, dimension (-1:mx+1) :: x, zb

    common /malh/ dx, dt, m
    common /vari/ U, x, zb
    common /grav/ g 
    
    open(unit = 1, file = 'db_Lax.out', status = 'unknown')
    do i = 1,m
        if (U(1,i).GT.0) then        
            write(1,*) x(i), U(1,i), U(2,i), U(2,i)/U(1,i)
        else
            write(1,*) x(i), U(1,i), U(2,i), 0
        end if
    end do

    close(1)
     
end subroutine output
