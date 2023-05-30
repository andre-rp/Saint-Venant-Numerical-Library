!   Método MUSCL-Hancock 
! 
!   Solução das equações de Saint-Venant utilizando esquema MUSCL para recons_
!   trução linear e o passo de Hancock para evoluir 0.5*dt, para depois 
!   avançar a solução no tempo.
!   
!   Os testes suportados por esse programa são: Escoamentos subcrítico, 
!   transcrítico e leito seco.

program db_MUSCLH
    implicit none
    
    real             :: tic, toc
    double precision :: t, tmax
    integer          :: k, ntmax
    
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
        
!       Reconstrução linear MUSCL
        call MUSCL

!       Calculo de vetor de fluxo numerico com o solver HLL
        call HLLRS

!       Avanço da solução
        call update

!       Acrescenta os termos fontes (splitting technique)
        call source
                               
!       Aplicação das condições de fronteira (fronteira simetrica)
        call boundary
      
        if (t.GE.tmax) then            
!           Escrita das informações ao final do tempo maximo de execução        
            call output
            
!           Executa o plot da solução no tempo final
            call execute_command_line('gnuplot -p plot_db_MUSCLH.plt')
        
            call cpu_time(toc)
            
            write(*,*) '--------------------------------'           
            write(*,*) 'Numero de iterações = ', k
            write(*,*) 'Tempo de execução   = ', toc-tic, 'segundos'
            write(*,*) '--------------------------------' 
            exit
        end if
    end do

end program db_MUSCLH
!--------------------------------------------------------------------
!Leitura dos dados de entrada
!--------------------------------------------------------------------
subroutine leitura
    implicit none
    
    double precision     :: cr, g, dx, dt, xmin, xmax, epsilon
    double precision     :: b, n, hu, hd, vu, vd, t, tmax, dam 
    integer :: m, ntmax, limiter 
    
    common /cons/ n, cr, b
    common /limi/ epsilon, limiter
    common /domi/ xmin, xmax, dam
    common /cini/ hu, hd, vu, vd    
    common /temp/ t, tmax, ntmax
    common /malh/ dx, dt, m
    common /grav/ g
    
    open(unit = 1, file = 'db_MUSCLH.ini', status = 'unknown')
    
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
    read(1,*) epsilon   ! Limitador do leito seco
    read(1,*) limiter   ! Escolha do limitador (0, 1, 2 ou 3)
    read(1,*) ntmax     ! Numero maximo de iterações no tempo
    read(1,*) m         ! Numero de pontos/celulas no dominio
    read(1,*) g         ! Aceleração da gravidade
    
    close(1)
    
    write(*,*) '--------------------------------'
    write(*,*) 'Dados iniciais:'
    write(*,*) '--------------------------------'
    write(*,*) 'xmin    = ', xmin
    write(*,*) 'xmax    = ', xmax
    write(*,*) 'dam     = ', dam
    write(*,*) 'b       = ', b
    write(*,*) 'hu      = ', hu
    write(*,*) 'hd      = ', hd
    write(*,*) 'vu      = ', vu
    write(*,*) 'vd      = ', vd
    write(*,*) 'tmax    = ', tmax
    write(*,*) 'Manng   = ', n
    write(*,*) 'CFL     = ', cr
    write(*,*) 'epsilon = ', epsilon
    write(*,*) 'limiter = ', limiter
    write(*,*) 'ntmax   = ', ntmax
    write(*,*) 'm       = ', m
    write(*,*) 'g       = ', g
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
    double precision, dimension (2,-1:mx+1) :: U, F, S
    double precision, dimension (-1:mx+1) :: x, zb, zL, zR
    
    common /domi/ xmin, xmax, dam
    common /cini/ hu, hd, vu, vd
    common /malh/ dx, dt, m
    common /vari/ U, x, F, S
    common /varb/ zb, zL, zR
    common /grav/ g
    
    dx = (xmax-xmin)/dble(m) 
    
    x(0) = xmin
    do i = 1,m
        x(i)  = x(i-1)+dx
        
        zL(i) = 0
        zR(i) = 0
        zb(i) = 0.5*(zL(i)+zR(i))
        
        if (x(i).LE.dam) then
            U(1,i) = hu
            U(2,i) = hu*vu
        else
            U(1,i) = hd
            U(2,i) = hd*vd
        end if
    end do
    
    U(:,0) = U(:,1)
    U(:,m+1) = U(:,m)

end subroutine start
!--------------------------------------------------------------------
!Atualização do passo de tempo
!--------------------------------------------------------------------
subroutine CFLcon(k)
    implicit none

    double precision     :: cr, g, dx, dt, n, b, t, tmax, epsilon
    double precision     :: speloc = 0, smax = 0  
    integer :: i, m, k, ntmax, limiter 
    integer, parameter :: mx = 5000
    double precision, dimension (2,-1:mx+1) :: U, F, S
    double precision, dimension (-1:mx+1) :: x
    
    common /limi/ epsilon, limiter
    common /cons/ n, cr, b
    common /temp/ t, tmax, ntmax
    common /malh/ dx, dt, m
    common /vari/ U, x, F, S
    common /grav/ g

    do i = 1,m
        if (U(1,i).LE.epsilon) then
            speloc = 0
        else
            speloc = abs(U(2,i)/U(1,i))+sqrt(g*U(1,i))
        end if
        
        if (speloc.GT.smax) then
            smax = speloc
        end if
    end do   
    
    dt = cr*dx/smax 

    if (k.LE.5) then
        dt = 0.2*dt
    end if

    if ((t+dt).GT.tmax) then
        dt = tmax-t
    end if
    
    t = t+dt

end subroutine CFLcon
!--------------------------------------------------------------------
!MUSCL
!--------------------------------------------------------------------
subroutine MUSCL
    implicit none
    
    double precision     :: g, dx, dt, epsilon
    double precision     :: dup, ddw, delta, mino
    integer :: i, j, m, limiter 
    integer, parameter :: mx = 5000
    double precision, dimension (2,2,-1:mx+1) :: BEXT
    double precision, dimension (2,-1:mx+1) :: U, F, S
    double precision, dimension (-1:mx+1) :: x, zb, zR, zL
    double precision, dimension (2) :: FL, FR, uL, uR, Si
    
    common /limi/ epsilon, limiter
    common /vari/ U, x, F, S
    common /varb/ zb, zL, zR
    common /malh/ dx, dt, m
    common /extr/ BEXT
    common /grav/ g

    do i = 1,m   
        do j = 1,2
            dup = U(j,i)-U(j,i-1)
            ddw = U(j,i+1)-U(j,i)

!           Limitador de inclinacão                      
            if ((dup*ddw).LE.0) then
                delta = 0
            else            
                mino = min(abs(ddw),abs(dup))
                delta = sign(dble(1),ddw)*min(mino,limiter*mino)
            end if

!           Variaves (h,q) extrapoladas para as fronteiras a esquerda e a
!           direita                       
            uL(j) = U(j,i)-0.5*delta 
            uR(j) = U(j,i)+0.5*delta
        end do

        call flux(FL,uL)
        call flux(FR,uR)
        call msor(Si,U(:,i),uL,uR,zR(i),zL(i))
        
!       Passo de Hancock para evoluir a solucao dt/2        
        do j = 1,2
            BEXT(1,j,i) = uL(j)+0.5*(dt/dx)*(FL(j)-FR(j))+0.5*dt*Si(j)
            BEXT(2,j,i) = uR(j)+0.5*(dt/dx)*(FL(j)-FR(j))+0.5*dt*Si(j)
        end do
        
        if (BEXT(1,1,i).LT.0) then
                BEXT(1,1,i) = 0
                BEXT(1,2,i) = 0
        end if
        if (BEXT(2,1,i).LT.0) then
                BEXT(2,1,i) = 0
                BEXT(2,2,i) = 0
        end if

        if (BEXT(1,1,i).LT.epsilon) then
                BEXT(1,2,i) = 0
        end if
        if (BEXT(2,1,i).LT.epsilon) then
                BEXT(2,2,i) = 0
         end if         
    end do

end subroutine MUSCL
!--------------------------------------------------------------------
!Fluxo numerico utilizado no MUSCL
!--------------------------------------------------------------------
subroutine flux(FK,Uk)
    implicit none

    double precision :: g, hk, vk, epsilon, limiter
    double precision, dimension (2) :: FK, Uk
    
    common /limi/ epsilon, limiter
    common /grav/ g

    hk = Uk(1)
    if (hk.LE.epsilon) then
        FK(1) = 0
        FK(2) = 0 
    else
        vk = Uk(2)/Uk(1)
        FK(1) = hk*vk
        FK(2) = hk*vk**2+0.5*g*hk**2
    end if

end subroutine flux
!--------------------------------------------------------------------
!Termos fontes utilizado no MUSCL
!--------------------------------------------------------------------
subroutine msor(Si,Ui,uL,uR,zR,zL)
    implicit none

    double precision :: g, hk, vk, sf, Rh, zR, zL, n, cr, b, dx, dt, m
    double precision :: epsilon, limiter
    double precision, dimension (2) :: Si, Ui, uL, uR
    integer :: i
    
    common /limi/ epsilon, limiter
    common /malh/ dx, dt, m
    common /cons/ n, cr, b
    common /grav/ g

    hk = Ui(1)
        
    if (hk.LE.epsilon) then                 
        Si(1) = 0
        Si(2) = 0
    else
        vk = Ui(2)/hk
        
        Rh = b*hk/(b+2*hk)
        sf = (abs(vk)*vk*n**2)/Rh**(4./3.)
         
        Si(1) = 0
        Si(2) = -g*hk*sf-g*(uL(1)+uR(1))*0.5*(zR-zL)/dx
    end if

end subroutine msor
!--------------------------------------------------------------------
!HLLRS
!--------------------------------------------------------------------
subroutine HLLRS
    implicit none
    
    double precision     :: g, dx, dt, epsilon
    double precision     :: hL, hR, hs, vL, vR, cL, cR, sL, sR, iswet
    integer :: i, j, m, limiter 
    integer, parameter :: mx = 5000
    double precision, dimension (2,2,-1:mx+1) :: BEXT
    double precision, dimension (2,-1:mx+1) :: U, F, S
    double precision, dimension (-1:mx+1) :: x
    double precision, dimension (2) :: FL, FR, FHLL, CDL, CDR
    
    common /limi/ epsilon, limiter
    common /vari/ U, x, F, S
    common /malh/ dx, dt, m
    common /extr/ BEXT
    common /grav/ g
    
    do i = 1,m-1
        do j = 1,2
            CDL(j) = BEXT(2,j,i)
            CDR(j) = BEXT(1,j,i+1)
        end do
        
        iswet = 1
        
        hL = CDL(1) 
        cL = sqrt(g*hL)
        if (hL.LE.epsilon) then
            vL = 0
        else
            vL = CDL(2)/hL
        end if
        
        hR = CDR(1) 
        cR = sqrt(g*hR)
        if (hR.LE.epsilon) then
            vR = 0
        else
            vR = CDR(2)/hR
        end if
        
        call flux(FL,CDL)
        call flux(FR,CDR)
       
!       Velocidade da crista e da cauda para as ondas no caso de leito 
!       seco        
        if (hR.LE.epsilon.AND.hL.GT.epsilon) then
!           Leito seco à direita 
            sL = vL-cL
            sR = vL+2*cL
            iswet = 0
        end if
        if (hL.LE.epsilon.AND.hR.GT.epsilon) then
!           Leito seco à esquerda         
            sL = vR-2*cR
            sR = vR+CR
            iswet = 0
        end if
        if (hL.LE.epsilon.AND.hR.LE.epsilon) then
!           Leito seco
            sL = -epsilon
            sR = epsilon
            iswet = 0
        end if
        
!       Velocidade da crista e da cauda para as ondas no caso de leito 
!       molhado        
        if (iswet.EQ.1) then
            hs = ((0.5*(cL+cR)+0.25*(vL-vR))**2)/g
            
            if (hs.LE.hL) then
                sL = vL-cL
            else
                sL = vL-(cL*sqrt(0.5*hs*(hs+hL)))/hL
            end if 
            
            if (hs.LE.hR) then
                sR = vR+cR
            else
                sR = vR+(cR*sqrt(0.5*hs*(hs+hR)))/hR
            end if
        end if
        
!       Calculo fluxo na interface (i,i+1) em x/t=0        
        if (sL.GE.0) then
            F(:,i) = FL
        end if
        
        if (sL.LE.0.AND.sR.GE.0) then
            FHLL = sR*FL-sL*FR+sL*sR*(CDR-CDL)
            F(:,i) = FHLL/(sR-sL)  
        end if
        
        if (sR.LE.0) then
            F(:,i) = FR
        end if    
    end do 
   
end subroutine HLLRS
!--------------------------------------------------------------------
!Atualização da solução
!--------------------------------------------------------------------
subroutine update
    implicit none
    
    double precision :: dx, dt, epsilon, limiter  
    integer :: i, m
    integer, parameter :: mx = 5000
    double precision, dimension (2,-1:mx+1) :: U, F, S
    double precision, dimension (-1:mx+1) :: x
    
    common /limi/ epsilon, limiter
    common /vari/ U, x, F, S
    common /malh/ dx, dt, m
    
    do i = 2,m-1
        U(:,i) = U(:,i)-(dt/dx)*(F(:,i)-F(:,i-1))
 
!       Conservacao da positividade na profundidade       
        if (U(1,i).LT.0) then
            U(:,i) = 0
        end if
        
!       Evitar altas velocidade        
        if (U(1,i).LE.epsilon) then
            U(2,i) = 0
        end if
    end do

end subroutine update
!--------------------------------------------------------------------
!Termos fontes
!--------------------------------------------------------------------
subroutine source
    implicit none
    
    double precision :: g, n, cr, b, dx, dt, epsilon, limiter
    double precision :: hk, vk, Rh, sf  
    integer :: i, m
    integer, parameter :: mx = 5000
    double precision, dimension (2,2,-1:mx+1) :: BEXT
    double precision, dimension (2,-1:mx+1) :: U, F, S
    double precision, dimension (-1:mx+1) :: x, zb, zL, zR
    
    common /limi/ epsilon, limiter
    common /vari/ U, x, F, S
    common /varb/ zb, zL, zR
    common /malh/ dx, dt, m
    common /cons/ n, cr, b
    common /extr/ BEXT
    common /grav/ g
    
    do i = 2,m-1
        hk = U(1,i)
        
        if (hk.LE.epsilon) then  
            S(1,i) = 0
            S(2,i) = 0
        else
            vk = U(2,i)/hk
            
            Rh = b*hk/(b+2*hk)
            sf = (abs(vk)*vk*n**2)/Rh**(4./3.)
            
            S(1,i) = 0
            S(2,i) = -g*hk*sf-g*(BEXT(1,1,i)+BEXT(2,1,i))*0.5*(zR(i)-zL(i))/dx
        end if
        
        U(2,i) = U(2,i)+dt*S(2,i)
    end do
        

end subroutine source
!--------------------------------------------------------------------
!Condição de fronteira
!--------------------------------------------------------------------
subroutine boundary
    implicit none
    
    double precision :: dx, dt
    integer :: m
    integer, parameter :: mx = 5000
    double precision, dimension (2,-1:mx+1) :: U, F, S
    double precision, dimension (-1:mx+1) :: x
    
    common /vari/ U, x, F, S
    common /malh/ dx, dt, m
    
    U(1,1) = U(1,2)
    U(2,1) = 0
    U(:,0) = U(:,1)
    
    U(1,m) = U(1,m-1)
    U(2,m) = 0   
    U(:,m+1) = U(:,m)
    
end subroutine boundary
!--------------------------------------------------------------------
!Geração do arquivo de saída
!--------------------------------------------------------------------
subroutine output
    implicit none
    
    double precision  :: g, dx, dt, epsilon
    integer :: i, m, limiter 
    integer, parameter :: mx = 5000
    double precision, dimension (2,-1:mx+1) :: U, F, S
    double precision, dimension (-1:mx+1) :: x, zb, zL, zR

    common /limi/ epsilon, limiter
    common /malh/ dx, dt, m
    common /vari/ U, x, F, S
    common /varb/ zb, zL, zR
    common /grav/ g 
    
    open(unit = 1, file = 'db_MUSCLH.out', status = 'unknown')
    do i = 1,m
        if (U(1,i).LE.epsilon) then        
            write(1,*) x(i), U(1,i), U(2,i), 0
        else
            write(1,*) x(i), U(1,i), U(2,i), U(2,i)/U(1,i)
        end if
    end do

    close(1)
     
end subroutine output
