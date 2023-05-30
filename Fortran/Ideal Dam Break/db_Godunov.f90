!   Método Godunov de primeira ordem 
! 
!   Solucao das equacoes de Saint-Venant utilizando esquema de Godunov 
!   de primeira ordem, que utiliza o solver de Riemann para determinar  
!   o fluxo numerico e depois avancar a solucao no tempo.
!   
!   Os testes suportados por esse programa são: Escoamentos subcrítico, 
!   transcrítico e leito seco.

program db_Godunov
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

!       Calculo de vetor de fluxo numerico com o solver exato de Riemann
        call GodunovRs

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
            call execute_command_line('gnuplot -p plot_db_Godunov.plt')
        
            call cpu_time(toc)
            
            write(*,*) '--------------------------------'           
            write(*,*) 'Numero de iterações = ', k
            write(*,*) 'Tempo de execução   = ', toc-tic, 'segundos'
            write(*,*) '--------------------------------' 
            exit
        end if
    end do

end program db_Godunov
!--------------------------------------------------------------------
!Leitura dos dados de entrada
!--------------------------------------------------------------------
subroutine leitura
    implicit none
    
    double precision     :: cr, g, dx, dt, xmin, xmax, epsilon, tol
    double precision     :: b, n, hu, hd, vu, vd, t, tmax, dam 
    integer :: m, ntmax, nint 
    
    common /cons/ n, cr, b
    common /limi/ epsilon
    common /domi/ xmin, xmax, dam
    common /cini/ hu, hd, vu, vd    
    common /temp/ t, tmax, ntmax
    common /tole/ tol, nint
    common /malh/ dx, dt, m
    common /grav/ g
    
    open(unit = 1, file = 'db_Godunov.ini', status = 'unknown')
    
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
    read(1,*) ntmax     ! Numero maximo de iterações no tempo
    read(1,*) nint      ! Numero maximo de iterações do metodo de Newton
    read(1,*) tol       ! Tolêrancia do metodo de Newton 
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
    write(*,*) 'ntmax   = ', ntmax
    write(*,*) 'nint    = ', nint
    write(*,*) 'tol     = ', tol
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
    integer :: i, m, k, ntmax 
    integer, parameter :: mx = 5000
    double precision, dimension (2,-1:mx+1) :: U, F, S
    double precision, dimension (-1:mx+1) :: x
    
    common /limi/ epsilon
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
!Fluxo Numerico utilizando os solver exatos do problema de Riemann
!--------------------------------------------------------------------
subroutine GodunovRs
    implicit none
    
    double precision :: g, dx, dt, hL, hR, vL, vR, cL, cR
    integer :: i, m
    integer, parameter :: mx = 5000
    double precision, dimension (2,-1:mx+1) :: U, F, S
    double precision, dimension (-1:mx+1) :: x
    
    common /ctes/ hL, hR, vL, vR, cL, cR
    common /vari/ U, x, F, S
    common /malh/ dx, dt, m
    common /grav/ g

    do i = 1,m-1   
        hL = U(1,i)
        cL = sqrt(g*hL)
        if (hL.LE.0) then
            hL = 0
            vL = 0
        else
            vL = U(2,i)/U(1,i)
        end if
        
        hR = U(1,i+1)
        cR = sqrt(g*hR)
        if (hR.LE.0) then
            hR = 0
            vR = 0
        else
            vR = U(2,i+1)/U(1,i+1)
        end if
        
        call RiemannSolver(F(:,i))
    end do
    
end subroutine GodunovRs
!--------------------------------------------------------------------
!Determinação do problema de Riemann
!--------------------------------------------------------------------
subroutine RiemannSolver(Fi)
    implicit none
    
    double precision :: h, v, hL, hR, vL, vR, cL, cR, cs
    double precision :: g, epsilon
    double precision, dimension (2) :: Fi
    
    common /limi/ epsilon
    common /ctes/ hL, hR, vL, vR, cL, cR
    common /grav/ g
    
    cs = (vR-vL)-2*(cL+cR)
    
    if (hL.LE.epsilon.OR.hR.LE.epsilon.OR.cs.GE.0) then
        call drybed(h,v)
    else
        call wetbed(h,v)
    end if
    
    Fi(1) = h*v
    Fi(2) = h*v**2+0.5*g*h**2    
    
end subroutine RiemannSolver
!--------------------------------------------------------------------
!Solver de Riemann para leito seco (drybed)
!--------------------------------------------------------------------
subroutine drybed(h,v)
    implicit none
    
    double precision :: h, v, s = 0
    double precision :: hL, hR, vL, vR, cL, cR, epsilon
    
    common /ctes/ hL, hR, vL, vR, cL, cR
    common /limi/ epsilon
        
    if (hL.LE.epsilon) then
        call drylef(h,v,s)
    else
        if (hR.LE.epsilon) then
            call dryrig(h,v,s)
        else
            call drymid(h,v,s)
        end if
    end if      
        
end subroutine drybed
!--------------------------------------------------------------------
!Construção da solução para o caso de leito seco à esquerda (dry left, hL = 0)
!--------------------------------------------------------------------
subroutine drylef(h,v,s)
    implicit none
    
    double precision :: g, h, v, s, sR, sRs
    double precision :: hL, hR, vL, vR, cL, cR
    
    common /ctes/ hL, hR, vL, vR, cL, cR
    common /grav/ g
    
    sR = vR+cR
    if (s.GE.sR) then
        h = hR
        v = vR
    else
        sRs = vR-2*cR
        if (s.GE.sRs) then
            h = ((-vR+2*cR+s)**2)/(9*g)
            v = (vR-2*cR+2*s)/3
        else
            h = hL
            v = vL
        end if
    end if
    
end subroutine drylef
!--------------------------------------------------------------------
!Construção da solução para o caso de leito seco à direita (dry right, hR = 0)
!--------------------------------------------------------------------
subroutine dryrig(h,v,s)
    implicit none
    
    double precision :: g, h, v, s, sL, sLs
    double precision :: hL, hR, vL, vR, cL, cR
    
    common /ctes/ hL, hR, vL, vR, cL, cR
    common /grav/ g
    
    sL = vL-cL
    if (s.LE.sL) then
        h = hL
        v = vL
    else
        sLs = vL+2*cL
        if (s.LE.sLs) then
            h = (vL+2*cL-s)**2/(9*g)
            v = (vL+2*cL+2*s)/3
        else
            h = hR
            v = vR
        end if
    end if
    
end subroutine dryrig
!--------------------------------------------------------------------
!Construção da solução para o caso de leito seco no meio (hs = 0)
!--------------------------------------------------------------------
subroutine drymid(h,v,s)
    implicit none
    double precision :: g, h, v, s, sL, sR, sLs, sRs
    double precision :: hL, hR, vL, vR, cL, cR
    
    common /ctes/ hL, hR, vL, vR, cL, cR
    common /grav/ g
    
    sL = vL-cL 
    sR = vR+cR 
    
    sLs = vL+2*cL
    sRs = vR-2*cR

    if (s.LE.sL) then
        h = hL
        v = vL
    end if
        
    if (s.GT.sL.AND.s.LE.sLs) then
        h = (vL+2*cL-s)**2/(9*g)
        v = (vL+2*cL+2*s)/3
    end if
        
    if (s.GT.sLs.AND.s.LE.sRs) then
        h = 0
        v = 0
    end if
        
    if (s.GT.sRs.AND.s.LE.sR) then
        h = (-vR+2*cR+s)**2/(9*g)
        v = (vR-2*cR+2*s)/3
    end if
        
    if (s.GT.sR) then
        h = hR
        v = vR
    end if

end subroutine drymid
!--------------------------------------------------------------------
!Solver de Riemann para leito molhado (wetbed)
!--------------------------------------------------------------------
subroutine wetbed(h,v)
    implicit none

    double precision :: g, tol, h, v, s = 0
    double precision :: f, df, hz, fL, fR, dfL, dfR
    double precision :: hL, hR, hs, vL, vR, vs, cL, cR, cs
    integer :: k, nint 
    
    common /ctes/ hL, hR, vL, vR, cL, cR
    common /star/ hs, vs, cs
    common /tole/ tol, nint
    common /grav/ g

!   Chute inicial para o metodo de Newton-Raphson     
    hs = 1/g*(0.5*(cL+cR)+0.25*(vL-vR))**2
    hz = hs
    
!   Aplicação do método de Newton para determinar as variaveis na região estrela    
    do k = 1,nint
    
        cs = sqrt(g*hs)
        
        call Nflux(fL,dfL,hL)
        call Nflux(fR,dfR,hR)
        
        f  = fR+fL+vR-vL
        df = dfL+dfR
        
        hs = hs-f/df
        
        if (abs(hs-hz).LE.tol) then
            exit
        else
            hz = hs
        end if
    end do
    
    cs = sqrt(g*hs)
    vs = 0.5*(vR+vL+fR-fL)

    call samwet(h,v,s)
  
end subroutine wetbed
!--------------------------------------------------------------------
!Avaliação da função f(h*) pelo método de Newton
!--------------------------------------------------------------------
subroutine Nflux(fK,dfK,hK)
    implicit none
    
    double precision :: g, fK, dfK, hK, cK
    double precision :: hs, vs, cs
    
    common /star/ hs, vs, cs
    common /grav/ g
    
    cK = sqrt(g*hK)
    
    if (hs.LE.hK) then
        fK = 2*(cs-cK)
        dfK = g/cs
    else
        fK = (hs-hK)*sqrt((hs+hK)*g/(2*hs*hK))
        dfK = sqrt((hs+hK)*g/(2*hs*hK))-(hs-hK)*g/(4*hs**2)*1/sqrt((hs+hK)*g/(2*hs*hK))
    end if 
    
end subroutine Nflux
!--------------------------------------------------------------------
!Construção da solução para o caso de leito molhado
!--------------------------------------------------------------------
subroutine samwet(h,v,s)
    implicit none
    
    double precision :: g, h, v, s, lbL, lbR, sL, sLs, sR, sRs
    double precision :: hL, hR, hs, vL, vR, vs, cL, cR, cs
    
    common /ctes/ hL, hR, vL, vR, cL, cR
    common /star/ hs, vs, cs
    common /grav/ g

    if (s.LE.vs) then
!       Onda á esquerda
        if (hL.LE.hs) then
!           Onda de choque
            lbL = sqrt(0.5*(((hs+hL)*hs)/hL**2))
            sL = vL-lbL*cL
            if (s.LE.sL) then 
                h = hL
                v = vL
            else
                h = hs
                v = vs
            end if
        else
!           Onda de rarefação
            sL = vL-cL
            if (s.LE.sL) then 
                h = hL
                v = vL
            else
                sLs = vs-cs
                if (s.LE.sLs) then
                    h = ((vL+2*cL-s)**2)/(9*g)
                    v = (vL+2*cL+2*s)/3
                else 
                    h = hs
                    v = vs
                end if
            end if
        end if    
    else
!       Onda á direita
        if (hR.LE.hs) then
!           Onda de choque        
            lbR =  sqrt(0.5*(((hs+hR)*hs)/hR**2))
            sR = vR+lbR*cR
            if (s.GE.sR) then 
                h = hR
                v = vR
            else
                h = hs
                v = vs
            end if
        else
!           Onde de rarefação        
            sR = vR+cR
            if (s.GE.sR) then 
                h = hR
                v = vR
            else
                sRs = vs+cs
                if (s.GE.sRs) then
                    h = ((-vR+2*cR+s)**2)/(9*g)
                    v = (vR-2*cR+2*s)/3
                else 
                    h = hs
                    v = vs
                end if
            end if
        end if
    end if
end subroutine samwet
!--------------------------------------------------------------------
!Atualização da solução
!--------------------------------------------------------------------
subroutine update
    implicit none
    
    double precision :: dx, dt  
    integer :: i, m
    integer, parameter :: mx = 5000
    double precision, dimension (2,-1:mx+1) :: U, F, S
    double precision, dimension (-1:mx+1) :: x
    
    common /vari/ U, x, F, S
    common /malh/ dx, dt, m
    
    do i = 2,m-1
        U(:,i) = U(:,i)-(dt/dx)*(F(:,i)-F(:,i-1))
    end do

end subroutine update
!--------------------------------------------------------------------
!Termos fontes
!--------------------------------------------------------------------
subroutine source
    implicit none
    
    double precision :: g, n, cr, b, dx, dt, epsilon
    double precision :: hk, vk, Rh, sf  
    integer :: i, m
    integer, parameter :: mx = 5000
    double precision, dimension (2,-1:mx+1) :: U, F, S
    double precision, dimension (-1:mx+1) :: x, zb, zL, zR
    
    common /limi/ epsilon
    common /vari/ U, x, F, S
    common /varb/ zb, zL, zR
    common /malh/ dx, dt, m
    common /cons/ n, cr, b
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
            S(2,i) = -g*hk*sf-g*hk*(zR(i)-zL(i))/dx
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
    integer :: i, m 
    integer, parameter :: mx = 5000
    double precision, dimension (2,-1:mx+1) :: U, F, S
    double precision, dimension (-1:mx+1) :: x, zb, zL, zR

    common /limi/ epsilon
    common /malh/ dx, dt, m
    common /vari/ U, x, F, S
    common /varb/ zb, zL, zR
    common /grav/ g 
    
    open(unit = 1, file = 'db_Godunov.out', status = 'unknown')
    do i = 1,m
        if (U(1,i).LE.epsilon) then        
            write(1,*) x(i), U(1,i), U(2,i), 0
        else
            write(1,*) x(i), U(1,i), U(2,i), U(2,i)/U(1,i)
        end if
    end do

    close(1)
     
end subroutine output
