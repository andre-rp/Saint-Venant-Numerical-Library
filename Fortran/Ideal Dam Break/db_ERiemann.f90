!   Solver Exato do problema de Riemann 
! 
!   Solução exata do problema de Riemann por amostragem, baseado na velocidade
!   das ondas que compõe o problema (choque ou rarefação)
! 
!   Os testes suportados por esse programa são: Escoamentos subcrítico, 
!   transcrítico e leito seco.

program DB_ERiemann
    implicit none
    
    real :: tic, toc
    
    call cpu_time(tic)
    
!   Leitura do arquivo .ini que contém as condições iniciais 
    call leitura

!   Determinação do tipo do problema e construção da solução
    call RiemannSolver

!   Escrita das solução obtida     
    call output
    
!   Executa o plot da solução 
    call execute_command_line('gnuplot -p plot_db_ERiemann.plt . &>/dev/null')
    
    call cpu_time(toc)
            
    write(*,*) '--------------------------------'           
    write(*,*) 'Tempo de execução   = ', toc-tic, 'segundos'
    write(*,*) '--------------------------------'    

end program DB_ERiemann
!--------------------------------------------------------------------
!Leitura dos dados de entrada
!--------------------------------------------------------------------
subroutine leitura
    implicit none
    
    double precision     :: g, xmin, xmax, dx, tol, dam
    double precision     :: hL, hR, vL, vR, cL, cR, tmax 
    integer :: m, nint
    
    common /cini/ hL, hR, vL, vR, cL, cR
    common /malh/ xmin, xmax, dx, dam, m
    common /tole/ tol, nint
    common /temp/ tmax
    common /grav/ g
    
    open(unit = 1, file = 'db_ERiemann.ini', status = 'unknown')
    
    read(1,*) xmin      ! Ponto inicial do canal
    read(1,*) xmax      ! Ponto final do canal 
    read(1,*) dam       ! Posição da barragem
    read(1,*) hL        ! Altura da lamina d'água à esquerda da descontinuidade   
    read(1,*) hR        ! Altura da lamina d'água à direita da descontinuidade       
    read(1,*) vL        ! Velocidade à esquerda da descontinuidade
    read(1,*) vR        ! Velocidade à direita da descontinuidade
    read(1,*) tmax      ! Instante do tempo para solução final
    read(1,*) tol       ! Tolerancia para o método de Newton
    read(1,*) nint      ! Numero maximo de iterações no método de Newton
    read(1,*) m         ! Numero de pontos/celulas no dominio
    read(1,*) g         ! Aceleração da gravidade
    
    close(1)
    
    write(*,*) '--------------------------------'
    write(*,*) 'Dados iniciais:'
    write(*,*) '--------------------------------'
    write(*,*) 'xmin  = ', xmin
    write(*,*) 'xmax  = ', xmax
    write(*,*) 'dam   = ', dam
    write(*,*) 'hL    = ', hL
    write(*,*) 'hR    = ', hR
    write(*,*) 'vL    = ', vL
    write(*,*) 'vR    = ', vR
    write(*,*) 'tmax  = ', tmax
    write(*,*) 'tol   = ', tol
    write(*,*) 'nint  = ', nint
    write(*,*) 'm     = ', m
    write(*,*) 'g     = ', g
    write(*,*) '--------------------------------'
    
end subroutine leitura
!--------------------------------------------------------------------
!Escolha do solver de Riemann
!--------------------------------------------------------------------
subroutine RiemannSolver
    implicit none

    double precision :: g, hL, hR, hs, vL, vR, vs, cL, cR, cs
       
    common /cini/ hL, hR, vL, vR, cL, cR
    common /star/ hs, vs, cs
    common /grav/ g
    
    cL = sqrt(g*hL)
    cR = sqrt(g*hR)
    cs = (vR-vL)-2*(cL+cR)
    
    if (hL.EQ.0.OR.hR.EQ.0.OR.cs.GE.0) then 
        call drybed
    else
        call wetbed
    end if
    
end subroutine RiemannSolver
!--------------------------------------------------------------------
!Solver de Riemann para leito seco (drybed)
!--------------------------------------------------------------------
subroutine drybed
    implicit none
    
    double precision :: g, xmin, xmax, dx, tmax, dam, x, h, v, s
    double precision :: hL, hR, vL, vR, cL, cR
    integer :: i, m 
    integer, parameter :: mx = 5000 
    double precision, dimension (2,mx+1) :: U
    
    common /cini/ hL, hR, vL, vR, cL, cR
    common /malh/ xmin, xmax, dx, dam, m
    common /temp/ tmax
    common /grav/ g
    common /vari/ U
    
    dx = (xmax-xmin)/dble(m)

!   Construção da solução      
    do i = 1,m+1
        x = xmin+(i-1)*dx-dam
        s = x/tmax
        
        if (hL.EQ.0) then
            call drylef(h,v,s)
        else
            if (hR.EQ.0) then
                call dryrig(h,v,s)
            else
                call drymid(h,v,s)
            end if
        end if
        
        U(1,i) = h
        U(2,i) = h*v
    end do
        
end subroutine drybed
!--------------------------------------------------------------------
!Construção da solução para o caso de leito seco à esquerda (dry left, hL = 0)
!--------------------------------------------------------------------
subroutine drylef(h,v,s)
    implicit none
    
    double precision :: g, h, v, s, sR, sRs
    double precision :: hL, hR, vL, vR, cL, cR
    
    common /cini/ hL, hR, vL, vR, cL, cR
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
    
    common /cini/ hL, hR, vL, vR, cL, cR
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
    
    common /cini/ hL, hR, vL, vR, cL, cR
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
subroutine wetbed
    implicit none

    double precision :: g, xmin, xmax, dx, tmax, dam, tol, x, h, v, s
    double precision :: f, df, hz, fL, fR, dfL, dfR
    double precision :: hL, hR, hs, vL, vR, vs, cL, cR, cs
    integer :: i, k, m, nint 
    integer, parameter :: mx = 5000 
    double precision, dimension (2,mx+1) :: U
    
    common /cini/ hL, hR, vL, vR, cL, cR
    common /malh/ xmin, xmax, dx, dam, m
    common /star/ hs, vs, cs
    common /tole/ tol, nint
    common /temp/ tmax
    common /grav/ g
    common /vari/ U
    
    dx = (xmax-xmin)/dble(m)
    
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
    
    write(*,*) '--------------------------------'           
    write(*,*) 'Numero de iterações = ', k
    write(*,*) '--------------------------------'
    
    cs = sqrt(g*hs)
    vs = 0.5*(vR+vL+fR-fL)

!   Construção da solução    
    do i = 1,m+1
        x = xmin+(i-1)*dx-dam
        s = x/tmax

        call samwet(h,v,s)
        
        U(1,i) = h
        U(2,i) = h*v
    end do
    
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
    
    common /cini/ hL, hR, vL, vR, cL, cR
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
!Geração do arquivo de saída
!--------------------------------------------------------------------
subroutine output
    implicit none
    
    double precision :: xmin, xmax, dx, dam, x
    integer :: i, m
    integer, parameter :: mx = 5000 
    double precision, dimension (2,mx+1) :: U
    
    common /malh/ xmin, xmax, dx, dam, m
    common /vari/ U
    
    dx = (xmax-xmin)/dble(m)
    
    open(unit = 1, file = 'db_ERiemann.out', status = 'unknown')
    do i = 1,m+1
        x = xmin+(i-1)*dx-dam
        
        if (U(1,i).GT.0) then        
            write(1,*) x, U(1,i), U(2,i), U(2,i)/U(1,i)
        else
            write(1,*) x, U(1,i), U(2,i), 0
        end if
    end do

    close(1)
     
end subroutine output
