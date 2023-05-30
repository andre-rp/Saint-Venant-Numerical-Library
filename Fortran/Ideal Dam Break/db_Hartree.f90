!   Método de Hartree ou Método dos intervalos especificos 
! 
!   Solução das equações de Saint-Venant utilizando as curvas caracteristicas
!   para avançar a solução no tempo.
!   
!   Os testes suportados por esse programa são: Escoamentos subcrítico e 
!   transcrítico.

program DB_Hartree
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
        
!       Aplicação do método de Hartree para interpolar e avançar a solução
        call Hartree
               
!       Aplicação das condições de fronteira (fronteira simetrica)
        call boundary
      
        if (t.GE.tmax) then
!           Escrita das informações ao final do tempo maximo de execução        
            call output
            
!           Executa o plot da solução no tempo final
            call execute_command_line('gnuplot -p plot_db_Hartree.plt')
        
            call cpu_time(toc)
            
            write(*,*) '--------------------------------'           
            write(*,*) 'Numero de iterações = ', k
            write(*,*) 'Tempo de execução   = ', toc-tic, 'segundos'
            write(*,*) '--------------------------------'
            exit
        end if        
    end do

end program DB_Hartree
!--------------------------------------------------------------------
!Leitura dos dados de entrada
!--------------------------------------------------------------------
subroutine leitura
    implicit none
    
    double precision     :: cr, g, dx, dt, xmin, xmax
    double precision     :: b, n, hu, hd, vu, vd, t, tmax, dam 
    integer :: m, ntmax
    
    common /domi/ xmin, xmax, dam
    common /cini/ hu, hd, vu, vd    
    common /temp/ t, tmax, ntmax
    common /malh/ dx, dt, m
    common /cons/ n, cr, b
    common /grav/ g
    
    open(unit = 1, file = 'db_Hartree.ini', status = 'unknown')
    
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
            U(1,i) = vu
            U(2,i) = sqrt(g*hu)
        else
            U(1,i) = vd
            U(2,i) = sqrt(g*hd)
        end if
    end do

end subroutine start
!--------------------------------------------------------------------
!Atualização do passo de tempo
!--------------------------------------------------------------------
subroutine CFLcon(k)
    implicit none

    double precision     :: cr, dx, dt, n, b, t, tmax 
    integer :: m, k, ntmax
    integer, parameter :: mx = 5000
    double precision, dimension (2,-1:mx+1) :: U
    double precision, dimension (-1:mx+1) :: x, zb
    
    common /temp/ t, tmax, ntmax
    common /malh/ dx, dt, m
    common /cons/ n, cr, b
    common /vari/ U, x, zb

    dt = cr*dx*minval(1/(abs(U(1,1:m))+U(2,1:m)))

    if (k.LE.5) then
        dt = 0.2*dt
    end if

    if ((t+dt).GT.tmax) then
        dt = tmax-t
    end if
    
    t = t+dt

end subroutine CFLcon
!--------------------------------------------------------------------
!Método de Hartree
!Interpolação dos valores intermediarios e avanço da solução
!--------------------------------------------------------------------
subroutine Hartree
    implicit none
    
    double precision     :: g, dx, dt
    double precision     :: va, vb, vc, ca, cb, cc 
    integer :: i, m
    integer, parameter :: mx = 5000
    double precision, dimension (mx) :: hL, hR, vL, vR, cL, cR, sL, sR, dxL, dxR
    double precision, dimension (2,-1:mx+1) :: U
    double precision, dimension (-1:mx+1) :: x, zb
    
    common /malh/ dx, dt, m
    common /vari/ U, x, zb
    common /grav/ g
    
    do i = 1,m
        dxL(i) = dt*(U(1,i)+U(2,i))
        dxR(i) = dt*(U(1,i)-U(2,i))
    end do

    do i = 2,m-1
        va = U(1,i-1)
        ca = U(2,i-1) 
        vb = U(1,i) 
        cb = U(2,i)
        
        vL(i) = (vb*dx+dt*(cb*va-ca*vb))/(dx+dt*(vb+cb-va-ca))
        cL(i) = (cb*dx+vL(i)*dt*(ca-cb))/(dx+dt*(cb-ca))
        hL(i) = (cL(i)**2)/g
        
        call source(hL(i),vL(i),sL(i))
    end do
    
    do i = 2,m-1
        if (dxR(i).LT.0) then
            vb = U(1,i)
            cb = U(2,i)
            vc = U(1,i+1)
            cc = U(2,i+1)

            vR(i) = (vb*dx+dt*(cb*vc-cc*vb))/(dx+dt*(vc-cc-vb+cb))
            cR(i) = (cb*dx+vR(i)*dt*(cb-cc))/(dx+dt*(cb-cc))
            hR(i) = (cR(i)**2)/g
        else
            va = U(1,i-1)
            ca = U(2,i-1)
            vb = U(1,i)
            cb = U(2,i)
                
            vR(i) = (vb*dx+dt*(ca*vb-cb*va))/(dx+dt*(vb-cb-va+ca))
            cR(i) = (cb*dx+vR(i)*dt*(ca-cb))/(dx+dt*(ca-cb))
            hR(i) = (cR(i)**2)/g;
        end if
        
        call source(hR(i),vR(i),sR(i))
    end do

    do i = 2,m-1
        U(1,i) = 0.5*(2*(cL(i)-cR(i))+vL(i)+vR(i)-g*dt*(sL(i)+sR(i)))
        U(2,i) = 0.25*(2*(cL(i)+cR(i))+vL(i)-vR(i)-g*dt*(sL(i)-sR(i)))
    end do

end subroutine Hartree
!--------------------------------------------------------------------
!Termo fonte - Equação de Manning + Derivada de zb
!--------------------------------------------------------------------
subroutine source(hk,vk,sK)
    implicit none
    
    double precision  :: cr, n, b
    double precision  :: hk, vk, sk, Rh

    common /cons/ n, cr, b

    Rh = (b*hk)/(b+2*hk)
    
    sk = (abs(vk)*vk*n**2)/(Rh**(4./3.))
    
end subroutine
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
    
    U(1,1) = 0
    U(2,1) = U(2,2)

    U(1,m) = 0
    U(2,m) = U(2,m-1)


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
    
    open(unit = 1, file = 'db_Hartree.out', status = 'unknown')
    do i = 1,m
        if (U(1,i).GT.0) then        
            write(1,*) x(i), (U(2,i)**2)/g, U(1,i)*(U(2,i)**2)/g, U(1,i)
        else
            write(1,*) x(i), (U(2,i)**2)/g, U(1,i)*(U(2,i)**2)/g, 0
        end if
    end do

    close(1)
     
end subroutine output
