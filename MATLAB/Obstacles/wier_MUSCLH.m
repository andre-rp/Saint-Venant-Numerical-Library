%   Metodo MUSCL-Hancock 
% 
%   Solucao das equacoes de Saint-Venant utilizando esquema MUSCL para 
%   reconstrucao linear e o passo de Hancock para evoluir 0.5*dt, com  
%   fluxo HLL e depois avancar a solucao no tempo.
%   
%   Os testes suportados por esse programa sao: Escoamento sobre acude,
%   teste estatico e os demais testes ja executados a menos de mudancas
%   nas condicoes inicias e de contorno. 
tic
clear all, close all, clc

%variaveis globais
global g dx dt m n cr b lim epsilon tstatic t tmax
global xmin xmax hs qc qin

%Constantes
g   = 9.8; cr      = 0.9; ntmax = 10e6; 
lim = 1;   epsilon = 1e-16; 

%malha
xmin = -1; xmax = 1; m = 200; 

%dados do problema
b   = 1; hs   = 0.3; qc = 0.01; qin = 0.111;
n   = 0.01;

t   = 0; tmax = 200; tstatic = 1; I = linspace(2,m-1,m-2);

[x,zb,zL,zR,U] = start;

for k = 1:ntmax 
    
    dt = CFLcon(U,k);
    
%   Reconstrucao da informacao dentro da celula.
    BEXT = MUSCL(U,zL,zR);
    
%   Formacao dos problemas de Riemann utilizando a extrapolacao para os
%   valores nas fronteiras de cada celula.
    F = HLLRS(BEXT);

%   Atualizacao do problema homogeneo.
    U = update(U,F);
    
%   Acrescimo dos termos fontes na solucao homogenea.
    U = source(U,BEXT,zL,zR);
    
    if ~tstatic
        U = boundary(U);
    end

    if t >= tmax
        for i = 1:m+2
            if U(1,i) <= epsilon
                v(i) = 0;
            else
                v(i) = U(2,i)/U(1,i);
            end
        end
        output(x(I),zb(I),U(:,I),v(I));
        break
    end
    k = k+1;
end

figure(1)
plot(x(I),U(1,I)+zb(I),'b',x(I),zb(I),'k')
xlim([xmin xmax]), ylim([0 max(U(1,:))+0.1])
xlabel('x (m)'), ylabel('h+zb (m)')
title([t k]), grid on

figure(2)
plot(x(I),v(I),'b')
xlim([xmin xmax]), ylim([min(v(I))-0.1 max(v(I))+0.1])
xlabel('x (m)'), ylabel('v (m/s)')
title([t k])

figure(3)
plot(x(I),U(2,I),'b')
xlim([xmin xmax]), ylim([min(U(2,:))-0.1 max(U(2,:))+0.1])
xlabel('x (m)'), ylabel('q (m^2/s)')
title([t k]), grid on 

toc
%-----------------------------------------------------------------------
%Inicializacao do vetor de variaveis x, zb e  U
%-----------------------------------------------------------------------
function [x,zb,zL,zR,U] = start
    global dx m tstatic
    global xmin xmax
    
    U = zeros(2,m+2); [x,zb,zL,zR] = deal(zeros(1,m+2)); 
     
    dx = (xmax-xmin)/m;
    
    x(1) = xmin;
    for i = 2:m+1
        x(i) = x(i-1)+dx;
        
        zL(i) = 0.2*exp(-8.68055*x(i-1)^2); 
        zR(i) = 0.2*exp(-8.68055*x(i)^2);
        
        xc    = x(i)-0.5*dx;
        zb(i) = 0.2*exp(-8.68055*xc^2); ;
    end
    
    if ~tstatic
        U = wstart(U);
    else
        U = sstart(U,zb);
    end
end
%-----------------------------------------------------------------------
%Inicializa a variavel U para o escoamento sobre o acude
%-----------------------------------------------------------------------
function U = wstart(U)
    global g dx m qc qin
    
    hc = (qc^2/g)^(1/3);
    
    hp = hc; hn = hc;
    xp = 0;  xn = 0;   
    aux = round(m/2);
    
    U(1,aux) = hc;
    U(2,aux) = qc;
    
    for i = aux+1:m
        hp = frk4(hp,xp,hc,dx);
        xp = xp+dx;
        
        U(1,i) = hp;
        U(2,i) = qc;
    end

    for i = aux-1:-1:1
        hn = frk4(hn,xn,hc,-dx);
        xn = xn-dx;
        
        U(1,i) = hn;
        U(2,i) = qc;
    end
end
%-----------------------------------------------------------------------
%Aplica o metodo de Runge-Kutta de 4 ordem para determinar h
%-----------------------------------------------------------------------
function h = frk4(hz,xz,hc,deltax)
     
    h = hz;
    x = xz;
    k1 = deriva(h,x,hc,deltax);
    
    h = hz+0.5*k1;
    x = xz+0.5*deltax;
    k2 = deriva(h,x,hc,deltax);
    
    h = hz+0.5*k2;
    x = xz+0.5*deltax;
    k3 = deriva(h,x,hc,deltax);
    
    h = hz+k3;
    x = xz+deltax;
    salto4 = deriva(h,x,hc,deltax);
    
    h = hz+(k1+2*k2+2*k3+salto4)/6;
    
end
%-----------------------------------------------------------------------
%Calcula a derivada para o metodo de Runge-Kutta
%-----------------------------------------------------------------------
function kn = deriva(h,x,hc,deltax)
    global g qc 

    zb  = 0.2*exp(-8.68055*x^2);
    dz  = 0.2*exp(-8.68055*x^2)*(-8.68055*2*x);
    ddz = dz*(-8.68055*2*x)+zb*(-8.68055*2);
    
    F = sqrt(qc^2/(g*h^3));
    
    if x == 0
        dhdx = -sqrt(-ddz*hc/3);
    else
        dhdx = -dz/(1-F^2);
    end
    kn = deltax*dhdx;
end
%-----------------------------------------------------------------------
%Calcula a derivada para o metodo de Runge-Kutta
%-----------------------------------------------------------------------
function U = sstart(U,zb)
    global m hs
    
    for i = 1:m+2
        U(1,i) = hs-zb(i);
        U(2,i) = 0;
    end
end
%-----------------------------------------------------------------------
%Escolha do passo de tempo
%-----------------------------------------------------------------------
function dt = CFLcon(U,k)
    global g dx dt m cr epsilon t tmax
    
    smax = 0; speloc = 0;
    
    for i = 2:m+1
        if U(1,i) <= epsilon
            speloc = 0;
        else
            speloc = abs(U(2,i)/U(1,i))+sqrt(g*U(1,i));
        end
        
        if speloc > smax
            smax = speloc;
        end         
    end
    
    dt = cr*dx/smax;
    
    if k <= 5
        dt = 0.2*dt;
    end
    
    if t+dt > tmax
        dt = tmax-t;
    end
    t = t+dt;
end
%-----------------------------------------------------------------------
%MUSCL
%-----------------------------------------------------------------------
function BEXT = MUSCL(U,zL,zR)
    global g dx dt m lim epsilon
    
    BEXT = zeros(2,2,m+2); [uL,uR] = deal(zeros(2,1));
    
    for i = 2:m+1
        for j = 1:2
            dup = U(j,i)-U(j,i-1);
            ddw = U(j,i+1)-U(j,i);

%           Limitador de inclinacao            
            if dup*ddw <= 0
                delta = 0;
            else
                mino = min([abs(ddw) abs(dup)]);
                delta = sign(ddw)*min([mino lim*mino]);
            end
                
%           Variaves (h,q) extrapoladas para as fronteiras a esquerda e a
%           direita.
            uL(j) = U(j,i)-0.5*delta; uR(j) = U(j,i)+0.5*delta;
        end
        
        fL = flux(uL);
        fR = flux(uR);
        Sc = msor(U(:,i),uL,uR,zR(i),zL(i));
        
%       Passo de Hancock para evoluir a solucao dt/2.
        for j = 1:2
            BEXT(1,j,i) = uL(j)+0.5*dt/dx*(fL(j)-fR(j))+0.5*dt*Sc(j);
            BEXT(2,j,i) = uR(j)+0.5*dt/dx*(fL(j)-fR(j))+0.5*dt*Sc(j);
        end
        
%       Conservacao da positividade na profundidade.
        if BEXT(1,1,i) < 0
            BEXT(1,1,i) = 0;
            BEXT(1,2,i) = 0;
        end
        if BEXT(2,1,i) < 0
            BEXT(2,1,i) = 0;
            BEXT(2,2,i) = 0;
        end
        
%       Evitar altas velocidade.
        if BEXT(1,1,i) < epsilon
            BEXT(1,2,i) = 0;
        end
        if BEXT(2,1,i) < epsilon
            BEXT(2,2,i) = 0;
        end
    end   
end
%-----------------------------------------------------------------------
%Fluxo numerico utilizado no MUSCL
%-----------------------------------------------------------------------
function fK = flux(uK)
    global g epsilon
    
    hK = uK(1); 
    if hK <= epsilon
        fK = [0; 0];
    else
        vK = uK(2)/hK;
        fK = [hK*vK; hK*vK^2+0.5*g*hK^2];
    end
end
%-----------------------------------------------------------------------
%Termos fontes para o passo de Hancock
%-----------------------------------------------------------------------
function Sc = msor(Uk,uL,uR,zR,zL)
    global g dx n b epsilon

    hk = Uk(1);
    if hk <= epsilon
        Sc = [0; 0];
    else
        vk = Uk(2)/hk;
        
        Rh = b*hk/(2*hk+b);
        sf = (abs(vk)*vk*n^2)/Rh^(4/3);
        
        Sc(1) = 0;
        Sc(2) = -g*hk*sf...
            -g*(uL(1)+uR(1))*0.5*(zR-zL)/dx;
    end
end
%-----------------------------------------------------------------------
%HLLRS
%-----------------------------------------------------------------------
function F = HLLRS(BEXT)
    global g m epsilon

    [CDL,CDR] = deal(zeros(2,1)); F = zeros(2,m+2);
    
    for i = 2:m
        for j = 1:2
            CDL(j) = BEXT(2,j,i);
            CDR(j) = BEXT(1,j,i+1);
        end

        iswet = 1;
        
        hL = CDL(1); cL = sqrt(g*hL);
        if hL <= epsilon
            vL = 0;
        else
            vL = CDL(2)/hL; 
        end
            
        hR = CDR(1); cR = sqrt(g*hR);
        if hR <= epsilon
            vR = 0;
        else
            vR = CDR(2)/hR; 
        end
        
        fL = flux(CDL); 
        fR = flux(CDR);
        
%       Velocidade da crista e da cauda para as ondas no caso de leito 
%       seco
        if hR < epsilon & hL > epsilon
%           Leito seco a direita 
            sL = vL-cL;
            sR = vL+2*cL;
            iswet = 0;
        end
        if hL < epsilon & hR > epsilon
%           Leito seco a esquerda
            sL = vR-2*cR;
            sR = vR+cR;
            iswet = 0;
        end
        if hL < epsilon & hR < epsilon
%           Leito seco
            sL = -epsilon;
            sR = epsilon;
            iswet = 0;
        end

%       Velocidade da crista e da cauda para as ondas no caso de leito 
%       molhado
        if iswet
            hs = 1/g*(0.5*(cL+cR)+0.25*(vL-vR))^2;
            
            if hs <= hL
                sL = vL-cL;
            else
                sL = vL-cL*sqrt(0.5*hs*(hs+hL))/hL;
            end
            
            if hs <= hR
                sR = vR+cR;
            else
                sR = vR+cR*sqrt(0.5*hs*(hs+hR))/hR;
            end
        end
        
%       Calculo fluxo na interface (i,i+1) em x/t=0
        if sL >= 0
            F(:,i) = fL;
        end
        
        if sL <= 0 & sR >= 0
            FHLL = sR*fL-sL*fR+sL*sR*(CDR-CDL);
            F(:,i) = FHLL/(sR-sL);
        end
        
        if sR <= 0
            F(:,i) = fR;
        end         
    end
end
%-----------------------------------------------------------------------
%Atualizacao da solucao
%-----------------------------------------------------------------------
function U = update(U,F)
    global dx dt m epsilon
    
    for i = 3:m
        U(:,i) = U(:,i)-dt/dx*(F(:,i)-F(:,i-1));
        
%       Conservacao da positividade na profundidade.
        if U(1,i) < 0
            U(:,i) = 0;
        end
%       Evitar altas velocidade.
        if U(1,i) <= epsilon
            U(2,i) = 0;
        end
    end
end
%-----------------------------------------------------------------------
%Termos fontes - Equacao de Manning + Derivada de zb
%-----------------------------------------------------------------------
function U = source(U,BEXT,zL,zR)
    global g dx dt m n b epsilon
    
    for i = 3:m
        hk = U(1,i);
        if hk <= epsilon
            S(1,i) = 0;
            S(2,i) = 0;
        else
            vk = U(2,i)/hk;
            
            Rh = b*hk/(2*hk+b);
            sf = (abs(vk)*vk*n^2)/Rh^(4/3);
            
            S(1,i) = 0;
            S(2,i) = -g*hk*sf...
                -g*(BEXT(1,1,i)+BEXT(2,1,i))*0.5*(zR(i)-zL(i))/dx;
        end
        U(2,i) = U(2,i)+dt*S(2,i);
    end
end
%-----------------------------------------------------------------------
%Condicoes de fronteira
%-----------------------------------------------------------------------
function U = boundary(U)
    global m qin
    
    U(1,2) = U(1,3);
    U(2,2) = qin;
    U(:,1) = U(:,2);
    
    U(:,m+1) = U(:,m);
    U(:,m+2) = U(:,m+1);
end
%-----------------------------------------------------------------------
%Geracao do arquivo de saida
%-----------------------------------------------------------------------
function output(x,zb,U,v)
    outfile = [x;U;v;zb];
    fileID = fopen('static_MUSCLH.out','w');
    fprintf(fileID,'%6.2f %12.8f %12.8f %12.8f %12.8f\r\n',outfile);
    fclose(fileID);
end