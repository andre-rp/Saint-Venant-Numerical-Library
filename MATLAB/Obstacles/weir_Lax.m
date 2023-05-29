%   Metodo de Lax difusivo 
% 
%   Solucao das equacoes de Saint-Venant utilizando o metodo de Lax 
%   difusivo para avancar a solucao no tempo.
%     
%   Os testes suportados por esse programa sao: Escoamento sobre acude,
%   teste estatico e os demais testes ja executados a menos de mudancas
%   nas condicoes inicias e de contorno. 
tic
clear all, close all, clc

%variaveis globais
global g dx dt m n cr b alpha tstatic t tmax 
global xmin xmax hs qc qin

%constantes
g = 9.8; cr = 0.9; alpha = 0.9; ntmax = 10e6;

%malha
xmin = -1; xmax = 1; m = 201; 

%dados do problema
b   = 1; hs   = 0.3; qc = 0.01; qin = 0.111;
n   = 0.01;

t   = 0; tmax = 200; tstatic = 1;

[x,zb,U] = start;

for k = 1:ntmax 
    
    dt = CFLcon(U,k);

    U = Lax(x,U);
    
    if ~tstatic
        U = boundary(U);
    end

    if t >= tmax
        output(x,zb,U);
        break
    end
    k = k+1;
end

figure(1)
plot(x,U(1,:)+zb,'b',x,zb,'k')
xlim([xmin xmax]), ylim([0 max(U(1,:))+0.1])
xlabel('x (m)'), ylabel('h+zb (m)')
title([t k]), grid on

figure(2)
plot(x,U(2,:)./U(1,:),'b')
xlim([xmin xmax]), ylim([min(U(2,:)./U(1,:))-0.1 max(U(2,:)./U(1,:))+0.1])
xlabel('x (m)'), ylabel('v (m/s)')
title([t k]), grid on 

figure(3)
plot(x,U(2,:),'b')
xlim([xmin xmax]), ylim([min(U(2,:))-0.1 max(U(2,:))+0.1])
xlabel('x (m)'), ylabel('q (m^2/s)')
title([t k]), grid on 

toc

%-----------------------------------------------------------------------
%Inicializacao do vetor de variaveis x, zb e  U
%-----------------------------------------------------------------------
function [x,zb,U] = start
    global dx m tstatic
    global xmin xmax
    
    U = zeros(2,m); [x,zb] = deal(zeros(1,m));    
    
    dx = (xmax-xmin)/(m-1);
    for i = 1:m
        x(i) = xmin+(i-1)*dx;
        zb(i) = 0.2*exp(-8.68055*x(i)^2);
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
    global g dx m qc 
    
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
    
    for i = 1:m
        U(1,i) = hs-zb(i);
        U(2,i) = 0;
    end
end
%-----------------------------------------------------------------------
%Escolha do passo de tempo
%-----------------------------------------------------------------------
function dt = CFLcon(U,k)
    global g dx dt cr alpha t tmax
    
    dt = (1-alpha)*cr*min(dx./(abs(U(2,:)./U(1,:))+sqrt(g*U(1,:))));
    
    if k <= 5
        dt = 0.2*dt;
    end
    
    if t+dt > tmax
        dt = tmax-t;
    end
    t = t+dt;
end
%-----------------------------------------------------------------------
%Metodo de Lax difusivo
%-----------------------------------------------------------------------
function U = Lax(x,U)
    global dx dt m alpha
    
    Uz = U;
    for i = 2:m-1
        U(:,i) = alpha*Uz(:,i)+(1-alpha)*0.5*(Uz(:,i+1)+Uz(:,i-1))...
            -0.5*dt/dx*(flux(Uz(:,i+1))-flux(Uz(:,i-1)))+dt*source(x(i),Uz(:,i));
    end  
end
%-----------------------------------------------------------------------
%Avaliacao do vetor de fluxo no ponto U_i no tempo t 
%-----------------------------------------------------------------------
function F = flux(U)
	global g
    
    hk = U(1);
    vk = U(2)/hk;
    F = [hk*vk; hk*vk^2+0.5*g*hk^2];
end
%-----------------------------------------------------------------------
%Termos fontes - Equacao de Manning + Derivada de zb
%-----------------------------------------------------------------------
function S = source(x,U)
    global g n b
    
    hk = U(1);
    vk = U(2)/hk;
    
    Rh = (b*hk)/(b+2*hk);
    dz = 0.2*exp(-8.68055*x^2)*(-8.68055*2*x);
    
    S = [0; -g*hk*(abs(vk)*vk*n^2)/Rh^(4/3)-g*hk*dz];
end
%-----------------------------------------------------------------------
%Condicoes de fronteira
%-----------------------------------------------------------------------
function U = boundary(U)
    global m qin
    
    U(1,1) = U(1,2);
    U(2,1) = qin;
    
    U(:,m) = U(:,m-1);
end
%-----------------------------------------------------------------------
%Geracao do arquivo de saida
%-----------------------------------------------------------------------
function output(x,zb,U)
    v = U(2,:)./U(1,:);
    outfile = [x;U;v;zb];
    fileID = fopen('static_Lax.out','w');
    fprintf(fileID,'%6.2f %12.8f %12.8f %12.8f %12.8f\r\n',outfile);
    fclose(fileID);
end