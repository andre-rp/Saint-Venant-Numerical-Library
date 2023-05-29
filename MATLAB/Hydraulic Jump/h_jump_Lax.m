%   Metodo de Lax difusivo 
% 
%   Solucao das equacoes de Saint-Venant utilizando o metodo de Lax difusivo 
%   para avancar a solucao no tempo.
%   
%   Os testes suportados por esse programa sao: Salto hidraulico e os
%   demais testes ja executados a menos de mudancas nas condicoes de
%   alteracoes nas condicoes de contorno.
tic
clear all, close all, clc

%variaveis globais
global g dx dt m n cr b alpha t tmax
global xmin xmax hu vu 

%constantes
g = 9.8; cr = 0.9; alpha = 0.9; ntmax = 10e6;

%malha
xmin = 0.305; xmax = 14; m = 151; 

%dados do problema
b   = 0.46; hu   = 0.031; Fr = 7; vu = Fr*sqrt(g*hu); 
n   = 0.008;

t   = 0; tmax = 400;

[x,zb,U] = start;

for k = 1:ntmax 
    
    dt = CFLcon(U,k);

    U = Lax(U);
    
    U = boundary(U);

    if t >= tmax
        output(x,U);
        break
    end
    k = k+1;
end

figure(1)
plot(x,U(1,:)+zb,'b')
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
%Inicialização do vetor de variaveis x, zb e  U
%-----------------------------------------------------------------------
function [x,zb,U] = start
    global dx m
    global xmin xmax hu vu
    
    U = zeros(2,m); [x,zb] = deal(zeros(1,m));    
    
    dx = (xmax-xmin)/(m-1);
    for i = 1:m
        x(i) = xmin+(i-1)*dx;
        zb(i) = 0;
        
        U(1,i) = hu;
        U(2,i) = hu*vu;
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
function U = Lax(U)
    global dx dt m alpha
    
    Uz = U;
    for i = 2:m-1
        U(:,i) = alpha*Uz(:,i)+(1-alpha)*0.5*(Uz(:,i+1)+Uz(:,i-1))...
            -0.5*dt/dx*(flux(Uz(:,i+1))-flux(Uz(:,i-1)))+dt*source(Uz(:,i));
    end  
end
%-----------------------------------------------------------------------
%Avaliação do vetor de fluxo no ponto U_i no tempo t 
%-----------------------------------------------------------------------
function F = flux(U)
	global g
    
    hk = U(1);
    if hk <= 0
        F = [0; 0]
    else
        vk = U(2)/hk;
        F = [hk*vk; hk*vk^2+0.5*g*hk^2];
    end
end
%-----------------------------------------------------------------------
%Termos fontes - Equacao de Manning + Derivada de zb
%-----------------------------------------------------------------------
function S = source(U)
    global g n b
    
    hk = U(1);
    if hk <= 0
        S = [0; 0];
    else
        vk = U(2)/hk;
    
        Rh = (b*hk)/(b+2*hk);
    
        S = [0; -g*hk*(abs(vk)*vk*n^2)/Rh^(4/3)];
    end
end
%-----------------------------------------------------------------------
%Condicoes de fronteira
%-----------------------------------------------------------------------
function U = boundary(U)
    global m t hu vu
    
    U(1,1) = hu;
    U(2,1) = hu*vu;
    
    U(1,m) = min([0.265 0.031+0.00468*t]);
    U(2,m) = U(2,m-1);
end
%-----------------------------------------------------------------------
%Geracao do arquivo de saida
%-----------------------------------------------------------------------
function output(x,U)
    v = U(2,:)./U(1,:);
    outfile = [x;U;v];
    fileID = fopen('h_jump_Lax.out','w');
    fprintf(fileID,'%6.2f %12.8f %12.8f %12.8f\r\n',outfile);
    fclose(fileID);
end