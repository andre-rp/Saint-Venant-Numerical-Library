%   Metodo de Lax difusivo 
% 
%   Solucao das equacoes de Saint-Venant utilizando o metodo de Lax difusivo 
%   para avancar a solucao no tempo.
%   
%   Os testes suportados por esse programa sao: Escoamentos subcritico e 
%   transcritico.
tic
clear all, close all, clc

%variaveis globais
global g dx dt m n cr b alpha t tmax
global xmin xmax hu hd vu vd dam

%constantes
g = 9.8; cr = 0.3; alpha = 0.9; ntmax = 10e6;

%malha
xmin = -4.65; xmax = 4.35; m = 901; 

%dados do problema
b   = 0.3; hu   = 0.25; hd = 0.10; vu = 0; vd = 0;
dam = 0;   n    = 0.005;

T = 8.9; t = 0; tmax = T/sqrt(g/hu);

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
plot(x/hu,(U(1,:)+zb)/hu,'b')
xlim([xmin/hu xmax/hu]), ylim([0 max(U(1,:)/hu)+0.1])
xlabel('x/h_u'), ylabel('h/h_u')
title([t k]), grid on

figure(2)
plot(x,U(2,:)./U(1,:),'b')
xlim([xmin xmax]), ylim([min(U(2,:)./U(1,:))-0.5 max(U(2,:)./U(1,:))+0.5])
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
    global xmin xmax hu hd vu vd dam
    
    U = zeros(2,m); [x,zb] = deal(zeros(1,m));    
    
    dx = (xmax-xmin)/(m-1);
    for i = 1:m
        x(i) = xmin+(i-1)*dx;
        zb(i) = 0;
        
        if x(i) <= dam
            U(1,i) = hu;
            U(2,i) = hu*vu;
        else
            U(1,i) = hd;
            U(2,i) = hd*vd;
        end
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
    vk = U(2)/hk;
    F = [hk*vk; hk*vk^2+0.5*g*hk^2];
end
%-----------------------------------------------------------------------
%Termos fontes - Equacao de Manning + Derivada de zb
%-----------------------------------------------------------------------
function S = source(U)
    global g n b
    
    hk = U(1);
    vk = U(2)/hk;
    
    Rh = (b*hk)/(b+2*hk);
    
    S = [0; -g*hk*(abs(vk)*vk*n^2)/Rh^(4/3)];
end
%-----------------------------------------------------------------------
%Condicoes de fronteira
%-----------------------------------------------------------------------
function U = boundary(U)
    global m
    
    U(1,1) = U(1,2);
    U(2,1) = 0;
    
    U(1,m) = U(1,m-1);
    U(2,m) = 0;
end
%-----------------------------------------------------------------------
%Geracao do arquivo de saida
%-----------------------------------------------------------------------
function output(x,U)
    global hu
    
    v = U(2,:)./U(1,:);
    outfile = [x/hu;U(1,:)/hu;U(2,:);v];
    fileID = fopen('db_Lax.out','w');
    fprintf(fileID,'%6.2f %12.8f %12.8f %12.8f\r\n',outfile);
    fclose(fileID);
end