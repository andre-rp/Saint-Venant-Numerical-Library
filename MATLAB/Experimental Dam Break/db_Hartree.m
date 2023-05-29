%   Metodo de Hartree ou Metodo dos Intervalos Especificos
%
%   Solucao das equacoes de Saint-Venant utilizando as curvas caracteristicas
%   para avancar a solucao no tempo.
%
%   Os testes suportados por esse programa sao: Escoamentos subcritico e
%   transcritico.
tic
clear all, close all, clc

%variaveis globais
global g dx dt m n cr b t tmax
global xmin xmax hu hd vu vd dam

%constantes
g = 9.8; cr = 0.9; ntmax = 10e6;

%malha
xmin = -4.65; xmax = 4.35; m = 901; 

%dados do problema
b   = 0.3; hu   = 0.25; hd = 0.10; vu = 0; vd = 0;
dam = 0;   n    = 0.005;

T = 8.9; t = 0; tmax = T/sqrt(g/hu);

[x,zb,U] = start;

for k = 1:ntmax 
    
    dt = CFLcon(U,k);

    U = hartree(U);
    
    U = boundary(U);
    
    if t >= tmax
        output(x/hu,U(1,:),[(U(2,:).^2/g)/hu; U(1,:).*U(2,:).^2/g]);
        break
    end
    k = k+1;
end

figure(1)
plot(x/hu,(zb+U(2,:).^2/g)/hu, 'b')
xlim([xmin/hu xmax/hu]), ylim([0 max((zb+U(2,:).^2/g)/hu)+0.1])
xlabel('x/h_u'), ylabel('h/h_u')
title([t k]), grid on 

figure(2)
plot(x, U(1,:), 'b')
xlim([xmin xmax]), ylim([min(U(1,:))-0.5 max(U(1,:))+0.5])
xlabel('x (m)'), ylabel('v (m/s)')
title([t k]), grid on 

figure(3)
plot(x, b*U(1,:).*U(2,:).^2/g, 'b')
xlim([xmin xmax]), ylim([min(U(1,:))-0.1 max(b*U(1,:).*U(2,:).^2/g)+0.1])
xlabel('x (m)'), ylabel('q (m^2/s)')
title([t k]), grid on 

toc

%-----------------------------------------------------------------------
%Inicialização do vetor de variaveis x, zb e  U
%-----------------------------------------------------------------------
function [x,zb,U] = start
    global xmin xmax hu hd vu vd dam
    global g dx m

    U = zeros(2,m); [x,zb] = deal(zeros(1,m));    
    
    dx = (xmax-xmin)/(m-1);
    
    for i = 1:m
        x(i)  = xmin+(i-1)*dx;
        zb(i) = 0;
        
        if x(i) <= dam
            U(1,i) = vu;
            U(2,i) = sqrt(g*hu);
        else
            U(1,i) = vd;
            U(2,i) = sqrt(g*hd);
        end
    end
end
%-----------------------------------------------------------------------
%Escolha do passo de tempo
%-----------------------------------------------------------------------
function dt = CFLcon(U,k)
    global g dx dt m n cr b t tmax

    dt = cr*min(dx./(abs(U(1,:))+U(2,:)));
   
    if k <= 5
        dt = 0.2*dt;
    end
    
    if t+dt > tmax
        dt = tmax-t;
    end
    t = t+dt;
end
%-----------------------------------------------------------------------
%Metodo de Hartree
%Interpolacao dos valores intermediarios e avanco da solucao
%-----------------------------------------------------------------------
function U = hartree(U)
    global g dx dt m 
    
    [vL,cL,hL,sL] = deal(zeros(1,m)); [vR,cR,hR,sR] = deal(zeros(1,m)); 

    dxL = dt*(U(1,:)+U(2,:)); dxR = dt*(U(1,:)-U(2,:));

%   Caracteristica Positiva 
    for i = 2:m-1
        va = U(1,i-1); ca = U(2,i-1); vb = U(1,i); cb = U(2,i);
        
        vL(i) = (vb*dx+dt*(cb*va-ca*vb))/(dx+dt*(vb+cb-va-ca));
        cL(i) = (cb*dx+vL(i)*dt*(ca-cb))/(dx+dt*(cb-ca));
        hL(i) = cL(i)^2/g;
        
        sL(i) = source(hL(i),vL(i)); 
    end

%   Caracteristica Negativa    
    for i = 2:m-1
        if dxR(i) < 0
            vb = U(1,i); cb = U(2,i); vc = U(1,i+1); cc = U(2,i+1);

            vR(i) = (vb*dx+dt*(cb*vc-cc*vb))/(dx+dt*(vc-cc-vb+cb));
            cR(i) = (cb*dx+vR(i)*dt*(cb-cc))/(dx+dt*(cb-cc));
            hR(i) = cR(i)^2/g;
        else
            va = U(1,i-1); ca = U(2,i-1); vb = U(1,i); cb = U(2,i);
            
            vR(i) = (vb*dx+dt*(ca*vb-cb*va))/(dx+dt*(vb-cb-va+ca));
            cR(i) = (cb*dx+vR(i)*dt*(ca-cb))/(dx+dt*(ca-cb));
            hR(i) = cR(i)^2/g;
        end
        sR(i) = source(hR(i),vR(i));
    end

%   Atualizacao da solucao
    for i = 2:m-1
        U(1,i) = 1/2*(2*(cL(i)-cR(i))+vL(i)+vR(i)-g*dt*(sL(i)+sR(i)));
        U(2,i) = 1/4*(2*(cL(i)+cR(i))+vL(i)-vR(i)-g*dt*(sL(i)-sR(i)));
    end
end
%-----------------------------------------------------------------------
%Termos fontes - Equacao de Manning + Derivada de zb
%-----------------------------------------------------------------------
function Sk = source(hk,vk)
    global n cr b 
    
    Rh = (b*hk)/(b+2*hk);
    
    Sk = (abs(vk)*vk*n^2)/Rh^(4/3);
end
%-----------------------------------------------------------------------
%Condicoes de fronteira
%-----------------------------------------------------------------------
function U = boundary(U)
    global m
    
    U(1,1) = 0;
    U(2,1) = U(2,2);
    
    U(1,m) = 0;
    U(2,m) = U(2,m-1);
end
%-----------------------------------------------------------------------
%Geracao do arquivo de saida
%-----------------------------------------------------------------------
function output(x,v,U)
    outfile = [x;U;v];
    fileID = fopen('db_Hartree.out','w');
    fprintf(fileID,'%6.2f %12.8f %12.8f %12.8f\r\n',outfile);
    fclose(fileID);
end