%   Metodo Godunov primeira ordem 
% 
%   Solucao das equacoes de Saint-Venant utilizando esquema de Godunov 
%   de primeira ordem, que utiliza o solver de Riemann para determinar  
%   o fluxo numerico e depois avancar a solucao no tempo.
%   
%   Os testes suportados por esse programa sao: Salto hidraulico e os
%   demais testes ja executados a menos de mudancas nas condicoes de
%   alteracoes nas condicoes de contorno.
tic
clear all, close all, clc

%variaveis globais
global g dx dt m n cr b nint tol epsilon t tmax
global xmin xmax hu vu

%Constantes
g   = 9.8;  cr   = 0.9; ntmax   = 10e6; 
tol = 1e-6; nint = 50;  epsilon = 1e-16;

%malha
xmin = 0.305; xmax = 14; m = 150; 

%dados do problema
b   = 0.46; hu   = 0.031; Fr = 7; vu = Fr*sqrt(g*hu); 
n   = 0.008;

t   = 0; tmax = 400; I = linspace(2,m-1,m-2);

[x,zb,zL,zR,U] = start;

for k = 1:ntmax 

    dt = CFLcon(U,k);
    
    F = GodunovRs(U);

%   Atualizacao do problema homogeneo.
    U = update(U,F);
    
%   Acrescimo dos termos fontes na solucao homogenea.
    U = source(U,zL,zR);
    
    U = boundary(U);

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
plot(x(I),U(1,I)+zb(I),'b')
xlim([xmin xmax]), ylim([0 max(U(1,:))+0.1])
xlabel('x (m)'), ylabel('h+zb (m)')
title([t k]), grid on

figure(2)
plot(x(I),v(I),'b')
xlim([xmin xmax]), ylim([min(v(I))-0.1 max(v(I))+0.1])
xlabel('x (m)'), ylabel('v (m/s)')
title([t k]), grid on 

figure(3)
plot(x(I),U(2,I),'b')
xlim([xmin xmax]), ylim([min(U(2,:))-0.1 max(U(2,:))+0.1])
xlabel('x (m)'), ylabel('q (m^2/s)')
title([t k]), grid on 

toc
%-----------------------------------------------------------------------
%Inicialização do vetor de variaveis x, zb e  U
%-----------------------------------------------------------------------
function [x,zb,zL,zR,U] = start
    global dx m
    global xmin xmax hu vu 
    
    U = zeros(2,m+2); [x,zb,zL,zR] = deal(zeros(1,m+2)); 
     
    dx = (xmax-xmin)/m;
    
    x(1) = xmin;
    for i = 2:m+1
        x(i) = x(i-1)+dx;
        
        zL(i) = 0; zR(i) = 0;       
        zb(i) = 0.5*(zL(i)+zR(i));
        
        U(1,i) = hu;
        U(2,i) = hu*vu;
    end
    
    U(:,1)   = U(:,2); 
    U(:,m+2) = U(:,m+1);
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
%Fluxo Numerico utilizando os solver exatos do problema de Riemann
%-----------------------------------------------------------------------
function F = GodunovRs(U)
    global m
    
    F = zeros(2,m+2);
    
    for i = 2:m
        hL = U(1,i);
        if hL <= 0
            hL = 0;
            vL = 0;           
        else
            vL = U(2,i)/U(1,i);
        end
        
        hR = U(1,i+1);
        if hR <= 0
            hR = 0;
            vR = 0;  
        else
            vR = U(2,i+1)/U(1,i+1);
        end

        F(:,i) = RiemannSolver(hL,hR,vL,vR);
    end
end
%-----------------------------------------------------------------------
%Determinacao do problema de Riemann
%-----------------------------------------------------------------------
function F = RiemannSolver(hL,hR,vL,vR)
    global g epsilon
    
    cL = sqrt(g*hL); 
    cR = sqrt(g*hR);
    cs = (vR-vL)-2*(cL+cR);
    
    if hL <= epsilon | hR <= epsilon | cs >= 0
%       call dry solver
        [h,v] = drybed(hL,hR,vL,vR,cL,cR);
    else
%       call wet solver
        [h,v] = wetbed(hL,hR,vL,vR,cL,cR);
    end
 
    F = [h*v; h*v^2+0.5*g*h^2];
end
%-----------------------------------------------------------------------
%Solver de Riemann para leito seco (drybed)
%-----------------------------------------------------------------------
function [h,v] = drybed(hL,hR,vL,vR,cL,cR)
    global epsilon
    
    if hL <= epsilon
%       Leito seco a esquerda
        [h,v] = drylef(hL,hR,vL,vR,cL,cR,0);
    else
        if hR <= epsilon
%           Leito seco a direita
            [h,v] = dryrig(hL,hR,vL,vR,cL,cR,0);
        else
%           Leito seco no meio
            [h,v] = drymid(hL,hR,vL,vR,cL,cR,0);
        end
    end
end
%-----------------------------------------------------------------------
%Construcao da solucao para leito seco a esquerda (dry left, hL = 0)
%-----------------------------------------------------------------------
function [h,v] = drylef(hL,hR,vL,vR,cL,cR,s)
    global g
    
    sR = vR+cR;
    if s >= sR
        h = hR;
        v = vR;
    else
        sRs = vR-2*cR;
        if s >= sRs
            h = (1/3*(-vR+2*cR+s))^2/g;
            v = 1/3*(vR-2*cR+2*s);
        else
            h = hL;
            v = vL;
        end
    end
end
%-----------------------------------------------------------------------
%Construcao da solucao para leito seco a direita (dry right, hR = 0)
%-----------------------------------------------------------------------
function [h,v] = dryrig(hL,hR,vL,vR,cL,cR,s)
    global g
    
    sL = vL-cL;
    if s <= sL
        h = hL;
        v = vL;
    else
        sLs = vL+2*cL;
        if s <= sLs
            h = (1/3*(vL+2*cL-s))^2/g;
            v = 1/3*(vL+2*cL+2*s);
        else
            h = hR;
            v = vR;
        end
    end
end
%-----------------------------------------------------------------------
%Construcao da solucao para leito seco no meio (hs = 0)
%-----------------------------------------------------------------------
function [h,v] = drymid(s,hL,hR,vL,vR,cL,cR)
    global g
    
    sL = vL-cL; sLs = vL+2*cL;
    sR = vR+cR; sRs = vR-2*cR;
    
    if s <= sL
        h = hL;
        v = vL;
    end
    
    if s > sL & s <= sLs
        h = (1/3*(vL+2*cL-s))^2/g;
        v = 1/3*(vL+2*cL+2*s);
    end
    
    if s > sLs & s <= sRs
        h = 0;
        v = 0;
    end
    
    if s > sRs & s <= sR
        h = (1/3*(-vR+2*cR+s))^2/g;
        v = 1/3*(vR-2*cR+2*s);
    end
    
    if s > sR
        h = hR;
        v = vR;
    end
end
%-----------------------------------------------------------------------
%Solver de Riemann para leito molhado (wetbed)
%-----------------------------------------------------------------------
function [h,v] = wetbed(hL,hR,vL,vR,cL,cR)
    global g nint tol

%   Chute inicial para o metodo de Newton-Raphson        
    hs = 1/g*(0.5*(cL+cR)+0.25*(vL-vR))^2;
    hz = hs;
    
    for i = 1:nint
        cs = sqrt(g*hs);
        
%       Fluxo numerico para o metodo de Newton-Raphson        
        [fL,dfL] = Nflux(hL,hs,vL,cL,cs);
        [fR,dfR] = Nflux(hR,hs,vR,cR,cs);
        
        f = fR+fL+vR-vL;
        df = dfL+dfR;
        hs = hs-f/df;
        
        if abs(hs-hz) <= tol
            break
        else
            hz = hs;
        end
    end
    
    cs = sqrt(g*hs);
    vs = 0.5*(vR+vL+fR-fL);
    
    [h,v] = samwet(hL,hR,hs,vL,vR,vs,cL,cR,cs,0);
end
%-----------------------------------------------------------------------
%Avaliacao da funcao f(h*) pelo metodo de Newton
%-----------------------------------------------------------------------
function [fK,dfK] = Nflux(hK,hs,vK,cK,cs)
    global g 
    
%   Onda K(Left or Right)
    if hs <= hK
        fK = 2*(cs-cK);
        dfK = g/cs;
    else
        fK = (hs-hK)*sqrt((hs+hK)*g/(2*hs*hK));
        dfK = sqrt((hs+hK)*g/(2*hs*hK))...
            -(hs-hK)*g/(4*hs^2)*1/sqrt((hs+hK)*g/(2*hs*hK));
    end
end
%-----------------------------------------------------------------------
%Construcao da solucao para leito molhado
%-----------------------------------------------------------------------
function [h,v] = samwet(hL,hR,hs,vL,vR,vs,cL,cR,cs,s)
    global g
    
    if s <= vs
%       Onda a esquerda
        if hL <= hs
%           Onda de choque
            lbL = sqrt(0.5*(hs+hL)*hs/hL^2);
            sL = vL-lbL*cL;
            if s <= sL
                h = hL;
                v = vL;
            else
                h = hs;
                v = vs;
            end
        else
%           Onda de rarefacao
            sL = vL-cL;
            if s <= sL
                h = hL;
                v = vL;
            else
                sLs = vs-cs;
                if s <= sLs
                    h = (1/3*(vL+2*cL-s))^2/g;
                    v = 1/3*(vL+2*cL+2*s);
                else
                    h = hs;
                    v = vs;
                end
            end
        end
    else
%       Onda a direita
        if hR <= hs
%           Onda de choque
            lbR = sqrt(0.5*(hs+hR)*hs/hR^2);
            sR = vR+lbR*cR;
            if s >= sR
                h = hR;
                v = vR;
            else
                h = hs;
                v = vs;
            end
        else
%           Onda de rarefacao
            sR = vR+cR;
            if s >= sR
                h = hR;
                v = vR;
            else
                sRs = vs+cs;
                if s >= sRs
                    h = (1/3*(-vR+2*cR+s))^2/g;
                    v = 1/3*(vR-2*cR+2*s);
                else
                    h = hs; 
                    v = vs;
                end
            end
        end
    end
end
%-----------------------------------------------------------------------
%Atualizacao da solucao
%-----------------------------------------------------------------------
function U = update(U,F)
    global dx dt m
    
    for i = 3:m
        U(:,i) = U(:,i)-dt/dx*(F(:,i)-F(:,i-1));
    end
end
%-----------------------------------------------------------------------
%Termos fontes - Equacao de Manning + Derivada de zb
%-----------------------------------------------------------------------
function U = source(U,zL,zR)
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
            S(2,i) = -g*hk*sf-g*hk*(zR(i)-zL(i))/dx;
        end
        if (U(2,i)+dt*S(2,i)) < 0
            sd = 1;
        end
        U(2,i) = U(2,i)+dt*S(2,i);
    end
end
%-----------------------------------------------------------------------
%Condicoes de fronteira
%-----------------------------------------------------------------------
function U = boundary(U)
    global m t hu vu
    
    U(1,2) = hu;
    U(2,2) = hu*vu;
    U(:,1) = U(:,2);
    
    U(1,m+1) = min([0.265 0.031+0.00468*t]);
    U(2,m+1) = U(2,m);
    U(:,m+2) = U(:,m+1);
end
%-----------------------------------------------------------------------
%Geracao do arquivo de saida
%-----------------------------------------------------------------------
function output(x,zb,U,v)
    outfile = [x;U;v;zb];
    fileID = fopen('h_jump_Godunov.out','w');
    fprintf(fileID,'%6.2f %12.8f %12.8f %12.8f %12.8f\r\n',outfile);
    fclose(fileID);
end