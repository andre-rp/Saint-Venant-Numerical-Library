%   Solver Exato do problema de Riemann
% 
%   Solucao exata do problema de Riemann por amostragem, baseado na velocidade
%   das ondas que compoe o problema (choque ou rarefacao).
%   
%   Os testes suportados por esse programa sao: Escoamentos subcritico, 
%   transcritico e leito seco.
tic
clear all, close all, clc

%variaveis globais
global g dx m nint tol tmax
global xmin xmax dam

%constantes
g = 9.8; tol = 1e-6; nint = 50;  

%malha
xmin = -30; xmax = 30; m = 600;

%dados do problema
hL   = 3; hR = 0; vL = 0; vR = 0;
dam  = 0;

tmax = 2;

U = RiemannSolver(hL,hR,vL,vR);

for i = 1:m+1
    x(i) = xmin+(i-1)*dx-dam;
    if U(1,i) == 0
        v(i) = 0;
    else
        v(i) = U(2,i)/U(1,i);
    end
end

output(x,U,v);

figure(1)
plot(x,U(1,:),'b')
xlim([xmin xmax]), ylim([0 max(U(1,:))+0.5])
xlabel('x (m)'), ylabel('h+zb (m)')
title([tmax]), grid on

figure(2)
plot(x,v,'b')
xlim([xmin xmax]), ylim([min(v)-0.5 max(v)+0.5])
xlabel('x (m)'), ylabel('v (m/s)')
title([tmax]), grid on 

figure(3)
plot(x,U(2,:),'b')
xlim([xmin xmax]), ylim([min(U(2,:))-0.1 max(U(2,:))+0.1])
xlabel('x (m)'), ylabel('q (m^2/s)')
title([tmax]), grid on 

toc

%-----------------------------------------------------------------------
%Escolha do solver de Riemann - Molhado(wet) ou Seco(dry)
%-----------------------------------------------------------------------
function U = RiemannSolver(hL,hR,vL,vR)
    global g
    
    cL = sqrt(g*hL); 
    cR = sqrt(g*hR);
    cs = (vR-vL)-2*(cL+cR);
    
    if hL == 0 | hR ==0 | cs >= 0
%       call dry solver
        U = drybed(hL,hR,vL,vR,cL,cR);
    else
%       call wet solver
        U = wetbed(hL,hR,vL,vR,cL,cR);
    end
end
%-----------------------------------------------------------------------
%Solver de Riemann para leito seco (drybed)
%-----------------------------------------------------------------------
function U = drybed(hL,hR,vL,vR,cL,cR)
    global g dx m tmax
    global xmin xmax dam
    
    U = zeros(2,m+1);
    
    dx = (xmax-xmin)/m;
    
    for i = 1:m+1
        x = xmin+(i-1)*dx-dam;
        s = x/tmax;
        
        if hL == 0
%           Leito seco a esquerda
            [h,v] = drylef(hL,hR,vL,vR,cL,cR,s);
        else
            if hR == 0
%               Leito seco a direita
                [h,v] = dryrig(hL,hR,vL,vR,cL,cR,s);
            else
%               Leito seco no meio
                [h,v] = drymid(hL,hR,vL,vR,cL,cR,s);
            end
        end
        U(1,i) = h;
        U(2,i) = h*v;
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
function U = wetbed(hL,hR,vL,vR,cL,cR)
    global g dx m nint tol tmax
    global xmin xmax dam
    
    U = zeros(2,m+1);
    
    dx = (xmax-xmin)/m;
    
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
    
%   Solucao por amostragem    
    for i = 1:m+1
        x = xmin+(i-1)*dx-dam;
        s = x/tmax;
        
        [h,v] = samwet(hL,hR,hs,vL,vR,vs,cL,cR,cs,s);
        
        U(1,i) = h;
        U(2,i) = h*v;
    end
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
%Geracao do arquivo de saida
%-----------------------------------------------------------------------
function output(x,U,v)
    outfile = [x;U;v];
    fileID = fopen('db_ERiemann.out','w');
    fprintf(fileID,'%6.2f %12.8f %12.8f %12.8f\r\n',outfile);
    fclose(fileID);
end