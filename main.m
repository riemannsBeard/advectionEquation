clc
clear all
close all

c = 1;
L = 1;
h = 1e-4;
x = 0:h:L;
N = length(x);
CFL = 0.5;
dt = CFL*h/c;
tf = 0.5;

%% Matrices de diferenciacion

A1 = speye(N);
A1(1,1) = 0;

A2 = spdiags([ones(N,1) -ones(N,1)], [0 -1], N, N);
A2(1,1) = 0;

%% Inicializacion

u0 = sparse(exp(-10*(4*x - 1).^2)');
u0(1) = 0;

%% Simulacion

% Primer paso temporal (metodo de Euler)
u1 = sparse(-c*dt/h*A2*u0 + u0);

uold2 = u0;
uold1 = u1;

figure(1)
plot(x, u1)
ylim([0 1])

for t = 2*dt:dt:tf
    
    u = A1*uold1 - 1.5*c*(dt/h)*A2*uold1 + 0.5*c*(dt/h)*A2*uold2;
    u(1) = 0;
    
    uold2 = uold1;
    uold1 = u;
    
    t = t + dt;
    
%     figure(1)
%     plot(x, u)
%     ylim([0 1])
%     
%     pause(0.0625)
%     drawnow
end

figure(1)
plot(x, u)
set(gca, 'FontSize', 16, 'fontname', 'times')
ylim([0 1])

%% Solucion teorica
uTeo = fTeo(x, c, t);

figure(1),
hold on
plot(x, uTeo)

h = legend('Adams-Bashforth', 'Theoretical', 'location', 'best');
h.FontSize = 16;