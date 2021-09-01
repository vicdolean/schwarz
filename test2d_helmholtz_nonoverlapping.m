clear all, close all, clc

k = 30;
L = 1;
delta = 0;
sigma = 0.1;
% sigma = 1;

n = 80;
freq = 1:floor(2*k*L/pi);
i = 1;
rhoestkt = 0;
for m = freq
    kt = m*pi/L;
    lambda = sqrt(kt^2-k^2+1i*k*sigma);
    Dt = (lambda + 1i*k)^2*exp(lambda*(2*delta+L))-(lambda-1i*k)^2*exp(-lambda*(2*delta+L));
    a = ((lambda + 1i*k )^2*exp(2*lambda*delta) - (lambda-1i*k)^2*exp(-2*lambda*delta))/Dt;
    b = -(lambda^2 -(1i*k)^2)*(exp(lambda*L) - exp(-lambda*L))/Dt;
    rhoestkt(i) = max(max(abs(a-b),abs(a+b)),abs(a));
    [M,S] = iteration_matrix(n,a,b);
    rho(i) = max(abs(S));
    i = i+1;
end

figure(1)
plot(freq*pi/L,rho,'bx-',freq*pi/L,rhoestkt,'r:','LineWidth',2,'MarkerSize',12)
legend({'Convergence factor','Theoretical limit'},'FontSize',32,'Interpreter','LaTeX','Location','SouthEast')
xlabel('Fourier number $\tilde k$','FontSize',36,'Interpreter','LaTeX')
ylabel('Convergence factor of the Fourier mode','FontSize',36,'Interpreter','LaTeX')
title('Convergence factor of the Schwarz algorithm','FontSize',36,'Interpreter','LaTeX')
grid on, axis square
set(gca,'FontSize',22)
saveas(gcf,'FourierModeConvergenceHelmholtz2D','epsc')