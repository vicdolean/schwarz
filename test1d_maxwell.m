clear all, close all, clc

k = 30;
L = 1;
delta = L/10;
sigma = 0.1;
% sigma = 5;

lambda = sqrt(1i*k*sigma-k^2);
Dt = (lambda +1i*k)^2*exp(lambda*(2*delta+L))-(lambda-1i* k)^2*exp(-lambda*(2*delta+L));
a = ((lambda+1i*k)^2*exp(2*lambda*delta) - (lambda - 1i*k)^2*exp(-2*lambda*delta))/Dt;
b = -1i*k*sigma*(exp(lambda*L) - exp(-lambda*L))/Dt;

theta = -pi:1/50:pi;
n = 5:5:160;
L1 = cos(theta)*a + sqrt(-sin(theta).^2*a^2 + b^2);
L2 = cos(theta)*a - sqrt(-sin(theta).^2*a^2 + b^2);
rhoest = max(abs(a-b),abs(a+b));

figure(10)
i = 1;
for ni = n
    [M,S] = iteration_matrix(ni,a,b);
    rho(i) = max(abs(S));
    i = i+1;
    plot(real(S),imag(S),'bx', real(L1),imag(L1),'r-',real(L2),imag(L2),'r-')
    legend('Spectrum of T1d','Theoretical estimate')
    title('Spectrum of the Schwarz iteration matrix vs. theoretical estimate')
    grid on
    pause(0.1)
end

figure(1)
plot(real(S),imag(S),'bx','MarkerSize',12), hold on
plot(real(L1),imag(L1),'r-',real(L2),imag(L2),'r-','LineWidth',2), hold off
legend({'Spectrum','Theoretical limit'},'FontSize',32,'Interpreter','LaTeX','Location','Best')
xlabel('$\mathrm{Re}(\lambda)$','FontSize',36,'Interpreter','LaTeX')
ylabel('$\mathrm{Im}(\lambda)$','FontSize',36,'Interpreter','LaTeX')
title('Spectrum of the Schwarz iteration matrix','FontSize',36,'Interpreter','LaTeX')
grid on, axis square
set(gca,'FontSize',22)
saveas(gcf,'Spectrum1D','epsc')

figure(2)
plot(n,rho,'bx-',[0,max(n)],[rhoest,rhoest],'r-','LineWidth',2,'MarkerSize',12)
legend({'Convergence factor','Theoretical limit'},'FontSize',32,'Interpreter','LaTeX','Location','SouthEast')
xlabel('Number of subdomains $N$','FontSize',36,'Interpreter','LaTeX')
ylabel('Convergence factor','FontSize',36,'Interpreter','LaTeX')
title('Convergence factor of the Schwarz algorithm','FontSize',36,'Interpreter','LaTeX')
grid on, axis square
set(gca,'FontSize',22)
saveas(gcf,'Convergence1D','epsc')