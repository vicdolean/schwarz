clear all, close all, clc

k = 30;
L = 1;
delta = L/10;
theta = -pi:1/50:pi;
sigma = 0.1;
% sigma = 1;

n = [2,5:5:80];
j = 1;
freq = 1:floor(2*k*L/pi);
for ni = n
    i = 1;
    rhoest = 0;
    for m = freq
        kt = m*pi/L;
        lambda = sqrt(kt^2-k^2+1i*k*sigma);
        Dt = (lambda + 1i*k + sigma)^2*exp(lambda*(2*delta+L))-(lambda-1i*k-sigma)^2*exp(-lambda*(2*delta+L));
        a = ((lambda + 1i*k + sigma)^2*exp(2*lambda*delta) - (lambda-1i*k-sigma)^2*exp(-2*lambda*delta))/Dt;
        b = -(lambda^2 -(1i*k+sigma)^2)*(exp(lambda*L) - exp(-lambda*L))/Dt;
        rhoestkt(i) = max(max(abs(a-b),abs(a+b)));
        rhoest = max(rhoest,max(abs(a-b),abs(a+b)));
        L1 = (cos(theta)*a + sqrt(-sin(theta).^2*a^2 + b^2));
        L2 = (cos(theta)*a - sqrt(-sin(theta).^2*a^2 + b^2));
        [M,S] = iteration_matrix(ni,a,b);
        rho(i) = max(abs(S));
        i = i+1;
%         figure(10)
%         plot(real(S),imag(S),'bx',real(L1),imag(L1),'r-',real(L2),imag(L2),'r-'), grid on
%         legend('Spectrum','Theoretical estimate')
%         pause(0.01)
    end
    rhomax(j) = max(rho);
    j = j+1;
end

figure(1)
plot(freq*pi/L,rho,'bx-','LineWidth',2,'MarkerSize',12)
% plot(freq*pi/L,rho,'bx-',freq*pi/L,rhoestkt,'rx-','LineWidth',2,'MarkerSize',12)
% legend({'Convergence factor','Theoretical limit'},'FontSize',32,'Interpreter','LaTeX','Location','NorthEast')
xlabel('Fourier number $\tilde k$','FontSize',36,'Interpreter','LaTeX')
ylabel('Convergence factor of the Fourier mode','FontSize',36,'Interpreter','LaTeX')
title('Convergence factor of the Schwarz algorithm','FontSize',36,'Interpreter','LaTeX')
grid on, axis square
set(gca,'FontSize',22)
saveas(gcf,'FourierModeConvergenceMaxwell2D','epsc')

figure(2)
plot(n,rhomax,'bx-',[0,max(n)],[rhoest,rhoest],'r-','LineWidth',2,'MarkerSize',12)
legend({'Convergence factor','Theoretical limit'},'FontSize',32,'Interpreter','LaTeX','Location','SouthEast')
xlabel('Number of subdomains $N$','FontSize',36,'Interpreter','LaTeX')
ylabel('Convergence factor','FontSize',36,'Interpreter','LaTeX')
title('Convergence factor of the Schwarz algorithm','FontSize',36,'Interpreter','LaTeX')
grid on, axis square
set(gca,'FontSize',22)
saveas(gcf,'ConvergenceMaxwell2D','epsc')