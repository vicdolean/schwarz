clear all, close all, clc

k = 30;
sigma = 0.1;

L = 1;
delta = L/10;

lambda = sqrt(1i*k*sigma-k^2);
Dt = (lambda +1i*k)^2*exp(lambda*(2*delta+L))-(lambda-1i* k)^2*exp(-lambda*(2*delta+L));
a= ((lambda+1i*k)^2*exp(2*lambda*delta) - (lambda - 1i*k)^2*exp(-2*lambda*delta))/Dt;
b= -1i*k*sigma*(exp(lambda*L) - exp(-lambda*L))/Dt;

theta = -pi:1/50:pi; n=5:10:160;
L1 = cos(theta)*a + sqrt(-sin(theta).^2*a^2 + b^2);
L2 = cos(theta)*a - sqrt(-sin(theta).^2*a^2 + b^2);
rhoest = max(abs(a-b),abs(a+b));

i=1;
for ni = n
    [M,S] = iteration_matrix(ni,a,b);
    rho(i) = max(abs(S));
    i = i+1;
    plot(real(S),imag(S),'bx', real(L1),imag(L1),'r-',real(L2),imag(L2),'r-')
    legend('Spectrum of T1d','Theoretical estimate')
    title('Spectrum of the Schwarz iteration matrix vs. theoretical estimate')
    grid on
end
figure(1)
plot(real(S),imag(S),'bx', real(L1),imag(L1),'r-',real(L2),imag(L2),'r-')
legend('Spectrum of T1d','Theoretical estimate')
title('Spectrum of the Schwarz iteration matrix vs. theoretical estimate')
grid on
saveas(gcf,'Spectrum_Schwarz','epsc')

figure(2)
plot(n,rho,'b*-',n,rhoest*ones(size(n)),'r-')
legend('Convergence factor','Limiting spectral radius')
xlabel('Number of subdomains')
ylabel('Convergence factor')
title('Convergence of the algorithm for different number of subdomains')
grid on
saveas(gcf,'Conv_Schwarz','epsc')