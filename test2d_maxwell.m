clear all, close all, clc
omega = 30;
sigma = 1;
L = 1;
delta = L/10;
theta = -pi:1/50:pi;

n=[2 5 10 20 40 80];
j=1;
freq = 0:2:20;
for ni = n
    i=1;
    ni
    rhoest = 0;
    for m= freq
        kt = m*pi/L;
        lambda=sqrt(kt^2-omega^2+1i*omega*sigma);
        Dt =(lambda + 1i*omega + sigma)^2*exp(lambda*(2*delta+L))-(lambda-1i*omega-sigma)^2*exp(-lambda*(2*delta+L));
        a= ((lambda + 1i*omega + sigma)^2*exp(2*lambda*delta) - (lambda-1i*omega-sigma)^2*exp(-2*lambda*delta))/Dt;
        b= -(lambda^2 -(1i*omega+sigma)^2)*(exp(lambda*L) - exp(-lambda*L))/Dt;
        rhoest = max(rhoest,max(abs(a-b),abs(a+b)));
        L1 = (cos(theta)*a + sqrt(-sin(theta).^2*a^2 + b^2));
        L2 = (cos(theta)*a - sqrt(-sin(theta).^2*a^2 + b^2));
        [M,S] = iteration_matrix(ni,a,b);
        S = S;
        rho(i) = max(abs(S));
        i = i+1;
        plot(real(S),imag(S),'bx',real(L1),imag(L1),'r-',real(L2),imag(L2),'r-'), grid on
        legend('Spectrum','Theoretical estimate')
    end
    rhomax(j) = max(rho);
    j = j+1;
end
figure(1)
plot(freq*pi/L,rho,'b*-'), grid on
xlabel('Frequency')
ylabel('Convergence factor')
saveas(gcf,'Covergence_factor','epsc')
figure(2)
plot(n,rhomax,'b*-',n,rhoest*ones(size(n)),'r-'), grid on
legend('Convergence factor','Limiting spectral radius')
xlabel('Number of domains')
ylabel('Convergence factor');
saveas(gcf,'Conv_Schwarz','epsc')