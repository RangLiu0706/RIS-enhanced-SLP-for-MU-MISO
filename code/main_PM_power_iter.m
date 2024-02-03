% This Matlab script can be used to generate Fig. 4 in the paper:
% R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “Joint symbol-level precoding and reflecting designs for IRS-enhanced MU-MISO systems,” IEEE Trans. Wireless Commun., vol. 20, no. 2, pp. 798-811, Feb. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9219206
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02

clear;
clc;

M = 6;
N = 64;
K = 3;
SNR = 10;
B1 = 1;
B2 = 2;
B3 = 3;
global sigma2;
sigma2 = 1e-11;

omega = 4;  %%%% modulation order

N_sim = 1000;
iter_range = (1:1:14); 

power_my = zeros(1,length(iter_range));
power_my_1 = zeros(1,length(iter_range));
power_my_2 = zeros(1,length(iter_range));
power_my_3 = zeros(1,length(iter_range));

d_ar = 50;
d_ru = 3;
belta1 = sqrt(10^(0.3)/(1+10^(0.3)));
belta2 = sqrt(1/(1+10^(0.3)));

Nmax = 14;
res_th = -1;

for sim = 1:N_sim
    d_au = (d_ar+d_ru) + 2*d_ru*rand(1,K);
    H_au = zeros(K,M);
    for uk = 1:K
        h_LOS = exp(1i*pi*sin(pi*rand-pi/2))*exp(1i*pi*(0:1:M-1)'*sin(pi*rand-pi/2))';
        H_au(uk,:) = sqrt(10^(-3)*d_au(uk)^(-3.5))*(belta1*h_LOS+belta2*(randn(1,M)+1i*randn(1,M))/sqrt(2));
    end
    H_LOS = exp(1i*pi*(0:1:N-1)'*sin(pi*rand-pi/2))*exp(1i*pi*(0:1:M-1)'*sin(pi*rand-pi/2))';    
    H_ar = sqrt(10^(-3)*d_ar^(-2.8))*(belta1*H_LOS+belta2*(randn(N,M)+1i*randn(N,M))/sqrt(2));
    H_LOS = exp(1i*pi*(0:1:K-1)'*sin(pi*rand-pi/2))*exp(1i*pi*(0:1:N-1)'*sin(pi*rand-pi/2))';  
    H_ru = sqrt(10^(-3)*d_ru^(-2.5))*(belta1*H_LOS+belta2*(randn(K,N)+1i*randn(K,N))/sqrt(2));
        
    [~,~,p_my] = getX_my_PM(H_au,H_ar,H_ru,SNR,Nmax,res_th,omega);
    power_my = power_my + 10*log10(1000.*p_my./omega^K);
    
    [~,~,p_my_1] = getX_my_PM_b(H_au,H_ar,H_ru,SNR,B1,Nmax,res_th,omega);
    power_my_1 = power_my_1 + 10*log10(1000.*p_my_1./omega^K);
    
    [~,~,p_my_2] = getX_my_PM_b(H_au,H_ar,H_ru,SNR,B2,Nmax,res_th,omega);
    power_my_2 = power_my_2 + 10*log10(1000.*p_my_2./omega^K);
    
    [~,~,p_my_3] = getX_my_PM_b(H_au,H_ar,H_ru,SNR,B3,Nmax,res_th,omega);
    power_my_3 = power_my_3 + 10*log10(1000.*p_my_3./omega^K);  

end
 
power_my = power_my/sim;
power_my_1 = power_my_1/sim;
power_my_2 = power_my_2/sim;
power_my_3 = power_my_3/sim;

figure
plot(iter_range,power_my,'-o','color',[0.85,0.1,0.1],'LineWidth',1.5)
hold on
plot(iter_range,power_my_1,'-d','color',[0.1,0.1,0.1],'LineWidth',1.5)
plot(iter_range,power_my_2,'-s','color',[0.54,0.3,0.35],'LineWidth',1.5)
plot(iter_range,power_my_3,'-^','color',[0,0.5,0],'LineWidth',1.5)
hold off
xlabel('Number of iterations');
ylabel('Average transmit power (dBm)');
grid on
legend('Proposed, B = \infty','Proposed, B = 1','Proposed, B = 2', ...
    'Proposed, B = 3');

