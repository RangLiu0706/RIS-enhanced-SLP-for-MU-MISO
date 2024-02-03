% This Matlab script can be used to generate Fig. 13 in the paper:
% R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “Joint symbol-level precoding and reflecting designs for IRS-enhanced MU-MISO systems,” IEEE Trans. Wireless Commun., vol. 20, no. 2, pp. 798-811, Feb. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9219206
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02

clear;
clc;

M = 6;
N = 64;
K = 3;
power_t = 15;
power = 10^(0.1*power_t-3);
global sigma2;
sigma2 = 1e-11;

omega = 4;  %%% modulation order

N_sim = 1000;
d_range = (45:1:55);

SER_my = zeros(1,length(d_range));
SER_my_1 = zeros(1,length(d_range));
SER_my_2 = zeros(1,length(d_range));
SER_my_3 = zeros(1,length(d_range));

B1 = 1;
B2 = 2;
B3 = 3;
d_ar = 50;
c_ru = 3;
belta1 = sqrt(10^(0.3)/(1+10^(0.3)));
belta2 = sqrt(1/(1+10^(0.3)));
dt = 0.5;

Nmax = 20;
res_th = 1e-3;

Ns = 10000;

for sim = 1:N_sim

    H_LOS = exp(1i*pi*(0:1:N-1)'*sin(pi*rand-pi/2))*exp(1i*pi*(0:1:M-1)'*sin(pi*rand-pi/2))';
    H_ar = sqrt(10^(-3)*d_ar^(-2.8))*(belta1*H_LOS+belta2*(randn(N,M)+1i*randn(N,M))/sqrt(2));

    for d_index = 1:length(d_range)
        d = d_range(d_index);

        d_ru = [3,3,sqrt(dt^2+(50-d)^2)];
        d_au = [(d_ar+c_ru)+2*c_ru*rand(1,K-1),sqrt(d^2+dt^2)];
        H_au = zeros(K,M);
        H_ru = zeros(K,N);
        for uk = 1:K
            h_LOS = exp(1i*pi*sin(pi*rand-pi/2))*exp(1i*pi*(0:1:M-1)'*sin(pi*rand-pi/2))';
            H_au(uk,:) = sqrt(10^(-3)*d_au(uk)^(-3.5))*(belta1*h_LOS+belta2*(randn(1,M)+1i*randn(1,M))/sqrt(2));
            h_LOS2 = exp(1i*pi*sin(pi*rand-pi/2))*exp(1i*pi*(0:1:N-1)'*sin(pi*rand-pi/2))';
            H_ru(uk,:) = sqrt(10^(-3)*d_ru(uk)^(-2.5))*(belta1*h_LOS2+belta2*(randn(1,N)+1i*randn(1,N))/sqrt(2));
        end

        Noise = sqrt(0.5*sigma2)*(randn(K,Ns)+1i*randn(K,Ns));
        S = exp(1i*2*pi*rand(K,Ns));
        [S_index,index] = get_adaptive_modulate(S,omega*ones(1,K));

        [X_my,theta_my,t_my] = getX_my_QoS(H_au,H_ar,H_ru,power,Nmax,res_th,omega);
        R_my = (H_au + H_ru*diag(theta_my)*H_ar)*X_my(:,index) + Noise;
        SER_my(d_index) = SER_my(d_index) + sum(get_SER(R_my,S_index,omega*ones(1,K)))/K;

        [X_my_1,theta_my_1,t_1] = getX_my_QoS_b(H_au,H_ar,H_ru,power,B1,Nmax,res_th,omega);
        R_my_1 = (H_au + H_ru*diag(theta_my_1)*H_ar)*X_my_1(:,index) + Noise;
        SER_my_1(d_index) = SER_my_1(d_index) + sum(get_SER(R_my_1,S_index,omega*ones(1,K)))/K;

        [X_my_2,theta_my_2,t_2] = getX_my_QoS_b(H_au,H_ar,H_ru,power,B2,Nmax,res_th,omega);
        R_my_2 = (H_au + H_ru*diag(theta_my_2)*H_ar)*X_my_2(:,index) + Noise;
        SER_my_2(d_index) = SER_my_2(d_index) + sum(get_SER(R_my_2,S_index,omega*ones(1,K)))/K;

        [X_my_3,theta_my_3,t_3] = getX_my_QoS_b(H_au,H_ar,H_ru,power,B3,Nmax,res_th,omega);
        R_my_3 = (H_au + H_ru*diag(theta_my_3)*H_ar)*X_my_3(:,index) + Noise;
        SER_my_3(d_index) = SER_my_3(d_index) + sum(get_SER(R_my_3,S_index,omega*ones(1,K)))/K;

    end
end

SER_my = SER_my/sim;
SER_my_1 = SER_my_1/sim;
SER_my_2 = SER_my_2/sim;
SER_my_3 = SER_my_3/sim;

figure
semilogy(d_range,SER_my,'-o','color',[0.85,0.1,0.1],'LineWidth',1.5)
hold on
semilogy(d_range,SER_my_1,'-x','color',[0.1,0.1,0.1],'LineWidth',1.5,'markersize',7)
semilogy(d_range,SER_my_2,'-s','color',[0.54,0.3,0.35],'LineWidth',1.5)
semilogy(d_range,SER_my_3,'-^','color',[0,0.5,0],'LineWidth',1.5)
hold off
xlabel('{\itd} (m)');
ylabel('Average SER');
grid on
legend('Proposed, B = \infty','Proposed, B = 1','Proposed, B = 2', 'Proposed, B = 3');
