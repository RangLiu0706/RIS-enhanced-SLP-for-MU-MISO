% Solve for the problem (6).
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, "Joint symbol-level precoding and reflecting designs for IRS-enhanced MU-MISO systems,” IEEE Trans. Wireless Commun., vol. 20, no. 2, pp. 798-811, Feb. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9219206
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02

function [X,theta,power] = getX_my_PM(H_au,H_ar,H_ru,SNR,Nmax,res_th,omega)
[K,~] = size(H_au);
[~,N] = size(H_ru);
a = zeros(K*omega^K,1);
B = zeros(K*omega^K,N);
phi_u = zeros(K,omega^K);
[Theta,theta] = get_theta_initial_RCG(H_ru,H_au,H_ar);
X = getX_wo_IRS(H_au+H_ru*Theta*H_ar,SNR,omega);

for k = 0:1:omega^K-1

    s = dec2bin(k,log2(omega)*K);
    for i = 1:1:K
        t_ind = (0:1:log2(omega)-1);
        temp = sum( 2.^(t_ind).*(s(log2(omega)*i-t_ind)-48) );
        phi_u(i,k+1) = pi/omega + temp*2*pi/omega;
    end
end

res = 1;
iter = 1;
p_n = norm(X,'fro')^2;
while iter <= Nmax && res >= res_th
    p_re = p_n;

    for ii = 1:omega^K
        for jj = 1:K
            a((ii-1)*K+jj,1) = H_au(jj,:)*X(:,ii)*exp(-1i*phi_u(jj,ii));
            B((ii-1)*K+jj,:) = H_ru(jj,:)*diag(H_ar*X(:,ii)*exp(-1i*phi_u(jj,ii)));
        end
    end

    [theta] = get_theta_sum(a,B,K,theta,omega);
    Theta = diag(theta);

    if isnan(theta)
        theta = exp(1i*2*pi*randn(N,1));
        Theta = diag(theta);
        iter = 1;
    end

    [X] = getX_wo_IRS(H_au+H_ru*Theta*H_ar,SNR,omega);

    p_n = norm(X,'fro')^2;
    res = abs(1-p_n/p_re);
    power(1,iter) = p_n;
    iter = iter + 1;
end
end
