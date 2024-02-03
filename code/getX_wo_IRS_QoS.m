% Solve for the symbol-level precoding without RIS.
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, "Joint symbol-level precoding and reflecting designs for IRS-enhanced MU-MISO systems,” IEEE Trans. Wireless Commun., vol. 20, no. 2, pp. 798-811, Feb. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9219206
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02

function [X,t] = getX_wo_IRS_QoS(H,P,omega)

[M,N] = size(H);
phi_u = zeros(M,omega^M);
SNR = 5;
t0 = sqrt(10^(-11)*10^(0.1*SNR));
A = zeros(omega^M,omega^M);
X = zeros(N,omega^M);
for k = 0:1:omega^M-1
    s = dec2bin(k,log2(omega)*M);
    for i = 1:1:M
        t_ind = (0:1:log2(omega)-1);
        temp = sum( 2.^(t_ind).*(s(log2(omega)*i-t_ind)-48) );
        phi_u(i,k+1) = pi/omega + temp*2*pi/omega;
    end
end  

[X1] = getX_wo_IRS(H,SNR,omega);
for i = 1:1:omega^M
    A(i,i) = t0/(norm(X1(:,i),2));
end
[y] = get_y_powerAllocate(A,t0.*ones(omega^M,1)); 
y = sqrt(P*omega^M)/norm(y,2).*y;
p = y.*y;
for i = 1:1:omega^M
    X(:,i) = (sqrt(p(i)))/(norm(X1(:,i),2)).*X1(:,i);
end

t = min(min(real(H*X.*exp(-1i*phi_u))-abs(imag(H*X.*exp(-1i*phi_u)))));

