% Solve for the symbol-level precoding without RIS.
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, "Joint symbol-level precoding and reflecting designs for IRS-enhanced MU-MISO systems,” IEEE Trans. Wireless Commun., vol. 20, no. 2, pp. 798-811, Feb. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9219206
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02

function [X] = getX_wo_IRS(H,SNR,omega)
[M,N] = size(H);
phi = pi/omega;
phi_u = zeros(M,omega^M);
% X = zeros(N,omega^M);
% g = sqrt(10^(-11)*10^(0.1*SNR))*ones(M,1);
% Delta1 = [eye(N),zeros(N,N);zeros(N,N),-eye(N)];
% Delta2 = [zeros(N,N),eye(N);eye(N),zeros(N,N)];
% H_bar = zeros(M,2*N);
% c = [g;g];
% Nmax = 10000;
% res_th = 1e-20;
for k = 0:1:omega^M-1

    s = dec2bin(k,log2(omega)*M);
    for i = 1:1:M
        t_ind = (0:1:log2(omega)-1);
        temp = sum( 2.^(t_ind).*(s(log2(omega)*i-t_ind)-48) );
        phi_u(i,k+1) = pi/omega + temp*2*pi/omega;
    end
    % % %     %%%%%%%%%%%%%%%%%gd
    %     for i = 1:M
    %         H_bar(i,:) = [real(H(i,:).*exp(-1i*phi_u(i,k+1))),imag(H(i,:).*exp(-1i*phi_u(i,k+1)))];
    %     end
    %     A = [H_bar*(Delta1.*tan(phi)-Delta2);H_bar*(Delta1.*tan(phi)+Delta2)];
    %
    %     %iter
    %     lambda = rand(2*M,1);
    %     iter = 0;
    %     res = 1;
    %     while iter < Nmax && res >= res_th
    %         pre = lambda;
    %         iter = iter + 1;
    %         gf = 0.5*A*A.'*lambda-c;
    %         a = (c.'*c-0.5*c.'*A*A.'*lambda-0.5*lambda.'*A*(A.'*c-0.5*A.'*A*A.'*lambda))/((-0.25*A.'*A*A.'*lambda+0.5*A.'*c).'*(A.'*c-0.5*A.'*A*A.'*lambda));
    %         lambda = max(lambda-a*gf,zeros(2*M,1));
    %         res = norm(lambda-pre,2)/norm(lambda,2);
    %     end
    %     p = 0.5*A.'*lambda;
    %     X(:,k+1) = p(1:N)+1i*p(N+1:2*N);
end
cvx_begin quiet
variable X(N,omega^M) complex
minimize square_pos(norm(X,'fro'))
subject to
r = H*X.*exp(-1i*phi_u)./sqrt(10^-11);
real(r)*tan(phi)-abs(imag(r))-sqrt(10^(0.1*SNR)) >= 0;
cvx_end
end
