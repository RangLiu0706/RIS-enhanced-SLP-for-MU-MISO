% Initialize phi using the RCG algorithm.
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, "Joint symbol-level precoding and reflecting designs for IRS-enhanced MU-MISO systems,” IEEE Trans. Wireless Commun., vol. 20, no. 2, pp. 798-811, Feb. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9219206
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02

function [Theta,theta] = get_theta_initial_RCG(H_ru,H_au,H_ar)
[K,N] = size(H_ru);

A1 = zeros(N,N,K);
A2 = zeros(K,N);
A3 = zeros(K,N);
A4 = zeros(N,N,K);
a5 = zeros(1,K);
for k = 1:1:K
    ak = H_au(k,:);
    Bk = diag(H_ru(k,:))*H_ar;
    A1(:,:,k) = real(Bk)*real(Bk)' + imag(Bk)*imag(Bk)';
    A2(k,:) = 2*real(ak)*real(Bk)' + 2*imag(ak)*imag(Bk)';
    A3(k,:) = 2*imag(ak)*real(Bk)' - 2*real(ak)*imag(Bk)';
    A4(:,:,k) = imag(Bk)*real(Bk)' - real(Bk)*imag(Bk)';
    a5(k) = norm(ak,2)^2;
end

% Create the problem structure.
manifold = obliquefactory(2, N);
problem.M = manifold;

problem.cost = @cost_t;
        function f = cost_t(x)
            fk = zeros(1,K);
            for ik = 1:1:K
                fk(ik) = x(1,:)*A1(:,:,ik)*x(1,:)' + x(2,:)*A1(:,:,ik)*x(2,:)' ...
                    + A2(ik,:)*x(1,:)' + A3(ik,:)*x(2,:)' + 2*x(1,:)*A4(:,:,ik)*x(2,:)' + a5(ik);
            end
            f = epsl*log(sum(exp(-fk./epsl)));
        end
    problem.grad = @(x) problem.M.egrad2rgrad(x,egrad(x));
        function g = egrad(x)
            fk = zeros(1,K);
            for ik = 1:1:K
                fk(ik) = x(1,:)*A1(:,:,ik)*x(1,:)' + x(2,:)*A1(:,:,ik)*x(2,:)' ...
                    + A2(ik,:)*x(1,:)' + A3(ik,:)*x(2,:)' + 2*x(1,:)*A4(:,:,ik)*x(2,:)' + a5(ik);
            end
            Fp1 = zeros(K,N);
            Fp2 = zeros(K,N);
            for ik = 1:1:K
                Fp1(ik,:) = 2*x(1,:)*A1(:,:,ik)' + A2(ik,:) + 2*x(2,:)*A4(:,:,ik);
                Fp2(ik,:) = 2*x(2,:)*A1(:,:,ik)' + A3(ik,:) + 2*x(1,:)*A4(:,:,ik);
            end
            f3 = -sum(repmat(exp(-fk'./epsl),1,N).*Fp1)/sum(exp(-fk./epsl));
            f4 = -sum(repmat(exp(-fk'./epsl),1,N).*Fp2)/sum(exp(-fk./epsl));
            g = [f3;f4];
        end
    
    % Execute the optimization
    epsl= 1e-10;  %%%%%%%%%%%
        
    theta = ones(N,1);
    options.tolgradnorm = 1e-20;
    options.maxiter = 1000; %%%
    options.minstepsize = 1e-20;
    options.verbosity = 0;

    [x,~,~] = conjugategradient(problem,[real(theta.');imag(theta.')],options);
    theta = (x(1,:)+1i*x(2,:)).';
    Theta = diag(theta);

end