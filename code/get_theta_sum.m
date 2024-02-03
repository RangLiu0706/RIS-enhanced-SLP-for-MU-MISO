% Solve for the theta.
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, "Joint symbol-level precoding and reflecting designs for IRS-enhanced MU-MISO systems,” IEEE Trans. Wireless Commun., vol. 20, no. 2, pp. 798-811, Feb. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9219206
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02

function [theta] = get_theta_sum(a,B,K,theta,omega)
    [~,N] = size(B);
    
    phi = pi/omega;
    % Generate the problem data.
    ar = real(a);
    ai = imag(a);
    BR=real(B.');
    BI=imag(B.');
    for m = 1:1:K*4^K
        d1(:,m)=BI(:,m)-BR(:,m).*tan(phi);
        d2(:,m)=BR(:,m)+BI(:,m).*tan(phi);
        d3(:,m) = -BI(:,m)-BR(:,m)*tan(phi);
        d4(:,m) = BI(:,m)*tan(phi)-BR(:,m);
        d5(1,m) = ai(m,1)-ar(m,1);
        d6(1,m) = -ai(m,1)-ar(m,1);
    end
    
    % Create the problem structure.
    manifold = obliquefactory(2, N);
    problem.M = manifold;    
    % Define the problem cost function and its gradient.

    problem.cost = @cost_t;
        function f=cost_t(x)
            f1 = d5 + x(1,:)*d1 + x(2,:)*d2;
            f2 = d6 + x(1,:)*d3 + x(2,:)*d4;
            f_temp = epsl*log(exp(f1/epsl)+exp(f2/epsl));
            f = sum(f_temp);
        end
    problem.grad = @(x) problem.M.egrad2rgrad(x,egrad(x));
        function g = egrad(x)
            f1 = d5 + x(1,:)*d1 + x(2,:)*d2;
            f2 = d6 + x(1,:)*d3 + x(2,:)*d4;
            f5 = 1./(exp(1/epsl.*f1.')+exp(1/epsl.*f2.'));
            f6 = bsxfun(@times,exp(1/epsl.*f1.'),d1')+bsxfun(@times,exp(1/epsl.*f2.'),d3');
            f3 = sum(bsxfun(@times,f5,f6));
            
            f5 = 1./(exp(1/epsl.*f1.')+exp(1/epsl.*f2.'));
            f6 = bsxfun(@times,exp(1/epsl.*f1.'),d2')+bsxfun(@times,exp(1/epsl.*f2.'),d4');
            f4 = sum(bsxfun(@times,f5,f6));
            g = [f3;f4];
        end
    
    % Execute the optimization
    epsl= 1e-7;   
        
    options.tolgradnorm = 1e-10;
    options.maxiter = 200;  
    options.minstepsize = 1e-10;
    options.verbosity = 0;
    [x,~,~] = conjugategradient(problem,[real(theta.');imag(theta.')],options);
    theta = (x(1,:)+1i*x(2,:)).';
end

