% Solve for the power allocation.
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, "Joint symbol-level precoding and reflecting designs for IRS-enhanced MU-MISO systems,” IEEE Trans. Wireless Commun., vol. 20, no. 2, pp. 798-811, Feb. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9219206
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02

function [y] = get_y_powerAllocate(A,c)
[a_row,~] = size(A);
Nmax = 10000;
res_th = 1e-10;

lambda = rand(a_row,1);
iter = 0;
res = 1;
while iter < Nmax && res >= res_th
    pre = lambda;
    iter = iter + 1;
    gf = 0.5*A*A.'*lambda-c;
    a = (c.'*c-0.5*c.'*A*A.'*lambda-0.5*lambda.'*A*(A.'*c-0.5*A.'*A*A.'*lambda))/((-0.25*A.'*A*A.'*lambda+0.5*A.'*c).'*(A.'*c-0.5*A.'*A*A.'*lambda));
    lambda = max(lambda-a*gf,zeros(a_row,1));
    res = norm(lambda-pre,2)/norm(lambda,2);
end
y = 0.5*A.'*lambda;

