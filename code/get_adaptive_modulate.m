% Modulated the symbols.
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, "Joint symbol-level precoding and reflecting designs for IRS-enhanced MU-MISO systems,” IEEE Trans. Wireless Commun., vol. 20, no. 2, pp. 798-811, Feb. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9219206
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02

function [S,index] = get_adaptive_modulate(S,omega)

[K,N_s] = size(S);
index = zeros(1,N_s);

for ik = 1:1:K
    mpsk = exp( 1i*(pi/omega(ik)+2*pi/omega(ik)*(0:1:omega(ik)-1)) );
    Distance = abs(repmat(S(ik,:),omega(ik),1)-repmat(mpsk.',1,N_s));
    [~,S(ik,:)] = min(Distance); 
end

for in = 1:1:N_s
    temp = [];
    for ik = 1:1:K
        temp = [temp dec2bin(S(ik,in)-1,log(omega(ik))/log(2))];
    end
    index(in) = bin2dec(temp) + 1;
end


