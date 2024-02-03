% Calculate the SER.
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, "Joint symbol-level precoding and reflecting designs for IRS-enhanced MU-MISO systems,” IEEE Trans. Wireless Commun., vol. 20, no. 2, pp. 798-811, Feb. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9219206
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02

function ser = get_SER(R,S,omega)

[K,N_s] = size(R);
ser = zeros(1,K);
for ik = 1:1:K
    mpsk = exp( 1i*(pi/omega(ik)+2*pi/omega(ik)*(0:1:omega(ik)-1)) );
    Distance = abs(repmat(R(ik,:),omega(ik),1)-repmat(mpsk.',1,N_s));
    [~,R(ik,:)] = min(Distance); 
    [~,ser(ik)] = symerr(R(ik,:),S(ik,:));
end




