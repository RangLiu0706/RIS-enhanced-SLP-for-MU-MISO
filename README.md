### About the paper
This is a code package for the paper: 
R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “Joint symbol-level precoding and reflecting designs for IRS-enhanced MU-MISO systems,” IEEE Trans. Wireless Commun., vol. 20, no. 2, pp. 798-811, Feb. 2021.

@ARTICLE{9219206,
  author={Liu, Rang and Li, Ming and Liu, Qian and Swindlehurst, A. Lee},
  journal={IEEE Transactions on Wireless Communications}, 
  title={Joint Symbol-Level Precoding and Reflecting Designs for IRS-Enhanced MU-MISO Systems}, 
  year={2021},
  volume={20},
  number={2},
  pages={798-811},
  keywords={Precoding;Wireless communication;Quality of service;Interference;Phase shift keying;Minimization;Intelligent reflecting surface (IRS);symbol-level precoding (SLP);multiuser multiple-input single-output (MU-MISO) systems;manifold optimization},
  doi={10.1109/TWC.2020.3028371}}


- If you use this simulation code package in any way, please cite the original paper above.
- All codes are contributed by Rang Liu (email: rangl2@uci.edu; website: https://rangliu0706.github.io/). 
   Please feel free to contact with her if you have any suggestions. 
- The link of this paper is: https://ieeexplore.ieee.org/document/9219206
- More information can be found at: https://www.minglabdut.com/resource.html
- Copyright Notice: This code is licensed for personal, non-commercial use only, specifically for academic purposes. Copyright reserved by the MingLab (led by Prof. Ming Li), School of Information and Communication Engineering, Dalian University of Technology, Dalian 116024, China. 


### Software platform
- Please note that the MATLAB2022b is used for this simulation code package, and there may be some imcompatibility problems among different sofrware versions. 
- To run those codes, please download and install [CVX](http://cvxr.com/cvx/) & [Manopt](https://www.manopt.org/)

### Content of this simulation code package
- The files "main_PM_power_iter", "main_PM_power_SNR", "main_PM_power_N", "main_QoS_SER_iter", "main_QoS_SER_power", "main_QoS_SER_N", "main_PM_power_d",  and "main_QoS_SER_d" are used to obtain Figs. 4-13, respectively.


Abstract of the paper: 
Intelligent reflecting surfaces (IRSs) have emerged as a revolutionary solution to enhance wireless communications by changing propagation environment in a cost-effective and hardware-efficient fashion. In addition, symbol-level precoding (SLP) has attracted considerable attention recently due to its advantages in converting multiuser interference (MUI) into useful signal energy. Therefore, it is of interest to investigate the employment of IRS in symbol-level precoding systems to exploit MUI in a more effective way by manipulating the multiuser channels. In this article, we focus on joint symbol-level precoding and reflecting designs in IRS-enhanced multiuser multiple-input single-output (MU-MISO) systems. Both power minimization and quality-of-service (QoS) balancing problems are considered. In order to solve the joint optimization problems, we develop an efficient iterative algorithm to decompose them into separate symbol-level precoding and block-level reflecting design problems. An efficient gradient-projection-based algorithm is utilized to design the symbol-level precoding and a Riemannian conjugate gradient (RCG)-based algorithm is employed to solve the reflecting design problem. Simulation results demonstrate the significant performance improvement introduced by the IRS and illustrate the effectiveness of our proposed algorithms.
