% This examples presents one of Nelder's favorite examples, as described in [1]
% (see Excercise 8.6.1).
%
% References:
%   [1] Kelley, Carl T. Iterative methods for optimization. Society for
%       Industrial and Applied Mathematics, 1999.
clc
close all
clear
global VDD VTH alp_p  fclk ceq gnd ag  gpu beta_m ku kd vgthid
VDD=1;
vgthid=0.5;
t_nm_edp=zeros(10,40);
BestCost_nm_edp=zeros(10,40,500);
BestPosition_nm_edp=zeros(10,40,7);
BestVgthre_nm_edp=zeros(10,40);
Best_denm_edp_itr=zeros(10,40,500);
BestVgthre_denm_edp_itr=zeros(10,40,500);
BestPosition_denm_edp_itr=zeros(10,40,500,7);

for ii=2:10
    VDD=0.1*ii;
    for jj=1:40
        vgthid=0.01*jj+0.2;
        alp_p=0.5;
        fclk=500e6;
        beta_m=1.4;
        m=3;
        H=[1.40E+00 1.40E+00 1.40E+00 1.40E+00 1.40E+00 1.40E+00 1.40E+00]'*1e-6;
        W=[8.00E-01 1.00E+00 1.40E+00 8.00E-01 1.00E+00 1.40E+00 6.00E-01]'*1e-6;
        trs=[2.86E-02 4.96E-02 7.45E-02 1.71E-02 1.96E-02 2.20E-02 1.39E-02]'*1e-9;
        tfs=[1.44E-02 1.66E-02 1.80E-02 2.33E-02 3.75E-02 5.47E-02 1.24E-02]'*1e-9;
        Ih=[1.79E+02 2.69E+02 3.59E+02 5.09E+01 3.57E+01 2.80E+01 8.98E+01]'*1e-6;
        Il=[2.08E+01 3.12E+01 4.16E+01 2.71E+00 1.50E+00 1.06E+00 1.04E+01]'*1e-6;
        tp=(trs+tfs)/2;
        req=VDD./(Il+Ih);
        
        geq=(Il+Ih)/VDD;
        ag=H.*W;
        cin=tp./req;
        geqp=cin./trs;
        geqn=cin./tfs;
        ceq=cin;
        gnd=geqn;
        gpu=geqp;
        kd=[0.5	0.33	0.25	4	3	2	1]';
        ku=[2	3	4	0.25	0.33	0.5	1]';
        
        VTH=0.335;
        etamax=0.1;
        taomax=0.2e-9;
        eta1=((VDD-VTH)/(sqrt(etamax)+vgthid-VTH))^2*kd'-beta_m*ku';
        eta2=-((VDD-VTH)/(-sqrt(etamax)+vgthid-VTH))^2*kd'+beta_m*ku';
        tao1=0.345*ceq'-taomax*geqp';
        tao2=0.345*ceq'-taomax*geqn';
        
        
        initial_point = ones(7, 1) * 7;
        initial_simplex = nelders_favorite_initial_simplex(initial_point,[1, 2, 3, 4,5,6,7]);
        func = @sqrtdeltavth1;
 
        options = xoptimset(                                                           ...
            'Display', 'iter',                                                         ...
            'InitialSimplex', initial_simplex                                          ...
            );
        tic
        [x, fval, exitflag, output, Best_cost] = fminsearch_nm_bd_dis(func, initial_point, options);
        t_nm_edp(ii,jj)=toc;
        BestCost_nm_edp(ii,jj,:)=Best_cost(1:500);
        BestPosition_nm_edp(ii,jj,:)=x;
        BestVgthre_nm_edp(ii,jj)= vgthre(x,VDD);
    end
end

save best_nm_edp.mat BestPosition_nm_edp BestVgthre_nm_edp BestCost_nm_edp t_nm_edp