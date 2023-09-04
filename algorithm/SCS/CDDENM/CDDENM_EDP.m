%% 
% CDDENM

clc;
clear;
close all;
%global VTH alp_p  fclk ceq gnd ag  gpu beta_m ku kd
global VDD vgthid ceq gnd gpu
kk=0;
BestCost_denm_edp=zeros(10,40,500);
BestPosition_denm_edp_itr=zeros(10,40,500,7);
BestPosition_denm_edp=zeros(10,40,7);
BestVgthre_denm_edp=zeros(10,40);
BestVgthre_denm_edp_itr=zeros(10,40,500);
t_denm_edp=zeros(0,40);
Best_denm_edp_itr=zeros(10,40,500);
% BestCost_nm=zeros(10,40,500);
% BestPosition_nm=zeros(10,40,7);
% BestVgthre_nm=zeros(10,40);
% for ii=1:6
%     VDD=0.2*ii;
%     for jj=1:5
%     vgthid=0.2*jj;
% CoreNum=10; %设定机器CPU核心数量
% if isempty(gcp('nocreate')) %如果并行未开启
%     parpool(CoreNum);
% end
% delete(gcp('nocreate'))

for ii=2:10
    VDD=0.1*ii;
    for jj=1:40
        vgthid=0.01*jj+0.2;
        kk=kk+1
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
        %% Problem Definition
        figure;
        
        CostFunction = @(x,VDD, vgthid) sqrtdeltavth(x,VDD, vgthid);    % Cost Function
        
        nVar = 7;            % Number of Decision Variables
        
        VarSize = [1 nVar];   % Decision Variables Matrix Size
        
        VarMin = 0;          % Lower Bound of Decision Variables
        VarMax = 10;          % Upper Bound of Decision Variables
        
        %% DE Parameters
        
        MaxIt = 300;      % Maximum Number of Iterations
        
        nPop = 50;        % Population Size
        
        beta_min = 0.2;   % Lower Bound of Scaling Factor
        beta_max = 0.8;   % Upper Bound of Scaling Factor
        
        pCR = 0.2;        % Crossover Probability
        
        %% Initialization
        empty_individual=struct();
        empty_individual.Position = [];
        empty_individual.Cost = [];
        
        BestSol.Cost = inf;
        initial_point = ones(7, 1) * 7;
        initial_simplex = nelders_favorite_initial_simplex(initial_point,[1, 2, 3, 4,5,6,7]);
        func = @sqrtdeltavth1;
        
        options = xoptimset(                                                           ...
            'Display', 'iter',                                                         ...
            'InitialSimplex', initial_simplex                                          ...
            );
        tic
        [x1, fval, exitflag, output, Best_cost] = fminsearch_nm_bd_dis(func, initial_point, options);
        BestSol.Cost=fval;
        BestSol.Position=x1';
        
        % 定义基础向量
        %x = [1 2 3 4 5];
        if BestSol.Cost< 1e-8
            BestCost_denm_edp(ii,jj,:) = Best_cost;
            BestPosition_denm_edp(ii,jj,:)=BestSol.Position;
            BestVgthre_denm_edp(ii,jj)= vgthre(BestSol.Position',VDD);
            continue% 跳出之前要存储
        end
        
        % 定义均值和方差
        mu = 0;
        sigma = 1;
        
        pop = repmat(empty_individual, nPop, 1);
        
        rng('shuffle');
        
        for i = 1:nPop
                        % 按照正态分布生成随机数
            rand_vec = round(normrnd(mu, sigma, 1, nVar));
            
            % 将负数转换为0
            rand_vec(rand_vec < 0) = 0;
            
            % 生成新向量
            new_vec = x1' + rand_vec;
            
            % 打印新向量
            fprintf('New vector %d: ', i);
            disp(new_vec);
            
            pop(i).Position = new_vec;
            
            pop(i).Cost = CostFunction(pop(i).Position',VDD,vgthid);
            
            if pop(i).Cost<BestSol.Cost
                BestSol = pop(i);
            end
            
        end
        
        %BestCost = zeros(MaxIt, 1);
        
        %% DE Main Loop
        
        for it = 1:MaxIt
            
            for i = 1:nPop
                
                x = pop(i).Position;
                
                A = randperm(nPop);
                
                A(A == i) = [];
                
                a = A(1);
                b = A(2);
                c = A(3);
                
                % Mutation
                %beta = unifrnd(beta_min, beta_max);
                beta = unifrnd(beta_min, beta_max, VarSize);
                y = pop(a).Position+beta.*(pop(b).Position-pop(c).Position);
                y = max(y, VarMin);
                y = min(y, VarMax);
                
                % Crossover
                z = zeros(size(x));
                j0 = randi([1 numel(x)]);
                for j = 1:numel(x)
                    if j == j0 || rand <= pCR
                        z(j) = y(j);
                    else
                        z(j) = x(j);
                    end
                end
                z=round(z);
                NewSol=struct();
                NewSol.Position = z;
                NewSol.Cost = CostFunction(NewSol.Position',VDD,vgthid);
                
                if NewSol.Cost<pop(i).Cost
                    pop(i) = NewSol;
                    
                    if pop(i).Cost<BestSol.Cost
                        BestSol = pop(i);
                    end
                end
                
            end
            
            % Update Best Cost
            BestCost_denm_edp(ii,jj,it) = BestSol.Cost;
            % Show Iteration Information
            %disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost_denm_edp(ii,jj,it))] );
            BestPosition_denm_edp_itr(ii,jj,it,:)=BestSol.Position;
            BestVgthre_denm_edp_itr(ii,jj,it)=vgthre(BestSol.Position',VDD);
            Best_denm_edp_itr(ii,jj,it)=getedp(BestSol.Position',VDD);
%             if BestSol.Position==zeros(1,1,7)
%                 BestSol.Position=BestPosition_denm_edp_itr(ii,jj,it-1,:);
%                 break;
%             end
        end
        t_denm_edp(ii,jj)=toc;
        %% Show Results
        BestPosition_denm_edp(ii,jj,:)=BestSol.Position;
        BestVgthre_denm_edp(ii,jj)= vgthre(BestSol.Position',VDD);
%         BestCost_nm(ii,jj,:)=Best_cost(1:500);
%         BestPosition_nm(ii,jj,:)=x1;
%         BestVgthre_nm(ii,jj)= vgthre(x1,VDD)
        %      plot(reshape(BestCost(ii,jj,:),[MaxIt,1]));
        %     semilogy(BestCost, 'LineWidth', 2);
        %     xlabel('Iteration');
        %     ylabel('Best Cost');
        %     grid on;
        %     hold on
        
    end
end

save best_denm_edp.mat Best_denm_edp_itr BestPosition_denm_edp BestVgthre_denm_edp BestVgthre_denm_edp_itr BestCost_denm_edp t_denm_edp BestPosition_denm_edp_itr