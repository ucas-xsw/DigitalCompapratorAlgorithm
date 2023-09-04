clc;
clear;
close all;
%global VTH alp_p  fclk ceq gnd ag  gpu beta_m ku kd
kk=0;
BestCost_de_edp_edp=zeros(10,40,500);
BestPosition_de_edp=zeros(10,40,7);
BestPosition_de_edp_itr=zeros(10,40,500,7);
BestVgthre_de_edp=zeros(10,40);
BestVgthre_de_edp_itr=zeros(10,40,500);
t_de_edp=zeros(10,40);
Best_de_edp_edp_itr=zeros(10,40,500);
% for ii=1:5
%     VDD=0.2*ii;
%     for jj=1:5
%     vgthid=0.2*jj;
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
        
        MaxIt = 500;      % Maximum Number of Iterations
        
        nPop = 50;        % Population Size
        
        beta_min = 0.2;   % Lower Bound of Scaling Factor
        beta_max = 0.8;   % Upper Bound of Scaling Factor
        
        pCR = 0.2;        % Crossover Probability
        
        %% Initialization
        
        empty_individual.Position = [];
        empty_individual.Cost = [];
        
        BestSol.Cost = inf;
        
        pop = repmat(empty_individual, nPop, 1);
        
        for i = 1:nPop
            
            pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
            
            pop(i).Cost = CostFunction(pop(i).Position',VDD,vgthid);
            
            if pop(i).Cost<BestSol.Cost
                BestSol = pop(i);
            end
            
        end
        tic
        %BestCost_de_edp_edp = zeros(MaxIt, 1);
        
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
            BestCost_de_edp_edp(ii,jj,it) = BestSol.Cost;
            % Show Iteration Information
            disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost_de_edp_edp(ii,jj,it))] );
            BestPosition_de_edp_itr(ii,jj,it,:)=BestSol.Position;
            % Update Best Cost
            BestCost_de_edp_edp(ii,jj,it) = BestSol.Cost;
            % Show Iteration Information
            %disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost_de_edp_edp(ii,jj,it))] );
            BestPosition_de_edp_itr(ii,jj,it,:)=BestSol.Position;
            BestVgthre_de_edp_itr(ii,jj,it)=vgthre(BestSol.Position',VDD);
            Best_de_edp_edp_itr(ii,jj,it)=getedp(BestSol.Position',VDD);
            %         if BestSol.Position==zeros(1,1,7)
            %             BestSol.Position=BestPosition_de_edp_itr(ii,jj,it-1,:);
            %             break;
            %         end
        end
        t_de_edp(ii,jj)=toc;
        % Show Results
        BestPosition_de_edp(ii,jj,:)=BestSol.Position;
        BestVgthre_de_edp(ii,jj)= vgthre(BestSol.Position',VDD);
        %      plot(reshape(BestCost_de_edp_edp(ii,jj,:),[MaxIt,1]));
        %     semilogy(BestCost_de_edp_edp, 'LineWidth', 2);
        %     xlabel('Iteration');
        %     ylabel('Best Cost');
        %     grid on;
        %     hold on
    end
end

save best_de_edp_edp.mat BestPosition_de_edp BestVgthre_de_edp BestCost_de_edp_edp t_de_edp