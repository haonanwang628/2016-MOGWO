%___________________________________________________________________%
%  Multi-Objective Grey Wolf Optimizer (MOGWO)                      %
%  Source codes demo version 1.0                                    %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper:                                                     %
%                                                                   %
%    S. Mirjalili, S. Saremi, S. M. Mirjalili, L. Coelho,           %
%    Multi-objective grey wolf optimizer: A novel algorithm for     %
%    multi-criterion optimization, Expert Systems with Applications,%
%    in press, DOI: http://dx.doi.org/10.1016/j.eswa.2015.10.039    %       %
%                                                                   %
%___________________________________________________________________%

% I acknowledge that this version of MOGWO has been written using
% a large portion of the following code:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MATLAB Code for                                                  %
%                                                                   %
%  Multi-Objective Particle Swarm Optimization (MOPSO)              %
%  Version 1.0 - Feb. 2011                                          %
%                                                                   %
%  According to:                                                    %
%  Carlos A. Coello Coello et al.,                                  %
%  "Handling Multiple Objectives with Particle Swarm Optimization," %
%  IEEE Transactions on Evolutionary Computation, Vol. 8, No. 3,    %
%  pp. 256-279, June 2004.                                          %
%                                                                   %
%  Developed Using MATLAB R2009b (Version 7.9)                      %
%                                                                   %
%  Programmed By: S. Mostapha Kalami Heris                          %
%                                                                   %
%         e-Mail: sm.kalami@gmail.com                               %
%                 kalami@ee.kntu.ac.ir                              %
%                                                                   %
%       Homepage: http://www.kalami.ir                              %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all
% clc

% drawing_flag = 1;

% TestProblem='UF1';
% nVar=10;
% 
% fobj = cec09(TestProblem);
% 
% xrange = xboundary(TestProblem, nVar);
% 
% % Lower bound and upper bound
% lb=xrange(:,1)';
% ub=xrange(:,2)';
% 
% VarSize=[1 nVar];

%%%% 设置实验参数范围，包括实验次数和实验函数范围
Num_Test=5;   %%%% 每个函数独立进行Num_Test轮?
Num_Experiment=30;   %%%% 函数是从F1-FNum_Functions
AlgorithmName='MOGWO'; %%% 控制函数名?

GreyWolves_num=100;
MaxIt=3000;  % Maximum Number of Iterations
Archive_size=100;   % Repository Size

% alpha=0.1;  % Grid Inflation Parameter
% nGrid=10;   % Number of Grids per each Dimension
% beta=4; %=4;    % Leader Selection Pressure Parameter
% gamma=2;    % Extra (to be deleted) Repository Member Selection Pressure

m=2;
ALLFunction_AllTest=[];

% for ff=[1:21];
% for ff=[10:21];
% for ff=[1:4,6:9];
for ff=[1];
    clearvars -except Num_Test Num_Experiment AlgorithmName ALLFunction_AllTest ff MaxIt GreyWolves_num Archive_size m AllTest_Results  problem_name 
    %%%%% 创建文件夹?
    string_0ALL=['00000\',AlgorithmName,'_5维目标800次迭代100种群实验20210923\'];
    dirname00=[string_0ALL,'\F',num2str(ff),'\'];
    display(['**********  ',AlgorithmName,'算法优化F',num2str(ff),'的 ', 'M',num2str(m), ' 维实验   **********']);
   for testi=1:Num_Test   %%%% 控制每次实验测试次数
       dirname0=[dirname00,'test',num2str(testi),'_F',num2str(ff)];
       system(['mkdir ' dirname0]) %创建主文件夹
       dirname1=[dirname0,'\F',num2str(ff),'_fig'];
       system(['mkdir ' dirname1]) %创建文件夹  等待保存实验图像
       dirname2=[dirname0,'\F',num2str(ff),'_data'];
       system(['mkdir ' dirname2]) %创建文件夹  等待保存实验图像
       for kk=1:30 %%%% 控制实验次数的循环
           display(['**********  ',AlgorithmName,'算法优化F',num2str(ff),'的  第  ', num2str(kk), ' 次实验   **********']);
            rand('state',sum(100*clock));
% Initialization

TestProblem='UF1';
nVar=10;

fobj = cec09(TestProblem);

xrange = xboundary(TestProblem, nVar);

% Lower bound and upper bound
lb=xrange(:,1)';
ub=xrange(:,2)';

VarSize=[1 nVar];

alpha=0.1;  % Grid Inflation Parameter
nGrid=10;   % Number of Grids per each Dimension
beta=4; %=4;    % Leader Selection Pressure Parameter
gamma=2;    % Extra (to be deleted) Repository Member Selection Pressure

GreyWolves=CreateEmptyParticle(GreyWolves_num);


for i=1:GreyWolves_num
    GreyWolves(i).Velocity=0;
    GreyWolves(i).Position=zeros(1,nVar);
    for j=1:nVar
        GreyWolves(i).Position(1,j)=unifrnd(lb(j),ub(j),1);
    end
    GreyWolves(i).Cost=fobj(GreyWolves(i).Position')';
    GreyWolves(i).Best.Position=GreyWolves(i).Position;
    GreyWolves(i).Best.Cost=GreyWolves(i).Cost;
end

GreyWolves=DetermineDomination(GreyWolves);

Archive=GetNonDominatedParticles(GreyWolves);

Archive_costs=GetCosts(Archive);
G=CreateHypercubes(Archive_costs,nGrid,alpha);

for i=1:numel(Archive)
    [Archive(i).GridIndex Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
end

% MOGWO main loop

for it=1:MaxIt
    a=2-it*((2)/MaxIt);
    for i=1:GreyWolves_num
        
        clear rep2
        clear rep3
        
        % Choose the alpha, beta, and delta grey wolves
        Delta=SelectLeader(Archive,beta);
        Beta=SelectLeader(Archive,beta);
        Alpha=SelectLeader(Archive,beta);
        
        % If there are less than three solutions in the least crowded
        % hypercube, the second least crowded hypercube is also found
        % to choose other leaders from.
        if size(Archive,1)>1
            counter=0;
            for newi=1:size(Archive,1)
                if sum(Delta.Position~=Archive(newi).Position)~=0
                    counter=counter+1;
                    rep2(counter,1)=Archive(newi);
                end
            end
            Beta=SelectLeader(rep2,beta);
        end
        
        % This scenario is the same if the second least crowded hypercube
        % has one solution, so the delta leader should be chosen from the
        % third least crowded hypercube.
        if size(Archive,1)>2
            counter=0;
            for newi=1:size(rep2,1)
                if sum(Beta.Position~=rep2(newi).Position)~=0
                    counter=counter+1;
                    rep3(counter,1)=rep2(newi);
                end
            end
            Alpha=SelectLeader(rep3,beta);
        end
        
        % Eq.(3.4) in the paper
        c=2.*rand(1, nVar);
        % Eq.(3.1) in the paper
        D=abs(c.*Delta.Position-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand(1, nVar)-a;
        % Eq.(3.8) in the paper
        X1=Delta.Position-A.*abs(D);
        
        
        % Eq.(3.4) in the paper
        c=2.*rand(1, nVar);
        % Eq.(3.1) in the paper
        D=abs(c.*Beta.Position-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand()-a;
        % Eq.(3.9) in the paper
        X2=Beta.Position-A.*abs(D);
        
        
        % Eq.(3.4) in the paper
        c=2.*rand(1, nVar);
        % Eq.(3.1) in the paper
        D=abs(c.*Alpha.Position-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand()-a;
        % Eq.(3.10) in the paper
        X3=Alpha.Position-A.*abs(D);
        
        % Eq.(3.11) in the paper
        GreyWolves(i).Position=(X1+X2+X3)./3;
        
        % Boundary checking
        GreyWolves(i).Position=min(max(GreyWolves(i).Position,lb),ub);
        
        GreyWolves(i).Cost=fobj(GreyWolves(i).Position')';
    end
    
    GreyWolves=DetermineDomination(GreyWolves);
    non_dominated_wolves=GetNonDominatedParticles(GreyWolves);
    
    Archive=[Archive
        non_dominated_wolves];
    
    Archive=DetermineDomination(Archive);
    Archive=GetNonDominatedParticles(Archive);
    
    for i=1:numel(Archive)
        [Archive(i).GridIndex Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
    end
    
    if numel(Archive)>Archive_size
        EXTRA=numel(Archive)-Archive_size;
        Archive=DeleteFromRep(Archive,EXTRA,gamma);
        
        Archive_costs=GetCosts(Archive);
        G=CreateHypercubes(Archive_costs,nGrid,alpha);
        
    end
    
%     disp(['In iteration ' num2str(it) ': Number of solutions in the archive = ' num2str(numel(Archive))]);
%     save results
    
    % Results
    
    costs=GetCosts(GreyWolves);
    Archive_costs=GetCosts(Archive);
    
%     if drawing_flag==1
%         hold off
%         plot(costs(1,:),costs(2,:),'k.');
%         hold on
%         plot(Archive_costs(1,:),Archive_costs(2,:),'rd');
%         legend('Grey wolves','Non-dominated solutions');
%         drawnow
%     end
    HisPF{it} = Archive_costs';
end

time=toc;
            PF = Archive_costs';
            cg_curve=HisPF; %%% 历史目标函数值
            Time(kk)=time;
            cc=strcat(dirname2,'\',AlgorithmName,'优化次数_',num2str(kk),'.mat');
            result.time=time;
           
%             true_PF=TPF(m,Archive_size, problem_name);
            true_PF=TPF(m,Archive_size, problem_name);

         %%% using the matlab codes for calculating metric values
         

            hv = HV(PF,true_PF);   %超体积?
            gd = GD(PF, true_PF);                       %世代距离

            sp = Spacing(PF, true_PF);                  %空间分布 
            igd = IGD(PF, true_PF);            %反向世代距离

            hvd(kk)=hv;
            gdd(kk)=gd;
            ssp(kk)=sp;
            igdd(kk)=igd;

             save(cc)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       end
       mean_IGD = mean(igdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的IGD平均值 : ', num2str(mean_IGD)]);
       std_IGD=std(igdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的IGD标准差 : ', num2str(std_IGD)]);
       max_IGD=max(igdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的IGD最大值 : ', num2str(max_IGD)]);
       min_IGD=min(igdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的IGD最小值 : ', num2str(min_IGD)]);
       display('******************************** ');
       
       
       mean_GD = mean(gdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的GD平均值 : ', num2str(mean_GD)]);
       std_GD=std(gdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的GD标准差 : ', num2str(std_GD)]);
       max_GD=max(gdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的GD最大值 : ', num2str(max_GD)]);
       min_GD=min(gdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的GD最小值 : ', num2str(min_GD)]);
       display('******************************** ');
      
       
       mean_HV = mean(hvd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的HV平均值 : ', num2str(mean_HV)]);
       std_HV=std(hvd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的HV标准差 : ', num2str(std_HV)]);
       max_HV=max(hvd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的HV最大值 : ', num2str(max_HV)]);
       min_HV=min(hvd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的HV最小值 : ', num2str(min_HV)]);
       display('******************************** ');
       
       mean_SP = mean(ssp);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的SP平均值 : ', num2str(mean_SP)]);
       std_SP=std(ssp);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的SP标准差 : ', num2str(std_SP)]);
       max_SP=max(ssp);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的SP最大值 : ', num2str(max_SP)]);
       min_SP=min(ssp);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的SP最小值 : ', num2str(min_SP)]);
       display('******************************** ');
       
       mean_time=mean(Time);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的运行时间平均值 : ', num2str(mean_time)]);
       std_time=std(Time);
       display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的运行时间标准差 : ', num2str(std_time)]);
       display('******************************** ');
%        mean_X=mean(Best_X);
%         display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的最优解平均值 : ', num2str(mean_X)]);
%         std_X=std(Best_X);
%         display([AlgorithmName,'_Functions_F',num2str(ff),'测试',num2str(kk),'次实验的最优解标准差 : ', num2str(std_X)]);
%         %%%%%%%%%%%%%%%%%%
        cd=strcat(dirname0,'\Result汇总结果.mat');
        Result.IGDmean=mean_IGD;
        Result.IGDstd=std_IGD;
        Result.IGDmax=max_IGD;
        Result.IGDmin=min_IGD;
      
        
        Result.GDmean=mean_GD;
        Result.GDstd=std_GD;
        Result.GDmax=max_GD;
        Result.GDmin=min_GD;
        
        Result.HVmean=mean_HV;
        Result.HVstd=std_HV;
        Result.HVmax=max_HV;
        Result.HVmin=min_HV;
        
        Result.SPmean=mean_SP;
        Result.SPstd=std_SP;
        Result.SPmax=max_SP;
        Result.SPmin=min_SP;
        
        Result.Tmean=mean_time;
        Result.Tstd=std_time;
        %         Result.Xmean=mean_X;
        %         Result.Xstd=std_X;
        %         Result.Best_Y=Best_Y;
        %         Result.Best_X=Best_X;
        Result.Time=Time;
        %         Result.ResultVector=[mean_IGD,std_IGD,max_IGD,min_IGD,mean_GD,std_GD,max_GD,min_GD,mean_time,std_time];
        Result.ResultVector=[mean_IGD,std_IGD,max_IGD,min_IGD,mean_GD,std_GD,max_GD,min_GD,mean_HV,std_HV,max_HV,min_HV,mean_SP,std_SP,max_SP,min_SP,mean_time,std_time];
        %         Result.Best_History_Y=History_Y;
        save(cd,'Result')
        
        
        %         AllTest_Results(testi,:)=[mean_IGD,std_IGD,max_IGD,min_IGD,mean_time,std_time];
        % AllTest_Results(testi,:)=[mean_IGD,std_IGD,max_IGD,min_IGD,mean_GD,std_GD,max_GD,min_GD,mean_time,std_time];
        AllTest_Results(testi,:)=[mean_IGD,std_IGD,mean_GD,std_GD,mean_HV,std_HV,mean_SP,std_SP,mean_time,std_time];
    end
    cd=strcat(dirname00,'Result_AllTest.mat');
    save(cd,'AllTest_Results')
    ALLFunction_AllTest=[ALLFunction_AllTest;AllTest_Results];
end


cd=strcat(string_0ALL,'ALLFunction_AllTest.mat');
save(cd,'ALLFunction_AllTest')


