clc
clear
close all
% DIM_RATE在文中是mixrate BSA交叉过程中的混合率参数(mixrate)通过使用控制在试验中突变的个体元素的数量
Maxiter=500; % 最大迭代次数
PS=30;
D=30;
Xmin=-100;
Xmax=100;
DIM_RATE=1;
% 对种群实现初始化操作
if numel(Xmin)==1
    Xmin=Xmin*ones(1,D);
    Xmax=Xmax*ones(1,D);
end
pop=GeneratePopulation(PS,D,Xmin,Xmax);
fitnesspop=zeros(1,PS);
fitnessoffsprings=zeros(1,PS);
for i=1:PS
    fitnesspop(i)=fun(pop(i,:));
end
% 适用于第一次迭代的时候，还没有历史的位置
historical_pop=GeneratePopulation(PS,D,Xmin,Xmax); % see Eq.2 in [1]
% historical_pop  is swarm-memory of BSA as mentioned in [1].
% ------------------------------------------------------------------------------------------
fitness_record=zeros(1,Maxiter);

for iter=1:Maxiter
    %SELECTION-I
    if rand<rand
        historical_pop=pop;
    end  % see Eq.3 in [1]
    % randperm(PS) 会产生PS个不相同的随机整数
    historical_pop=historical_pop(randperm(PS),:); % see Eq.4 in [1]
    % F是幅度参数
    F=get_scale_factor; % see Eq.5 in [1], you can other F generation strategies

    % map是一个N*D的矩阵，if是让算法在每一行的ceil(DIM_RATE*rand*D)元素为1，
    % 后续的(ceil(DIM_RATE*rand*D)+1)到D的元素为0
    map=zeros(PS,D); % see Algorithm-2 in [1]
    if rand<rand
        for i=1:PS
            u=randperm(D);
            map(i,u(1:ceil(DIM_RATE*rand*D)))=1;
        end
    else
        for i=1:PS
            map(i,randi(D))=1; % map矩阵中第i行中的D个元素中随机选一个元素为1
        end
    end
    % RECOMBINATION (MUTATION+CROSSOVER)
    % MUTATION
    % 突变操作的实质是对部分维度（由前面的map数组确定）中的元素采取干扰操作，
    offsprings=pop+(map.*F).*(historical_pop-pop);   % see Eq.5 in [1]
    % 对产生的offsprings数组进行边界处理，防止可行解跃出可行边界
    offsprings=BoundaryControl(offsprings,Xmin,Xmax); % see Algorithm-3 in [1]
    % SELECTON-II 按照offsprings和pop的适应值来筛选留下来的个体
    for i=1:PS
        fitnessoffsprings(i)=fun(offsprings(i,:));
    end
    ind=fitnessoffsprings<fitnesspop;
    fitnesspop(ind)=fitnessoffsprings(ind);
    pop(ind,:)=offsprings(ind,:);
    [globalminimum,ind]=min(fitnesspop);
    globalminimizer=pop(ind,:);
    fitness_record(iter)=globalminimum;
end
plot(fitness_record,'r')

%% 部分函数
function pop=GeneratePopulation(popsize,dim,low,up)
pop=ones(popsize,dim);
for i=1:popsize
    for j=1:dim
        pop(i,j)=rand*(up(j)-low(j))+low(j);
    end
end
end
%---------------------------------------------------------------------------------------------
function pop=BoundaryControl(pop,low,up)
[popsize,dim]=size(pop);
for i=1:popsize
    for j=1:dim
        k=rand<rand; % you can change boundary-control strategy
        if pop(i,j)<low(j), if k, pop(i,j)=low(j); else pop(i,j)=rand*(up(j)-low(j))+low(j); end, end
        if pop(i,j)>up(j),  if k, pop(i,j)=up(j);  else pop(i,j)=rand*(up(j)-low(j))+low(j); end, end
    end
end
end


%---------------------------------------------------------------------------------------------
function F=get_scale_factor % you can change generation strategy of scale-factor,F
F=3*randn;                 % STANDARD brownian-walk
% F=4*randg;               % brownian-walk
% F=lognrnd(rand,5*rand);  % brownian-walk
% F=1/normrnd(0,5);        % pseudo-stable walk (levy-like)
% F=1./gamrnd(1,0.5);      % pseudo-stable walk (levy-like, simulates inverse gamma distribution; levy-distiribution)
end

%---------------------------------------------------------------------------------------------
function f=fun(x)
f=sum(x.^2);
end