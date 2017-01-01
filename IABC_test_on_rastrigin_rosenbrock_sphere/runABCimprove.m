%%%%%ARTIFICIAL BEE COLONY ALGORITHM%%%%

%Artificial Bee Colony Algorithm was developed by Dervis Karaboga in 2005 
%by simulating the foraging behaviour of bees.

%Copyright ?2008 Erciyes University, Intelligent Systems Research Group, The Dept. of Computer Engineering

%Contact:
%Dervis Karaboga (karaboga@erciyes.edu.tr )
%Bahriye Basturk Akay (bahriye@erciyes.edu.tr)

function GlobalMins=runABCimprove(b)

%******************928加入了Employed_Best和Employed_SecondBest************
% clear all
% close all
% clc



% Set ABC Control Parameters
ABCOpts = struct( 'ColonySize',  20, ...   % Number of Employed Bees+ Number of Onlooker Bees 
    'MaxCycles', 2000,...   % Maximum cycle number in order to terminate the algorithm
    'ErrGoal',     1e-20, ...  % Error goal in order to terminate the algorithm (not used in the code in current version)
    'Dim',       2 , ... % Number of parameters of the objective function   
    'Limit',   100, ... % Control paramter in order to abandone the food source 
    'lb',  -3, ... % Lower bound of the parameters to be optimized
    'ub',  3, ... %Upper bound of the parameters to be optimized
    'ObjFun' , 'griewank', ... %Write the name of the objective function you want to minimize
    'RunTime',1); % Number of the runs 



GlobalMins=zeros(ABCOpts.RunTime,ABCOpts.MaxCycles);

for r=1:ABCOpts.RunTime
    
% Initialise population
Range = repmat((ABCOpts.ub-ABCOpts.lb),[ABCOpts.ColonySize ABCOpts.Dim]);
Lower = repmat(ABCOpts.lb, [ABCOpts.ColonySize ABCOpts.Dim]);
Colony = rand(ABCOpts.ColonySize,ABCOpts.Dim) .* Range + Lower;%生成初始Colony，其中ColonySize行，Dim列，10*5
%zj先初始化种群规模。。。这个就是算法中式子：x(j)i=x(j)min+rand(0,1)(x(j)max-x(j)min)


ObjEmp=feval(ABCOpts.ObjFun,Colony);
FitEmp=calculateFitness(ObjEmp);
FitEmpArray=[FitEmp';1:length(Colony)]';%928除了Sphere，其他测试函数在FitEmp都要转置即FitEmp'
FitEmpSort=sortrows(FitEmpArray);
Employed=Colony(FitEmpSort(length(Colony)/2+1:length(Colony),2),:);
FitEmp=FitEmpSort(length(Colony)/2+1:length(Colony),1)';


%Employed=Colony(1:(ABCOpts.ColonySize/2),:);%前一半为引领蜂或食物源，5*5
%zj再将种群的前一半作为引领蜂规模

         %***********全局***********
        FitArrayCounterNumber=[FitEmp;1:length(Colony)/2]';%928除了Sphere，其他测试函数在FitEmp都要转置即FitEmp'
        maxFitArray=sortrows(FitArrayCounterNumber);
        Employed_Best=Employed(maxFitArray(length(maxFitArray),2),:);
        Employed_SecondBest=Employed(maxFitArray(length(maxFitArray)-1,2),:);

%evaluate and calculate fitness
ObjEmp=feval(ABCOpts.ObjFun,Employed);%计算引领蜂组Employed的每一行（每一个食物源Xi）的目标函数值，1*5
FitEmp=calculateFitness(ObjEmp);%计算食物源的适应度值，1*5

%set initial values of Bas
Bas=zeros(1,(ABCOpts.ColonySize/2));%1*5，每一列对应每一个Xi没有被更新的次数


GlobalMin=ObjEmp(find(ObjEmp==min(ObjEmp),end));%不带end结果一样，因为ObjEmp为一个列向量，找出目标函数ObjEmp中最小的值给了GlobalMin，1*1
GlobalParams=Employed(find(ObjEmp==min(ObjEmp),end),:);%GlobalParams存放目标函数ObjEmp最小时对应的解,1*5

Cycle=1;
while ((Cycle <= ABCOpts.MaxCycles)),%开始循环
    
    %%%%% Employed phase
    Employed2=Employed;%Employed2的每一行对应于Employed每一行的邻域搜索值，即Employed每一行代表Xij，Employed2每一行代表Vij
    for i=1:ABCOpts.ColonySize/2%对Xij的每一个i，只有一个j改变。
        Param2Change=fix(rand*ABCOpts.Dim)+1;%对应于j
        %zj这就是算法中要找的xij中的j，其是Dim中的一个随机数
        neighbour=fix(rand*(ABCOpts.ColonySize/2))+1;%对应于k
        neighbour1=fix(rand*(ABCOpts.ColonySize/2))+1;%对应于k
           
            while(neighbour==i||neighbour1==i||neighbour==neighbour1)
                neighbour=fix(rand*(length(Colony)/2))+1;
                neighbour1=fix(rand*(length(Colony)/2))+1;
            end;
        Employed2(i,Param2Change)=Employed_Best(1,Param2Change)+(Employed(neighbour,Param2Change)-Employed(neighbour1,Param2Change))*(rand-0.5)*2+(Employed_SecondBest(1,Param2Change)-Employed(i,Param2Change))*(rand-0.5)*2;
        %Employed2(i,Param2Change)=Employed(i,Param2Change)+(Employed(neighbour,Param2Change)-Employed(neighbour1,Param2Change))*(rand-0.5)*2+(Employed_Best(1,Param2Change)-Employed(i,Param2Change))*(rand-0.5)*2;
        %上面有两种位置更新公式，经测试发现似乎上面的带有Employed_Best和Employed_SecondBest要比下面的位置更新
        %公式在收敛时速度更快。
        
        
        
        %Employed2(i,Param2Change)=Employed(i,Param2Change)+(Employed(i,Param2Change)-Employed(neighbour,Param2Change))*(rand-0.5)*2;
        %不能超过上下界
       %zj这个就是算法中式子：vij=xij+suijishu(xij-xkj),上式中(rand-0.5)*2就是为了限制suijishu的范
       %围控制在（-1,1）之间的
        if (Employed2(i,Param2Change)<ABCOpts.lb)
             Employed2(i,Param2Change)=ABCOpts.lb;
         end;
        if (Employed2(i,Param2Change)>ABCOpts.ub)
            Employed2(i,Param2Change)=ABCOpts.ub;
        end;
         %zj可能出于规范的考虑，前面设置了参数lb和ub就是用来限定新的位置不应该越界，如果比下限小，则将下限赋给它，如果比上限大
        %则将上线赋给它
end;   

    ObjEmp2=feval(ABCOpts.ObjFun,Employed2);%计算每一个Vi的目标函数值，5*1
    FitEmp2=calculateFitness(ObjEmp2);%计算每一个Vi的适应度值，5*1
    [Employed ObjEmp FitEmp Bas]=GreedySelection(Employed,Employed2,ObjEmp,ObjEmp2,FitEmp,FitEmp2,Bas,ABCOpts);
     %zj  Employed：之前的位置数组；   Employed2:更新位置后的数组；   
    %ObjEmp：之前的目标函数值（每只蜜蜂）；    ObjEmp2：更新位置后的目标函数值（即每只蜜蜂的函数值）
    %Bas：表示被开采的次数；   ABCOpts：整个结构体zj
    
    %贪婪原则选择，Employed ObjEmp FitEmp分别对应选择后的食物源，函数值和适应度
    %Normalize
    NormFit=FitEmp/sum(FitEmp);%50*1，每一行对应于Xi和Vi较优解的Pi
    
     %***********全局***********
        FitArrayCounterNumber=[FitEmp';1:length(Colony)/2]';%928除了Sphere，其他测试函数在FitEmp都要转置即FitEmp'
        maxFitArray=sortrows(FitArrayCounterNumber);
        Employed_Best=Employed(maxFitArray(length(maxFitArray),2),:);
        Employed_SecondBest=Employed(maxFitArray(length(maxFitArray)-1,2),:);
    
    %%% Onlooker phase  
Employed2=Employed;
%zj本来这句是Employed2=Employed；但是通过本人阅读论文来看，一般情况是总样本数=引领蜂数量+跟随蜂数量；引领蜂数量=跟随蜂数量；并
%且，引领蜂是前50%，跟随蜂是后50%zj
i=1;
t=0;
while(t<ABCOpts.ColonySize/2) 
    %轮盘赌
    if(rand<NormFit(i))%NormFit(i)越大，被选中的次数越多
        t=t+1;%所有跟随蜂都必须选择引领蜂进行跟随
        Param2Change=fix(rand*ABCOpts.Dim)+1;
        neighbour=fix(rand*(ABCOpts.ColonySize/2))+1;
        
            while(neighbour==i||neighbour1==i||neighbour==neighbour1)
                neighbour=fix(rand*(length(Colony)/2))+1;
                neighbour1=fix(rand*(length(Colony)/2))+1;
            end;
        Employed2(i,Param2Change)=Employed_Best(1,Param2Change)+(Employed(neighbour,Param2Change)-Employed(neighbour1,Param2Change))*(rand-0.5)*2+(Employed_SecondBest(1,Param2Change)-Employed(i,Param2Change))*(rand-0.5)*2;
       % Employed2(i,Param2Change)=Employed(i,Param2Change)+(Employed(neighbour,Param2Change)-Employed(neighbour1,Param2Change))*(rand-0.5)*2+(Employed_Best(1,Param2Change)-Employed(i,Param2Change))*(rand-0.5)*2; 
        %上面有两种位置更新公式，经测试发现似乎上面的带有Employed_Best和Employed_SecondBest要比下面的位置更新
        %公式在收敛时速度更快。
        
        
        %Employed2(i,Param2Change)=Employed(i,Param2Change)+(Employed(i,Param2Change)-Employed(neighbour,Param2Change))*(rand-0.5)*2;
         if (Employed2(i,Param2Change)<ABCOpts.lb)
             Employed2(i,Param2Change)=ABCOpts.lb;
         end;
        if (Employed2(i,Param2Change)>ABCOpts.ub)
            Employed2(i,Param2Change)=ABCOpts.ub;
         end;
    ObjEmp2=feval(ABCOpts.ObjFun,Employed2);%计算跟随蜂邻域搜索解的目标函数
    FitEmp2=calculateFitness(ObjEmp2);%计算跟随蜂邻域搜索解的适应度
    [Employed ObjEmp FitEmp Bas]=GreedySelection(Employed,Employed2,ObjEmp,ObjEmp2,FitEmp,FitEmp2,Bas,ABCOpts,i);
   
   end;
    
    i=i+1;
    if (i==(ABCOpts.ColonySize/2)+1)  %如果超出范围，将i至1
        i=1;
    end;   
end;
    
    
    %%%Memorize Best
 CycleBestIndex=find(FitEmp==max(FitEmp));
 CycleBestIndex=CycleBestIndex(end);%我认为可以不要
 CycleBestParams=Employed(CycleBestIndex,:);%原本注释是“求每次循环中适度值最小所对应的食物源（解）”我认为是“求每次循环中适度值最大所对应的食物源（解）”
 CycleMin=ObjEmp(CycleBestIndex);%原本这句注释是“求每次循环中适度值最小所对应的函数值”，我认为是“求每次循环中适度值最大所对应的函数值”
 
 if CycleMin<GlobalMin %和全局最小进行比较
       GlobalMin=CycleMin;
       GlobalParams=CycleBestParams;
 end
 
 GlobalMins(r,Cycle)=GlobalMin;%记录每次循环所对应的全局最小值，1*2000
 
 %% Scout phase
 ind=find(Bas==max(Bas));%找到没有被更新次数最多的那个食物源Xi，并把次数和limit比较
ind=ind(end);
if (Bas(ind)>ABCOpts.Limit)
Bas(ind)=0;
%Employed(ind,:)=(ABCOpts.ub-ABCOpts.lb)*(0.5-rand(1,ABCOpts.Dim))*2+ABCOpt
%s.lb;
Employed(ind,:)=(ABCOpts.ub-ABCOpts.lb)*(0.5-rand(1,ABCOpts.Dim))*2;%+ABCOpts.lb;
%        Employed2(i,Param2Change)=Employed(i,Param2Change)+(Employed(i,Par
%        am2Change)-Employed(neighbour,Param2Change))*(rand-0.5)*2;
%message=strcat('burada',num2str(ind))
end;
ObjEmp=feval(ABCOpts.ObjFun,Employed);
FitEmp=calculateFitness(ObjEmp);
    


    fprintf('Cycle=%d ObjVal=%g\n',Cycle,GlobalMin);
    
    Cycle=Cycle+1;

end % End of ABC

end; %end of runs
% if ABCOpts.RunTime==1
%     semilogy(GlobalMins,'b');
% else
%     semilogy(mean(GlobalMins),'b');%若多次执行，求均值
% end
% %semilogy(mean(GlobalMins))
% title('Mean of Best function values');
% xlabel('cycles');
% ylabel('error');
fprintf('Mean =%g Std=%g\n',mean(GlobalMins(:,end)),std(GlobalMins(:,end)));%输出GlobalMins的均值和方差
  