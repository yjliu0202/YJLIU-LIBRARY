load('origineffmonth') %load efficiency

k=0; 
tic
 
for g=1:180
data=cat(2,G(g).DMU.D2); 

corbootmean=cat(2,G(g).DMU.cor); %读取修正后的bootstrap效率

sum0=0;
for i=1:20
    for j=1:20
             for b=1:500
                 sum0=sum0+(data(b,i)-corbootmean(1,i))*(data(b,j)-corbootmean(1,j));
             end
      cov2(i,j)=sum0/500;  
      sum0=0;
   end
end

bootcov=cov2; %计算效率的协方差

lj=0.05;uj=0.25;U=10; %边界和股票数量参数设置

x = sdpvar(1,20);
x2 = sdpvar(1,1);
P = binvar(1,20);  

%设置约束条件
    set=[];
    for j=1:20
        set=[set;x(j)>=P(j)*lj]; %下边界约束
    end
    for j=1:20
        set=[set;x(j)<=P(j)*uj]; %上边界约束
    end
    sum1=0;
    for i=1:20
         sum1=sum1+x(i)*corbootmean(1,i);  %组合效率
    end
    sum2=0;
    for i=1:20
        for j=1:20
          sum2=sum2+cov2(i,j)*x(i)*x(j);  %组合方差
        end
    end
    sum4=0;
    for j=1:20
        sum4=sum4+x(j);  
    end
    sum5=0;
    for j=1:20
        sum5=sum5+P(j);
    end 
set=[set,sum4==1];   %投资比例加和等于1
set=[set,sum5<=U];   %基数约束

ops = sdpsettings('solver', 'Gurobi+', 'verbose', 2, 'debug', 1);
ops.gurobi.NonConvex = 2;

%求最大化收益
g1=sum1;
optimize(set,-g1,ops);
optimalweight_maxret(g,:)=double(x);%输出最优权重
optimalweight_maxret_P(g,:)=double(P); %输出最优0-1变量
sum_optimalweight_maxret_P(g)=sum(optimalweight_maxret_P(g,:));
optimalweight_maxret_fit(g)=double(g1);% 输出目标函数值

%求最小化收益
g2=sum1;
optimize(set,g2,ops);
optimalweight_minret(g,:)=double(x); %输出最优权重
optimalweight_minret_P(g,:)=double(P);%输出最优0-1变量
sum_optimalweight_minret_P(g)=sum(optimalweight_minret_P(g,:));
optimalweight_minret_fit(g)=double(g2);% 输出目标函数值

%求最大化风险
g3=sum2;
optimize(set,-g3,ops);
optimalweight_maxrisk(g,:)=double(x);%输出最优权重
optimalweight_maxrisk_P(g,:)=double(P);%输出最优0-1变量
sum_optimalweight_maxrisk_P(g)=sum(optimalweight_maxrisk_P(g,:));
optimalweight_maxrisk_fit(g)=double(g3);% 输出目标函数值

g4=sum2;
optimize(set,g4,ops);
optimalweight_minrisk(g,:)=double(x);%输出最优权重
optimalweight_minrisk_P(g,:)=double(P);%输出最优0-1变量
sum_optimalweight_minrisk_P(g)=sum(optimalweight_minrisk_P(g,:));
optimalweight_minrisk_fit(g)=double(g4);% 输出目标函数值

minrisk(g)=optimalweight_minrisk_fit(g);
minret(g)=optimalweight_minret_fit(g);
maxrisk(g)=optimalweight_maxrisk_fit(g);
maxret(g)=optimalweight_maxret_fit(g);

set=[set (sum1-minret(g))/(maxret(g)-minret(g))>=x2(1)];
set=[set (maxrisk(g)-sum2)/(maxrisk(g)-minrisk(g))>=x2(1)];

%求多目标
g5=x2(1);
optimize(set,-g5,ops); 
optimalweight(g,:)=double(x);%输出最优权重
P=double(P);%输出最优0-1变量
HN(g,:)=P;
sumHN(g)=sum(HN(g,:));  % 输出最优股票数量
optimalfit(g)=double(g5);% 输出目标函数值

%% 下面是OR模型求解

p = sdpvar(20,1); %p表示Pj
y = binvar(20,1); % y表示yj

% 添加约束条件
C=[];

    sum1=0;
    for j=1:20
        sum1=sum1+p(j);
    end
    C=[C sum1==1];  %(14b)式

    sum2=0;
    for j=1:20
        sum2=sum2+p(j)*corbootmean(1,j);
    end
    C=[C sum2>=mean(corbootmean)];  %(14c)式

    for j=1:20  
      C=[C p(j)>=lj*y(j)];   %(16)式
    end
 
    for j=1:20
      C=[C p(j)<=y(j)*uj];   %(16)式
    end

    sum3=0;
    for j=1:20
        sum3=sum3+y(j);
    end

    C=[C sum3==U];   %(16)式

ops = sdpsettings('solver', 'Gurobi+', 'verbose', 2, 'debug', 1);
ops.gurobi.NonConvex = 2;

    % 目标函数
    sum25=0;
    for i=1:20
        for j=1:20
          sum25=sum25+bootcov(i,j)*p(i)*p(j); %计算组合方差
        end
    end

z = sum25;

% 求解
reuslt = optimize(C,z,ops);
if reuslt.problem == 0 % problem =0 代表求解成功
    H1=value(p) ;  % 最优投资比例1/∑zi π=x(1)=1/∑zi
    H2=value(y) ;  % 最优0-1决策变量
    H3=value(z) ;  % 最优目标函数值
else
    disp('求解出错');
end

 optimalweight2(g,:)=H1;
 HN2(g,:)=H2;
 sumHN2(g)=sum(HN2(g,:));
 optimalfit2(g,:)=H3;

k=k+1
end

%% 下面是计算投资策略绩效
ori=xlsread('20只能源股收益20060101_20211231month.xlsx',2,'B1:U3649'); %样本外全市场
%% 计算组合收益
% 第一个期限 
    for g=1:20
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(1,j)*ori(g,j);
            sum3=sum3+optimalweight(1,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
    end
    %计算累计财富
    wealth1=1;wealth2=1;wealth3=1;
    for g=1:20
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end

%第二个期限 
   for g=21:35
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(2,j)*ori(g,j);
            sum3=sum3+optimalweight(2,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=21:35
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
 % 第三个期限    
     for g=36:57
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(3,j)*ori(g,j);
            sum3=sum3+optimalweight(3,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=36:57
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
    % 第四个期限   
     for g=58:78
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(4,j)*ori(g,j);
            sum3=sum3+optimalweight(4,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=58:78
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
        % 第五个期限    
     for g=79:96
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(5,j)*ori(g,j);
            sum3=sum3+optimalweight(5,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=79:96
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
       % 第六个期限    
     for g=97:117
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(6,j)*ori(g,j);
            sum3=sum3+optimalweight(6,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=97:117
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
    % 第七个期限  
    for g=118:139
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(7,j)*ori(g,j);
            sum3=sum3+optimalweight(7,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=118:139
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end
    
     % 第八个期限  
    for g=140:162
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(8,j)*ori(g,j);
            sum3=sum3+optimalweight(8,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=140:162
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end
    
      % 第九个期限 
    for g=163:182
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(9,j)*ori(g,j);
            sum3=sum3+optimalweight(9,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=163:182
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end
    
      % 第十个期限  
    for g=183:200
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(10,j)*ori(g,j);
            sum3=sum3+optimalweight(10,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=183:200
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end
    
   % 第十一个期限  
    for g=201:222
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(11,j)*ori(g,j);
            sum3=sum3+optimalweight(11,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=201:222
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end
    
       % 第十二个期限 
    for g=223:242
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(12,j)*ori(g,j);
            sum3=sum3+optimalweight(12,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=223:242
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end
    
      % 第十三个期限 
    for g=243:264
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(13,j)*ori(g,j);
            sum3=sum3+optimalweight(13,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=243:264
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
    
         % 第十四个期限 
    for g=265:280
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(14,j)*ori(g,j);
            sum3=sum3+optimalweight(14,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=265:280
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
    
     % 第十五个期限  
    for g=281:301
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(15,j)*ori(g,j);
            sum3=sum3+optimalweight(15,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=281:301
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
    % 第十六个期限  
    for g=302:322
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(16,j)*ori(g,j);
            sum3=sum3+optimalweight(16,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=302:322
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end    
    
     % 第十七个期限  
    for g=323:342
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(17,j)*ori(g,j);
            sum3=sum3+optimalweight(17,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=323:342
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end    
    
        % 第十八个期限  
    for g=343:362
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(18,j)*ori(g,j);
            sum3=sum3+optimalweight(18,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=343:362
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     
   
   % 第十九个期限  
    for g=363:385
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(19,j)*ori(g,j);
            sum3=sum3+optimalweight(19,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=363:385
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end       
    
     % 第20个期限  
    for g=386:406
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(20,j)*ori(g,j);
            sum3=sum3+optimalweight(20,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=386:406
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
    
     % 第21个期限  
    for g=407:425
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(21,j)*ori(g,j);
            sum3=sum3+optimalweight(21,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=407:425
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end    
    
    % 第22个期限  
    for g=426:445
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(22,j)*ori(g,j);
            sum3=sum3+optimalweight(22,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=426:445
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end    
    
     % 第23个期限  
    for g=446:465
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(23,j)*ori(g,j);
            sum3=sum3+optimalweight(23,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=446:465
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end    
    
  
  % 第24个期限  
    for g=466:488
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(24,j)*ori(g,j);
            sum3=sum3+optimalweight(24,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=466:488
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
    
     % 第25个期限  
    for g=489:503
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(25,j)*ori(g,j);
            sum3=sum3+optimalweight(25,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=489:503
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
    
    % 第26个期限 
    for g=504:523
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(26,j)*ori(g,j);
            sum3=sum3+optimalweight(26,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=504:523
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
      
    % 第27个期限 
    for g=524:545
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(27,j)*ori(g,j);
            sum3=sum3+optimalweight(27,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=524:545
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     
    
    % 第28个期限  
    for g=546:566
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(28,j)*ori(g,j);
            sum3=sum3+optimalweight(28,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=546:566
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     
    
     % 第29个期限  
    for g=567:584
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(29,j)*ori(g,j);
            sum3=sum3+optimalweight(29,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=567:584
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     
    
       % 第30个期限 
    for g=585:606
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(30,j)*ori(g,j);
            sum3=sum3+optimalweight(30,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=585:606
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     
    
        % 第31个期限  
    for g=607:629
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(31,j)*ori(g,j);
            sum3=sum3+optimalweight(31,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=607:629
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     
    
         % 第32个期限  
    for g=630:650
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(32,j)*ori(g,j);
            sum3=sum3+optimalweight(32,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=630:650
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end       
    
        % 第33个期限 
    for g=651:672
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(33,j)*ori(g,j);
            sum3=sum3+optimalweight(33,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=651:672
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end       
    
          % 第34个期限  
    for g=673:688
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(34,j)*ori(g,j);
            sum3=sum3+optimalweight(34,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=673:688
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end       

   % 第35个期限 
    for g=689:709
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(35,j)*ori(g,j);
            sum3=sum3+optimalweight(35,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=689:709
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end       
  
   % 第36个期限  
    for g=710:732
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(36,j)*ori(g,j);
            sum3=sum3+optimalweight(36,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=710:732
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     
    
      % 第37个期限  
    for g=733:752
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(37,j)*ori(g,j);
            sum3=sum3+optimalweight(37,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=733:752
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     
    
          % 第38个期限 
    for g=753:767
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(38,j)*ori(g,j);
            sum3=sum3+optimalweight(38,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=753:767
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
         % 第39个期限 
    for g=768:790
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(39,j)*ori(g,j);
            sum3=sum3+optimalweight(39,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=768:790
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
       % 第40个期限 
    for g=791:811
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(40,j)*ori(g,j);
            sum3=sum3+optimalweight(40,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=791:811
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
   % 第41个期限 
    for g=812:831
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(41,j)*ori(g,j);
            sum3=sum3+optimalweight(41,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=812:831
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
       % 第42个期限 
    for g=832:850
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(42,j)*ori(g,j);
            sum3=sum3+optimalweight(42,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=832:850
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
          % 第43个期限 
    for g=851:872
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(43,j)*ori(g,j);
            sum3=sum3+optimalweight(43,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=851:872
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
        % 第44个期限 
    for g=873:894
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(44,j)*ori(g,j);
            sum3=sum3+optimalweight(44,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=873:894
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
    % 第45个期限 
    for g=895:913
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(45,j)*ori(g,j);
            sum3=sum3+optimalweight(45,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=895:913
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
       % 第46个期限 
    for g=914:929
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(46,j)*ori(g,j);
            sum3=sum3+optimalweight(46,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=914:929
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
     % 第47个期限 
    for g=930:951
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(47,j)*ori(g,j);
            sum3=sum3+optimalweight(47,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=930:951
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
        % 第48个期限 
    for g=952:974
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(48,j)*ori(g,j);
            sum3=sum3+optimalweight(48,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=952:974
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
           % 第49个期限 
    for g=975:994
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(49,j)*ori(g,j);
            sum3=sum3+optimalweight(49,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=975:994
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
               % 第50个期限 
    for g=995:1009
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(50,j)*ori(g,j);
            sum3=sum3+optimalweight(50,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=995:1009
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
      % 第51个期限 
    for g=1010:1032
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(51,j)*ori(g,j);
            sum3=sum3+optimalweight(51,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1010:1032
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
       % 第52个期限 
    for g=1033:1051
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(52,j)*ori(g,j);
            sum3=sum3+optimalweight(52,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1033:1051
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
        % 第53个期限 
    for g=1052:1072
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(53,j)*ori(g,j);
            sum3=sum3+optimalweight(53,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1052:1072
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
     % 第54个期限 
    for g=1073:1093
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(54,j)*ori(g,j);
            sum3=sum3+optimalweight(54,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1073:1093
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
      % 第55个期限 
    for g=1094:1114
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(55,j)*ori(g,j);
            sum3=sum3+optimalweight(55,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1094:1114
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
       % 第56个期限 
    for g=1115:1137
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(56,j)*ori(g,j);
            sum3=sum3+optimalweight(56,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1115:1137
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
       % 第57个期限 
    for g=1138:1158
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(57,j)*ori(g,j);
            sum3=sum3+optimalweight(57,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1138:1158
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
        % 第58个期限 
    for g=1159:1174
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(58,j)*ori(g,j);
            sum3=sum3+optimalweight(58,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1159:1174
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
         % 第59个期限 
    for g=1175:1196
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(59,j)*ori(g,j);
            sum3=sum3+optimalweight(59,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1175:1196
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
    % 第60个期限 
    for g=1197:1218
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(60,j)*ori(g,j);
            sum3=sum3+optimalweight(60,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1197:1218
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
     % 第61个期限 
    for g=1219:1233
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(61,j)*ori(g,j);
            sum3=sum3+optimalweight(61,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1219:1233
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
      % 第62个期限 
    for g=1234:1254
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(62,j)*ori(g,j);
            sum3=sum3+optimalweight(62,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1234:1254
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
      % 第63个期限 
    for g=1255:1276
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(63,j)*ori(g,j);
            sum3=sum3+optimalweight(63,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1255:1276
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
  
        % 第64个期限 
    for g=1277:1293
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(64,j)*ori(g,j);
            sum3=sum3+optimalweight(64,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1277:1293
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
    % 第65个期限 
    for g=1294:1315
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(65,j)*ori(g,j);
            sum3=sum3+optimalweight(65,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1294:1315
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
      % 第66个期限 
    for g=1316:1335
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(66,j)*ori(g,j);
            sum3=sum3+optimalweight(66,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1316:1335
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
        % 第67个期限 
    for g=1336:1357
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(67,j)*ori(g,j);
            sum3=sum3+optimalweight(67,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1336:1357
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end    
    
         % 第68个期限 
    for g=1358:1380
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(68,j)*ori(g,j);
            sum3=sum3+optimalweight(68,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1358:1380
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
    
       % 第69个期限 
    for g=1381:1400
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(69,j)*ori(g,j);
            sum3=sum3+optimalweight(69,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1381:1400
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end       
    
       % 第70个期限 
    for g=1401:1418
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(70,j)*ori(g,j);
            sum3=sum3+optimalweight(70,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1401:1418
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     
    
     % 第71个期限 
    for g=1419:1440
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(71,j)*ori(g,j);
            sum3=sum3+optimalweight(71,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1419:1440
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     

      % 第72个期限 
    for g=1441:1461
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(72,j)*ori(g,j);
            sum3=sum3+optimalweight(72,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1441:1461
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
    
       % 第73个期限 
    for g=1462:1481
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(73,j)*ori(g,j);
            sum3=sum3+optimalweight(73,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1462:1481
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     
    
    % 第74个期限 
    for g=1482:1496
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(74,j)*ori(g,j);
            sum3=sum3+optimalweight(74,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1482:1496
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     
 
      % 第75个期限 
    for g=1497:1517
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(75,j)*ori(g,j);
            sum3=sum3+optimalweight(75,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1497:1517
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  

       % 第76个期限 
    for g=1518:1535
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(76,j)*ori(g,j);
            sum3=sum3+optimalweight(76,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1518:1535
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  

        % 第77个期限 
    for g=1536:1557
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(77,j)*ori(g,j);
            sum3=sum3+optimalweight(77,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1536:1557
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
      % 第78个期限 
    for g=1558:1574
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(78,j)*ori(g,j);
            sum3=sum3+optimalweight(78,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1558:1574
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end    
    
        % 第79个期限 
    for g=1575:1597
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(79,j)*ori(g,j);
            sum3=sum3+optimalweight(79,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1575:1597
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
      % 第80个期限 
    for g=1598:1619
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(80,j)*ori(g,j);
            sum3=sum3+optimalweight(80,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1598:1619
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
   % 第81个期限 
    for g=1620:1638
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(81,j)*ori(g,j);
            sum3=sum3+optimalweight(81,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1620:1638
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
  % 第82个期限 
    for g=1639:1656
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(82,j)*ori(g,j);
            sum3=sum3+optimalweight(82,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1639:1656
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 

   % 第83个期限 
    for g=1657:1677
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(83,j)*ori(g,j);
            sum3=sum3+optimalweight(83,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1657:1677
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end    

     % 第84个期限 
    for g=1678:1699
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(84,j)*ori(g,j);
            sum3=sum3+optimalweight(84,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1678:1699
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
 
      % 第85个期限 
    for g=1700:1720
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(85,j)*ori(g,j);
            sum3=sum3+optimalweight(85,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1700:1720
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
 
        % 第86个期限 
    for g=1721:1736
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(86,j)*ori(g,j);
            sum3=sum3+optimalweight(86,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1721:1736
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 

          % 第87个期限 
    for g=1737:1757
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(87,j)*ori(g,j);
            sum3=sum3+optimalweight(87,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1737:1757
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
      % 第88个期限 
    for g=1758:1778
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(88,j)*ori(g,j);
            sum3=sum3+optimalweight(88,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1758:1778
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     
    
      % 第89个期限 
    for g=1779:1798
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(89,j)*ori(g,j);
            sum3=sum3+optimalweight(89,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1779:1798
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     
 
    % 第90个期限 
    for g=1799:1818
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(90,j)*ori(g,j);
            sum3=sum3+optimalweight(90,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1799:1818
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end        

   % 第91个期限 
    for g=1819:1841
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(91,j)*ori(g,j);
            sum3=sum3+optimalweight(91,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1819:1841
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     

  % 第92个期限 
    for g=1842:1862
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(92,j)*ori(g,j);
            sum3=sum3+optimalweight(92,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1842:1862
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end      
 
  % 第93个期限 
    for g=1863:1883
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(93,j)*ori(g,j);
            sum3=sum3+optimalweight(93,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1863:1883
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end      
    
    % 第94个期限 
    for g=1884:1901
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(94,j)*ori(g,j);
            sum3=sum3+optimalweight(94,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1884:1901
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
  
 % 第95个期限 
    for g=1902:1921
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(95,j)*ori(g,j);
            sum3=sum3+optimalweight(95,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1902:1921
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   

    % 第96个期限 
    for g=1922:1944
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(96,j)*ori(g,j);
            sum3=sum3+optimalweight(96,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1922:1944
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
 
   % 第97个期限 
    for g=1945:1964
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(97,j)*ori(g,j);
            sum3=sum3+optimalweight(97,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1945:1964
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end      

     % 第98个期限 
    for g=1965:1979
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(98,j)*ori(g,j);
            sum3=sum3+optimalweight(98,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1965:1979
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
   
      % 第99个期限 
    for g=1980:2001
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(99,j)*ori(g,j);
            sum3=sum3+optimalweight(99,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=1980:2001
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end    

        % 第100个期限 
    for g=2002:2022
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(100,j)*ori(g,j);
            sum3=sum3+optimalweight(100,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2002:2022
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
 
       % 第101个期限 
    for g=2023:2042
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(101,j)*ori(g,j);
            sum3=sum3+optimalweight(101,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2023:2042
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
 
      % 第102个期限 
    for g=2043:2063
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(102,j)*ori(g,j);
            sum3=sum3+optimalweight(102,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2043:2063
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
       % 第103个期限 
    for g=2064:2086
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(103,j)*ori(g,j);
            sum3=sum3+optimalweight(103,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2064:2086
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
    
    % 第104个期限 
    for g=2087:2107
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(104,j)*ori(g,j);
            sum3=sum3+optimalweight(104,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2087:2107
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
  
      % 第105个期限 
    for g=2108:2127
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(105,j)*ori(g,j);
            sum3=sum3+optimalweight(105,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2108:2127
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end    
    
     % 第106个期限 
    for g=2128:2144
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(106,j)*ori(g,j);
            sum3=sum3+optimalweight(106,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2128:2144
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end       
    
    % 第107个期限 
    for g=2145:2165
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(107,j)*ori(g,j);
            sum3=sum3+optimalweight(107,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2145:2165
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end

      % 第108个期限 
    for g=2166:2188
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(108,j)*ori(g,j);
            sum3=sum3+optimalweight(108,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2166:2188
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end
 
    % 第109个期限 
    for g=2189:2208
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(109,j)*ori(g,j);
            sum3=sum3+optimalweight(109,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2189:2208
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
    
  % 第110个期限 
    for g=2209:2224
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(110,j)*ori(g,j);
            sum3=sum3+optimalweight(110,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2209:2224
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end      
 
 % 第111个期限 
    for g=2225:2247
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(111,j)*ori(g,j);
            sum3=sum3+optimalweight(111,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2225:2247
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end       

% 第112个期限 
    for g=2248:2267
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(112,j)*ori(g,j);
            sum3=sum3+optimalweight(112,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2248:2267
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     

    % 第113个期限 
    for g=2268:2288
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(113,j)*ori(g,j);
            sum3=sum3+optimalweight(113,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2268:2288
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   

   % 第114个期限 
    for g=2289:2308
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(114,j)*ori(g,j);
            sum3=sum3+optimalweight(114,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2289:2308
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end       

      % 第115个期限 
    for g=2309:2329
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(115,j)*ori(g,j);
            sum3=sum3+optimalweight(115,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2309:2329
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
       % 第116个期限 
    for g=2330:2352
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(116,j)*ori(g,j);
            sum3=sum3+optimalweight(116,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2330:2352
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
 
    % 第117个期限 
    for g=2353:2372
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(117,j)*ori(g,j);
            sum3=sum3+optimalweight(117,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2353:2372
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 

      % 第118个期限 
    for g=2373:2388
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(118,j)*ori(g,j);
            sum3=sum3+optimalweight(118,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2373:2388
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
 
      % 第119个期限 
    for g=2389:2410
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(119,j)*ori(g,j);
            sum3=sum3+optimalweight(119,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2389:2410
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
      % 第120个期限 
    for g=2411:2432
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(120,j)*ori(g,j);
            sum3=sum3+optimalweight(120,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2411:2432
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
 
         % 第121个期限 
    for g=2433:2450
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(121,j)*ori(g,j);
            sum3=sum3+optimalweight(121,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2433:2450
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
  
          % 第122个期限 
    for g=2451:2468
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(122,j)*ori(g,j);
            sum3=sum3+optimalweight(122,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2451:2468
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 

    % 第123个期限 
    for g=2469:2491
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(123,j)*ori(g,j);
            sum3=sum3+optimalweight(123,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2469:2491
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
 
      % 第124个期限 
    for g=2492:2509
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(124,j)*ori(g,j);
            sum3=sum3+optimalweight(124,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2492:2509
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
  
        % 第125个期限 
    for g=2510:2529
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(125,j)*ori(g,j);
            sum3=sum3+optimalweight(125,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2510:2529
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 

          % 第126个期限 
    for g=2530:2551
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(126,j)*ori(g,j);
            sum3=sum3+optimalweight(126,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2530:2551
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
 
   % 第127个期限 
    for g=2552:2572
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(127,j)*ori(g,j);
            sum3=sum3+optimalweight(127,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2552:2572
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
 
      % 第128个期限 
    for g=2573:2595
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(128,j)*ori(g,j);
            sum3=sum3+optimalweight(128,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2573:2595
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 

       % 第129个期限 
    for g=2596:2616
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(129,j)*ori(g,j);
            sum3=sum3+optimalweight(129,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2596:2616
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 

     % 第130个期限 
    for g=2617:2633
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(130,j)*ori(g,j);
            sum3=sum3+optimalweight(130,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2617:2633
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end
  
     % 第131个期限
    for g=2634:2655
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(131,j)*ori(g,j);
            sum3=sum3+optimalweight(131,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2634:2655
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end
     
     % 第132个期限
    for g=2656:2676
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(132,j)*ori(g,j);
            sum3=sum3+optimalweight(132,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2656:2676
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end
    
  % 第133个期限
    for g=2677:2698
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(133,j)*ori(g,j);
            sum3=sum3+optimalweight(133,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2677:2698
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end    
     
   % 第134个期限
    for g=2699:2713
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(134,j)*ori(g,j);
            sum3=sum3+optimalweight(134,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2699:2713
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end       
    
  % 第135个期限
    for g=2714:2735
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(135,j)*ori(g,j);
            sum3=sum3+optimalweight(135,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2714:2735
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
    
 % 第136个期限
    for g=2736:2753
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(136,j)*ori(g,j);
            sum3=sum3+optimalweight(136,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2736:2753
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
    
% 第137个期限
    for g=2754:2775
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(137,j)*ori(g,j);
            sum3=sum3+optimalweight(137,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2754:2775
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end       
   
 % 第138个期限
    for g=2776:2795
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(138,j)*ori(g,j);
            sum3=sum3+optimalweight(138,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2776:2795
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end       
 
    % 第139个期限
    for g=2796:2817
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(139,j)*ori(g,j);
            sum3=sum3+optimalweight(139,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2796:2817
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
   
   % 第140个期限
    for g=2818:2840
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(140,j)*ori(g,j);
            sum3=sum3+optimalweight(140,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2818:2840
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end    
  
  % 第141个期限
    for g=2841:2859
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(141,j)*ori(g,j);
            sum3=sum3+optimalweight(141,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2841:2859
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end       
 
     % 第142个期限
    for g=2860:2877
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(142,j)*ori(g,j);
            sum3=sum3+optimalweight(142,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2860:2877
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
 
     % 第143个期限
    for g=2878:2899
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(143,j)*ori(g,j);
            sum3=sum3+optimalweight(143,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2878:2899
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
      % 第144个期限
    for g=2900:2919
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(144,j)*ori(g,j);
            sum3=sum3+optimalweight(144,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2900:2919
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
    
    % 第145个期限
    for g=2920:2941
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(145,j)*ori(g,j);
            sum3=sum3+optimalweight(145,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2920:2941
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     
 
  % 第146个期限
    for g=2942:2956
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(146,j)*ori(g,j);
            sum3=sum3+optimalweight(146,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2942:2956
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
    
 % 第147个期限
    for g=2957:2977
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(147,j)*ori(g,j);
            sum3=sum3+optimalweight(147,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2957:2977
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
 % 第148个期限
    for g=2978:2998
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(148,j)*ori(g,j);
            sum3=sum3+optimalweight(148,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2978:2998
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
    % 第149个期限
    for g=2999:3018
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(149,j)*ori(g,j);
            sum3=sum3+optimalweight(149,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=2999:3018
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
   % 第150个期限
    for g=3019:3037
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(150,j)*ori(g,j);
            sum3=sum3+optimalweight(150,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3019:3037
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
  
    % 第151个期限
    for g=3038:3060
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(151,j)*ori(g,j);
            sum3=sum3+optimalweight(151,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3038:3060
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
 
     % 第152个期限
    for g=3061:3082
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(152,j)*ori(g,j);
            sum3=sum3+optimalweight(152,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3061:3082
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
  
      % 第153个期限
    for g=3083:3102
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(153,j)*ori(g,j);
            sum3=sum3+optimalweight(153,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3083:3102
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
  
      % 第154个期限
    for g=3103:3120
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(154,j)*ori(g,j);
            sum3=sum3+optimalweight(154,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3103:3120
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
 
      % 第155个期限
    for g=3121:3141
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(155,j)*ori(g,j);
            sum3=sum3+optimalweight(155,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3121:3141
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end    
    
    % 第156个期限
    for g=3142:3163
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(156,j)*ori(g,j);
            sum3=sum3+optimalweight(156,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3142:3163
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end    
  
      % 第157个期限
    for g=3164:3179
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(157,j)*ori(g,j);
            sum3=sum3+optimalweight(157,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3164:3179
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end    
 
        % 第158个期限
    for g=3180:3199
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(158,j)*ori(g,j);
            sum3=sum3+optimalweight(158,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3180:3199
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
   
         % 第159个期限
    for g=3200:3221
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(159,j)*ori(g,j);
            sum3=sum3+optimalweight(159,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3200:3221
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
 
    % 第160个期限
    for g=3222:3242
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(160,j)*ori(g,j);
            sum3=sum3+optimalweight(160,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3222:3242
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
 
    % 第161个期限
      for g=3243:3260
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(161,j)*ori(g,j);
            sum3=sum3+optimalweight(161,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3243:3260
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     
 
     % 第162个期限
      for g=3261:3280
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(162,j)*ori(g,j);
            sum3=sum3+optimalweight(162,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3261:3280
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end    
    
    % 第163个期限
      for g=3281:3303
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(163,j)*ori(g,j);
            sum3=sum3+optimalweight(163,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3281:3303
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     
   
      % 第164个期限
      for g=3304:3324
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(164,j)*ori(g,j);
            sum3=sum3+optimalweight(164,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3304:3324
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end     
   
      % 第165个期限
      for g=3325:3346
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(165,j)*ori(g,j);
            sum3=sum3+optimalweight(165,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3325:3346
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   

     % 第166个期限
      for g=3347:3362
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(166,j)*ori(g,j);
            sum3=sum3+optimalweight(166,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3347:3362
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
 
     % 第167个期限
      for g=3363:3383
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(167,j)*ori(g,j);
            sum3=sum3+optimalweight(167,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3363:3383
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
 
      % 第168个期限
      for g=3384:3406
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(168,j)*ori(g,j);
            sum3=sum3+optimalweight(168,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3384:3406
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
 
     % 第169个期限
     for g=3407:3426
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(169,j)*ori(g,j);
            sum3=sum3+optimalweight(169,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3407:3426
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end   
  
     % 第170个期限
     for g=3427:3441
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(170,j)*ori(g,j);
            sum3=sum3+optimalweight(170,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3427:3441
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
  
     % 第171个期限
     for g=3442:3464
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(171,j)*ori(g,j);
            sum3=sum3+optimalweight(171,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3442:3464
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
  
     % 第172个期限
     for g=3465:3485
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(172,j)*ori(g,j);
            sum3=sum3+optimalweight(172,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3465:3485
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
  
      % 第173个期限
     for g=3486:3503
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(173,j)*ori(g,j);
            sum3=sum3+optimalweight(173,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3486:3503
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end  
   
     % 第174个期限
     for g=3504:3524
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(174,j)*ori(g,j);
            sum3=sum3+optimalweight(174,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3504:3524
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
   % 第175个期限
     for g=3525:3546
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(175,j)*ori(g,j);
            sum3=sum3+optimalweight(175,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3525:3546
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end    

     % 第176个期限
     for g=3547:3568
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(176,j)*ori(g,j);
            sum3=sum3+optimalweight(176,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3547:3568
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
  % 第177个期限
     for g=3569:3588
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(177,j)*ori(g,j);
            sum3=sum3+optimalweight(177,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3569:3588
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 

     % 第178个期限
       for g=3589:3604
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(178,j)*ori(g,j);
            sum3=sum3+optimalweight(178,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3589:3604
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
 
     % 第179个期限
       for g=3605:3626
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(179,j)*ori(g,j);
            sum3=sum3+optimalweight(179,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3605:3626
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
    
    % 第180个期限
       for g=3627:3649
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori(g,j);
            sum2=sum2+optimalweight2(180,j)*ori(g,j);
            sum3=sum3+optimalweight(180,j)*ori(g,j);
        end
        Pret1(g,1)=sum1;  
        Pret2(g,1)=sum2;  
        Pret3(g,1)=sum3;  
   end
    %计算累计财富
    for g=3627:3649
        wealth1(g+1,1)=Pret1(g,1)*wealth1(g,1);
        wealth2(g+1,1)=Pret2(g,1)*wealth2(g,1);
        wealth3(g+1,1)=Pret3(g,1)*wealth3(g,1);
    end 
   
    
figure(1)  
plot(wealth1(2:end),'black','linewidth',1.0)
hold on
plot(wealth2(2:end),'blue','linewidth',1.0)
hold on
plot(wealth3(2:end),'red','linewidth',1.0)
legend('EW','OR','MSGCE'); 
ylabel('Achieved wealth') %y轴坐标描述
title('Out of sample') %y轴坐标描述

%计算最大回撤
wealth_out=[wealth1 wealth2 wealth3];
MaxDD_out= maxdrawdown(wealth_out);

%计算组合收益
Pret_out=[Pret1 Pret2 Pret3];
%计算组合平均收益
meanPret_out=mean(Pret_out);

%计算方差 偏度和峰度
varPret_out=var(Pret_out);

skewPret_out=skewness(Pret_out);

kurtPret_out=kurtosis(Pret_out);

bench_out=1.000210998-1;%上证指数样本外的平均收益率  

%计算夏普比率
stdPret_out=std(Pret_out);

SR_out=(meanPret_out-1)./stdPret_out;

%计算Carmar比率
Carmar_out=(meanPret_out-1)./MaxDD_out;


%计算omega比率和索提诺比率
 for m=1:3
       
    for g=1:3649
        if  Pret_out(g,m)-1<bench_out
            a1_outsample(g,m)=bench_out-(Pret_out(g,m)-1);
            a2_outsample(g,m)=(Pret_out(g,m)-1-bench_out)^2;
        else   
            a1_outsample(g,m)=0;
            a2_outsample(g,m)=0;
        end
    end
    
 
       
     for g=1:3649
        if Pret_out(g,m)-1>bench_out
           b1(g,m)=Pret_out(g,m)-1-bench_out;
        else
           b1(g,m)=0;
        end

        if Pret_out(g,m)-1<bench_out
           b2(g,m)=(Pret_out(g,m)-1-bench_out)^2;
        else
           b2(g,m)=0;
        end
    end
    UPR_out(1,m)=mean(b1(:,m))/(mean(b2(:,m))^0.5);
    
    
OmegaR_out(1,m)=(meanPret_out(1,m)-1-bench_out)/mean(a1_outsample(:,m))+1;

D2_out(1,m)=sum(a2_outsample(:,m))/3649; %组合收益半方差  
SortinoR_out(1,m)=(meanPret_out(1,m)-1-bench_out)/(D2_out(1,m)^0.5);

 end
 
performance.outsample=[(meanPret_out(1,:)-1)*100;varPret_out(1,:)*100;skewPret_out(1,:);kurtPret_out(1,:);SR_out(1,:);SortinoR_out(1,:);OmegaR_out(1,:);MaxDD_out(1,:);Carmar_out(1,:)];


%% 08金融危机
%样本外
ori_crisis=xlsread('20只能源股收益20060101_20070629month.xlsx',2,'B1:U428'); %样本外全市场
%% 计算组合收益
% 第一个期限 
    for g=1:22
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(7,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(7,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
    end
    %计算累计财富
    wealth1_crisis=1;wealth2_crisis=1;wealth3_crisis=1;
    for g=1:22
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end

%第二个期限 
   for g=23:45
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(8,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(8,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=23:45
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end 
    
 % 第三个期限    
     for g=46:65
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(9,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(9,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=46:65
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end 
    
    % 第四个期限   
     for g=66:83
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(10,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(10,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
    end
    %计算累计财富
    for g=66:83
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end 
    
        % 第五个期限    
     for g=84:105
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(11,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(11,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=84:105
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end 
    
       % 第六个期限    
     for g=106:125
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(12,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(12,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=106:125
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end  
    
    % 第七个期限  
    for g=126:147
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(13,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(13,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=126:147
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end
    
     % 第8个期限  
    for g=148:163
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(14,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(14,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=148:163
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end 
    
    % 第9个期限  
    for g=164:184
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(15,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(15,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=164:184
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end 
    
   % 第10个期限  
    for g=185:205
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(16,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(16,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=185:205
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end 
    
    % 第11个期限  
    for g=206:225
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(17,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(17,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=206:225
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end 
  
     % 第12个期限
  for g=226:245
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(18,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(18,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=226:245
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end   
 
   % 第13个期限
  for g=246:268
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(19,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(19,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=246:268
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end   
    
     % 第14个期限
  for g=269:289
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(20,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(20,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=269:289
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end   
 
  % 第15个期限
  for g=290:308
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(21,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(21,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=290:308
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end   
 
    % 第16个期限
  for g=309:328
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(22,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(22,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=309:328
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end  
 
      % 第17个期限
  for g=329:348
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(23,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(23,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=329:348
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end  
   
        % 第18个期限
  for g=349:371
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(24,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(24,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=349:371
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end  
  
   % 第19个期限
  for g=372:386
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(25,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(25,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=372:386
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end  
    
     % 第20个期限
     for g=387:406
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(26,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(26,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=387:406
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end  
    
      % 第21个期限
     for g=407:428
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_crisis(g,j);
            sum2=sum2+optimalweight2(27,j)*ori_crisis(g,j);
            sum3=sum3+optimalweight(27,j)*ori_crisis(g,j);
        end
        Pret1_crisis(g,1)=sum1;  
        Pret2_crisis(g,1)=sum2;  
        Pret3_crisis(g,1)=sum3;  
   end
    %计算累计财富
    for g=407:428
        wealth1_crisis(g+1,1)=Pret1_crisis(g,1)*wealth1_crisis(g,1);
        wealth2_crisis(g+1,1)=Pret2_crisis(g,1)*wealth2_crisis(g,1);
        wealth3_crisis(g+1,1)=Pret3_crisis(g,1)*wealth3_crisis(g,1);
    end  
    
figure(2)  
plot(wealth1_crisis(2:end),'black','linewidth',1.0)
hold on
plot(wealth2_crisis(2:end),'blue','linewidth',1.0)
hold on
plot(wealth3_crisis(2:end),'red','linewidth',1.0)
legend('EW','OR','MSGCE'); 
ylabel('Achieved wealth') %y轴坐标描述
title('Out of sample-crisis') %y轴坐标描述

%计算最大回撤
wealth_out_crisis=[wealth1_crisis wealth2_crisis wealth3_crisis];
MaxDD_out_crisis= maxdrawdown(wealth_out_crisis);

%计算组合收益
Pret_out_crisis=[Pret1_crisis Pret2_crisis Pret3_crisis];
%计算组合平均收益
meanPret_out_crisis=mean(Pret_out_crisis);

%计算方差 偏度和峰度
varPret_out_crisis=var(Pret_out_crisis);
skewPret_out_crisis=skewness(Pret_out_crisis);
kurtPret_out_crisis=kurtosis(Pret_out_crisis);

bench_out_crisis=0.998286136-1;%上证指数样本外的平均收益率  

%计算夏普比率
stdPret_out_crisis=std(Pret_out_crisis);

SR_out_crisis=(meanPret_out_crisis-1)./stdPret_out_crisis;

%计算Carmar比率
Carmar_out_crisis=(meanPret_out_crisis-1)./MaxDD_out_crisis;

%计算omega比率和索提诺比率
 for m=1:3
    
    for g=1:428
        if  Pret_out_crisis(g,m)-1<bench_out_crisis
            a1_outsample_crisis(g,m)=bench_out_crisis-(Pret_out_crisis(g,m)-1);
            a2_outsample_crisis(g,m)=(Pret_out_crisis(g,m)-1-bench_out_crisis)^2;
        else   
            a1_outsample_crisis(g,m)=0;
            a2_outsample_crisis(g,m)=0;
        end
    end
    
     %计算upside potential ratio比率
    for g=1:428
        if Pret_out_crisis(g,m)-1>bench_out_crisis
           b1(g,m)=Pret_out_crisis(g,m)-1-bench_out_crisis;
        else
           b1(g,m)=0;
        end

        if Pret_out_crisis(g,m)-1<bench_out_crisis
           b2(g,m)=(Pret_out_crisis(g,m)-1-bench_out_crisis)^2;
        else
           b2(g,m)=0;
        end
    end
    UPR_out_crisis(1,m)=mean(b1(:,m))/(mean(b2(:,m))^0.5);
    OmegaR_out_crisis(1,m)=(meanPret_out_crisis(1,m)-1-bench_out_crisis)/mean(a1_outsample_crisis(:,m))+1;

D2_out_crisis(1,m)=sum(a2_outsample_crisis(:,m))/428; %组合收益半方差  
SortinoR_out_crisis(1,m)=(meanPret_out_crisis(1,m)-1-bench_out_crisis)/(D2_out_crisis(1,m)^0.5);

 end
 
performance.outsample_crisis=[(meanPret_out_crisis(1,:)-1)*100;varPret_out_crisis(1,:)*100;skewPret_out_crisis(1,:);kurtPret_out_crisis(1,:);SR_out_crisis(1,:);SortinoR_out_crisis(1,:);OmegaR_out_crisis(1,:);MaxDD_out_crisis(1,:);Carmar_out_crisis(1,:)];

%% 疫情时间段
%样本外
ori_COVID=xlsread('20只能源股收益20200101_20200430month.xlsx',2,'B1:U79'); %样本外全市场
%% 计算组合收益
% 第一个期限 
    for g=1:16
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_COVID(g,j);
            sum2=sum2+optimalweight2(157,j)*ori_COVID(g,j);
            sum3=sum3+optimalweight(157,j)*ori_COVID(g,j);
        end
        Pret1_COVID(g,1)=sum1;  
        Pret2_COVID(g,1)=sum2;  
        Pret3_COVID(g,1)=sum3;  
    end
    %计算累计财富
    wealth1_COVID=1;wealth2_COVID=1;wealth3_COVID=1;
    for g=1:16
        wealth1_COVID(g+1,1)=Pret1_COVID(g,1)*wealth1_COVID(g,1);
        wealth2_COVID(g+1,1)=Pret2_COVID(g,1)*wealth2_COVID(g,1);
        wealth3_COVID(g+1,1)=Pret3_COVID(g,1)*wealth3_COVID(g,1);
    end

  % 第2个期限 
    for g=17:36
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_COVID(g,j);
            sum2=sum2+optimalweight2(158,j)*ori_COVID(g,j);
            sum3=sum3+optimalweight(158,j)*ori_COVID(g,j);
        end
        Pret1_COVID(g,1)=sum1;  
        Pret2_COVID(g,1)=sum2;  
        Pret3_COVID(g,1)=sum3;  
    end
    %计算累计财富
    for g=17:36
        wealth1_COVID(g+1,1)=Pret1_COVID(g,1)*wealth1_COVID(g,1);
        wealth2_COVID(g+1,1)=Pret2_COVID(g,1)*wealth2_COVID(g,1);
        wealth3_COVID(g+1,1)=Pret3_COVID(g,1)*wealth3_COVID(g,1);
    end
  
   % 第3个期限 
    for g=37:58
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_COVID(g,j);
            sum2=sum2+optimalweight2(159,j)*ori_COVID(g,j);
            sum3=sum3+optimalweight(159,j)*ori_COVID(g,j);
        end
        Pret1_COVID(g,1)=sum1;  
        Pret2_COVID(g,1)=sum2;  
        Pret3_COVID(g,1)=sum3;  
    end
    %计算累计财富
    for g=37:58
        wealth1_COVID(g+1,1)=Pret1_COVID(g,1)*wealth1_COVID(g,1);
        wealth2_COVID(g+1,1)=Pret2_COVID(g,1)*wealth2_COVID(g,1);
        wealth3_COVID(g+1,1)=Pret3_COVID(g,1)*wealth3_COVID(g,1);
    end 
    
  % 第4个期限 
    for g=59:79
        sum1=0;sum2=0;sum3=0;
        for j=1:20   
            sum1=sum1+0.05*ori_COVID(g,j);
            sum2=sum2+optimalweight2(160,j)*ori_COVID(g,j);
            sum3=sum3+optimalweight(160,j)*ori_COVID(g,j);
        end
        Pret1_COVID(g,1)=sum1;  
        Pret2_COVID(g,1)=sum2;  
        Pret3_COVID(g,1)=sum3;  
    end
    %计算累计财富
    for g=59:79
        wealth1_COVID(g+1,1)=Pret1_COVID(g,1)*wealth1_COVID(g,1);
        wealth2_COVID(g+1,1)=Pret2_COVID(g,1)*wealth2_COVID(g,1);
        wealth3_COVID(g+1,1)=Pret3_COVID(g,1)*wealth3_COVID(g,1);
    end 
    
figure(3)
plot(wealth1_COVID(2:end),'black','linewidth',1.0)
hold on
plot(wealth2_COVID(2:end),'blue','linewidth',1.0)
hold on
plot(wealth3_COVID(2:end),'red','linewidth',1.0)
legend('EW','OR','MSGCE'); 
ylabel('Achieved wealth') %y轴坐标描述
title('Out of sample-COVID') %y轴坐标描述

%计算最大回撤
wealth_out_COVID=[wealth1_COVID wealth2_COVID wealth3_COVID];
MaxDD_out_COVID= maxdrawdown(wealth_out_COVID);

%计算组合收益
Pret_out_COVID=[Pret1_COVID Pret2_COVID Pret3_COVID];
%计算组合平均收益
meanPret_out_COVID=mean(Pret_out_COVID);

%计算方差 偏度和峰度
varPret_out_COVID=var(Pret_out_COVID);
skewPret_out_COVID=skewness(Pret_out_COVID);
kurtPret_out_COVID=kurtosis(Pret_out_COVID);

bench_out_COVID=0.999312944-1;%上证指数样本外的平均收益率  

%计算夏普比率
stdPret_out_COVID=std(Pret_out_COVID);
SR_out_COVID=(meanPret_out_COVID-1)./stdPret_out_COVID;

%计算Carmar比率
Carmar_out_COVID=(meanPret_out_COVID-1)./MaxDD_out_COVID;


%计算omega比率和索提诺比率
 for m=1:3
    
    for g=1:79
        if  Pret_out_COVID(g,m)-1<bench_out_COVID
            a1_outsample_COVID(g,m)=bench_out_COVID-(Pret_out_COVID(g,m)-1);
            a2_outsample_COVID(g,m)=(Pret_out_COVID(g,m)-1-bench_out_COVID)^2;
        else   
            a1_outsample_COVID(g,m)=0;
            a2_outsample_COVID(g,m)=0;
        end
    end
    
          %计算upside potential ratio比率
    for g=1:79
        if Pret_out_COVID(g,m)-1>bench_out_COVID
           b1(g,m)=Pret_out_COVID(g,m)-1-bench_out_COVID;
        else
           b1(g,m)=0;
        end

        if Pret_out_COVID(g,m)-1<bench_out_COVID
           b2(g,m)=(Pret_out_COVID(g,m)-1-bench_out_COVID)^2;
        else
           b2(g,m)=0;
        end
    end
    UPR_out_COVID(1,m)=mean(b1(:,m))/(mean(b2(:,m))^0.5);
    
OmegaR_out_COVID(1,m)=(meanPret_out_COVID(1,m)-1-bench_out_COVID)/mean(a1_outsample_COVID(:,m))+1;

D2_out_COVID(1,m)=sum(a2_outsample_COVID(:,m))/79; %组合收益半方差  
SortinoR_out_COVID(1,m)=(meanPret_out_COVID(1,m)-1-bench_out_COVID)/(D2_out_COVID(1,m)^0.5);

 end
 
performance.outsample_COVID=[(meanPret_out_COVID(1,:)-1)*100;varPret_out_COVID(1,:)*100;skewPret_out_COVID(1,:);kurtPret_out_COVID(1,:);SR_out_COVID(1,:);SortinoR_out_COVID(1,:);OmegaR_out_COVID(1,:);MaxDD_out_COVID(1,:);Carmar_out_COVID(1,:)];

toc