clear
load('数据.mat')
% 函数声明
global fun;global violate_fun;global D;
fun='obj_fun'; violate_fun='vio_g'; D=31;
global Q l L lowb upb MR Opt_sol w0 y1 y2 y3 c1 c2 sn

Q=data{1};Q{2}=Q{2}*1.1249;Q{3}=Q{3}*1.4802;
l=data{2};L=data{3};


MR=0.85;
limit=500;SPP=500;
sn=100;                              
MCN=2000;                            
Opt_sol=zeros(1,D);                  
Optimum=zeros(1,MCN);                
lowb=zeros(1,D);             
upb=[ones(1,D-11),ones(1,11)*7000];
w0=100;                              
y1=[2.1,3.15,10.5,5.25,21,15.75,7.35,5.25,10.5,18.9,5.25,7.35,15.75,36.75,26.25,15.75,17.85,21,36.75,52.5];
y2=[31.5,52.5,78.75,68.25,89.25,42,36.75,57.75,78.75,89.25,36.75,42,99.75,126,73.5,99.75,84,94.5,126,136.5];
y3=[31.5,47.25,105,178.5,210,47.25,31.5,52.5,84,126,31.5,63,189,294,105,157.5,189,231,262.5,367.5];
cmax=2;
cmin=0.5;
ispositive=1;    

vio=zeros(MCN,sn/2); 


max_yun=1;   
Opt=zeros(max_yun,D+3); 
for yun=1:max_yun
    tic
    failure=zeros(1,sn/2);  
    Col_fitness=zeros(1,sn/2);
    Col_f=zeros(1,sn/2);
    Col_violation=zeros(1,sn/2);
    Colony=zeros(sn/2,D);
    w=[w0,zeros(1,11)];
    for i=1:sn/2
        while 1
        Colony(i,1)=rand;
        for j=2:D-11
            Colony(i,j)=4*Colony(i,j-1)*(1-Colony(i,j-1));
        end
        for j=1:D-11
            if Colony(i,j)<=0.5
                Colony(i,j)=0;
            else
                Colony(i,j)=1;
            end
        end
        x=Colony(i,1:20);y=zeros(1,11);
        y11=y1*x'/3;y21=y2*x'*1.1249/7;y31=y3*x'*1.4802;
        c_cost(1:3)=y11;c_cost(4:10)=y21;c_cost(11)=y31;
        for j=1:11
            if w(j)>c_cost(j)
                y(j)=0;
                w(j+1)=(w(j)-c_cost(j))*1.025;
            else
                if j==1 || y(j-1)==0
                    y0=rand;
                else
                    y0=4*y0*(1-y0);
                end
                y(j)=y0*upb(20+j);
            end    
        end
        Colony(i,D-11+1:D)=y;
        [Col_fitness(i),Col_f(i),Col_violation(i),ispositive]=evaluation(Colony(i,:));
        if ispositive==1
            break;
        else
            disp('不满足非负');
        end
        end
    end    

    [max_fit,n]=max(Col_fitness);
    if max_fit>0
        Opt_sol=Colony(n,:); 
        min_vio=0;
        best_f=Col_f(n);
    else
        [min_vio,m]=min(Col_violation);max_fit=0;
        Opt_sol=Colony(m,:);
        best_f=Col_f(m);
    end

Cycle=1;
while Cycle<=MCN
    c1=(cmin-cmax)*Cycle/MCN+cmax;
    c2=(cmax-cmin)*Cycle/MCN+cmin;
    for i=1:sn/2
        while 1       
            v=modify(Colony,i);
            [v_fitness,v_f,v_violation,ispositive]=evaluation(v);
            if ispositive==1
                break;
            else
                disp('v不满足非负');
            end
        end
        isupdate=selection_Deb(Col_violation(i),v_violation,Col_f(i),v_f);
        if isupdate==0
            failure(i)=failure(i)+1;
        else
            Colony(i,:)=v;
            Col_fitness(i)=v_fitness;
            Col_f(i)=v_f;
            Col_violation(i)=v_violation;
            failure(i)=0;
        end
    end
    
p=zeros(1,sn/2);
for i=1:sn/2
    p(i)=calculate_pro(i,Col_fitness,Col_violation);
end
    

t=0;i=1;
while t<sn/2
    if i==0
        i=1;  
    end
    if rand<p(i)
        t=t+1;
        while 1
            v=modify(Colony,i);
            [v_fitness,v_f,v_violation,ispositive]=evaluation(v);
            if ispositive==1
                break;
            else 
                disp('v不满足非负');
            end
        end
        isupdate=selection_Deb(Col_violation(i),v_violation,Col_f(i),v_f);
        if isupdate==0
            failure(i)=failure(i)+1;
        else
            Colony(i,:)=v;
            Col_fitness(i)=v_fitness;
            Col_f(i)=v_f;
            Col_violation(i)=v_violation;
            failure(i)=0;
        end
    end
    i=i+1;
    i=mod(i,sn/2+1);
end


if mod(Cycle,SPP)==0
    for i=1:sn/2
    if failure(i)>limit
        while 1
            Colony(i,1:D-11)=round(rand(1,D-11));
            x=Colony(i,1:20);y=zeros(1,11);
            y11=y1*x'/3;y21=y2*x'*1.1249/7;y31=y3*x'*1.4802;
            c_cost(1:3)=y11;c_cost(4:10)=y21;c_cost(11)=y31;
            for j=1:11
                if w(j)>c_cost(j)
                    y(j)=0;
                    w(j+1)=(w(j)-c_cost(j))*1.025;
                else
                    y(j)=rand*upb(20+j);
                end
            end
            Colony(i,D-11+1:D)=y;
            failure(i)=0;
            [Col_fitness(i),Col_f(i),Col_violation(i),ispositive]=evaluation(Colony(i,:));
            if ispositive==1
                break;
            else
                disp('不满足非负');
            end
        end
    end
    end
end


vio(Cycle,:)=Col_violation;  


[max_fit1,n]=max(Col_fitness);
if max_fit1>0
    bestpop=Colony(n,:); 
    min_vio1=0;
    best_ff=Col_f(n);
else
    [min_vio1,m]=min(Col_violation);max_fit1=0;
    bestpop=Colony(m,:);
    best_ff=Col_f(m);
end

isupdate=selection_Deb(min_vio,min_vio1,best_f,best_ff);
if isupdate==1
    Opt_sol=bestpop;
    max_fit=max_fit1;min_vio=min_vio1;
end

if max_fit>0
spp=200;limit=200;
end
Optimum(Cycle)=max_fit;


Cycle=Cycle+1;
end

figure;
plot(Optimum,'b');
xlabel('Generation');
ylabel('fitness');
x=Opt_sol(1:20);y=Opt_sol(D-11+1:D);
y11=y1*x'/3;y21=y2*x'*1.1249/7;w=zeros(1,11);
w(1)=w0;
for i=2:4
    w(i)=(w(i-1)+y(i-1)-y11)*1.025-y(i-1)*1.045;
end
for i=5:11
    w(i)=(w(i-1)+y(i-1)-y21)*1.025-y(i-1)*1.045;
end
Opt_val=feval(fun,Opt_sol,w(11));
wealth=w(11)+4.0456*(52.5*x(1)+ 105*x(2)+ 210*x(3)+ 210*x(4)+ 630*x(5)+ 105*x(6)+ 84*x(7)+ 105*x(8)+ ...
    189*x(9)+ 399*x(10)+ 84*x(11)+ 105*x(12)+ 42*x(13)+ 735*x(14)+ 525*x(15)+ 315*x(16)+ 367.5*x(17)+...
    577.5*x(18)+ 840*x(19)+ 1207.5*x(20))- (31.5*x(1)+ 47.25*x(2)+ 105*x(3)+ 178.5*x(4)+ 210*x(5)+...
    47.25*x(6)+ 31.5*x(7)+ 52.5*x(8)+ 84*x(9)+ 126*x(10)+ 31.5*x(11)+ 63*x(12)+ 189*x(13)+...
    294*x(14)+ 105*x(15)+ 157.5*x(16)+ 189*x(17)+ 231*x(18)+ 262.5*x(19)+ 367.5*x(20))*1.4802;

disp('x的最优解=');disp(num2str(x));
disp('s的最优解=');disp(num2str(x));
disp('y的最优解=');disp(num2str(y));
disp('y的最优解的和=');disp(num2str(sum(y)));
disp('最优值=');disp(num2str(Opt_val));
disp('期末财富值=');disp(num2str(wealth));
yun
Opt(yun,:)=[Opt_sol,Opt_val,wealth,min_vio];
toc
sum(x)
end


    
            
                
                
    
        
    
            
        
                
            
        
            
        
                
            
            



        
