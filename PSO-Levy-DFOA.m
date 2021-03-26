function [bestPos,cg_curve]=PSO-Levy-DFOA(N,MaxFEs,lb,ub,dim,fobj)

%% 果蝇的相关参数
temp = 10;%果蝇寻优步长
noP=N;
wMax=0.99;%最大惯性权重值
wMin=0.1;%最小惯性权重值
c1=2;%学习因子c1
c2=2;%学习因子c2

%% 初始果蝇群体位置及相关参数
velx=rands(noP,dim);
vely=rands(noP,dim);
pBestScore=zeros(noP);
pBest=zeros(noP,dim);
gBest=zeros(1,dim);
fit=zeros(noP,1);
isFight=zeros(noP,1);
cg_curve=[];

%在全局范围内初始化种群个体
pos=initialization(noP,dim,ub,lb);
x_foa=1/sqrt(2)./pos;
y_foa=1/sqrt(2)./pos;
for i=1:noP
    pBestScore(i)=inf;
end
 gBestScore=inf;
it=1;    
FEs=0;
%纠正越界个体
Flag4ub=pos(i,:)>ub;
Flag4lb=pos(i,:)<lb;
pos(i,:)=(pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
%% 果蝇寻优开始
while  FEs < MaxFEs
    
    for i=1:size(pos,1)     
        %计算每个个体的目标函数适应度值，并更新最优个体信息
        if FEs<MaxFEs
            FEs=FEs+1;    
            fitness=fobj(pos(i,:));
            fit(i)=fitness;
            if(pBestScore(i)>fitness)
                pBestScore(i)=fitness;
                pBest(i,:)=pos(i,:);
                pBest_xfoa(i,:)=x_foa(i,:);
                pBest_yfoa(i,:)=y_foa(i,:);
            end
            if(gBestScore>fitness)
                gBestScore=fitness;
                gBest=pos(i,:);
                gBest_xfoa=x_foa(i,:);
                gBest_yfoa=x_foa(i,:);
            end
        else
            break;
        end
    end
    %寻找最优个体与最差个体
    fit_best=min(fit);
    fit_worst=max(fit);
    for i=1:noP
        %根据当前个体与最优个体和最差个体之间的差值来判定下次飞行采用PSO搜索策略还是Levy Fight
        if abs(fit(i)-fit_worst)>abs(fit(i)-fit_best)
            isFight(i)=1;
        end
    end
    %更新惯性权重 w 使其随迭代次数增加而减小
    w=wMax-(wMax-wMin)*MaxFEs;
    %生成Levy Fight 的飞行步长
    LF=LevyFight(noP,dim);
    for i=1:size(pos,1)
        for j=1:size(pos,2)   
            %较差子群采用PSO搜索策略进行位置更新
            if isFight(i)==0
                velx(i,j)=w*velx(i,j)+c1*rand()*(pBest_xfoa(i,j)- x_foa(i,j))+c2*rand()*(gBest_xfoa(j)-x_foa(i,j));
                vely(i,j)=w*vely(i,j)+c1*rand()*(pBest_yfoa(i,j)- y_foa(i,j))+c2*rand()*(gBest_yfoa(j)-y_foa(i,j));
                x_foa(i,j)=x_foa(i,j)+velx(i,j);
                y_foa(i,j)=y_foa(i,j)+vely(i,j);
                pos(i,j)=1/sqrt(x_foa(i,j)^2+y_foa(i,j)^2);
                if pos(i,j)>ub
                    pos(i,j)=ub;
                end
                if pos(i,j)<lb
                    pos(i,j)=lb;
                end
            else
                %较优子群采用Levy Fight进行位置更新
                x_foa(i,j)=gBest_xfoa(j)+temp*LF(i,j);
                y_foa(i,j)=gBest_yfoa(j)+temp*LF(i,j);
                pos(i,j)=1/sqrt(x_foa(i,j)^2+y_foa(i,j)^2);
            end
        end
    end
    %纠正越界个体
    pos(i,pos(i,:)>ub)=ub;
    pos(i,pos(i,:)<lb)=lb;
    cg_curve(it)=gBestScore;
    it=it+1;
    bestPos=gBest;
end

end