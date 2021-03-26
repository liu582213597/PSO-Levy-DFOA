function [bestPos,cg_curve]=PSO-Levy-DFOA(N,MaxFEs,lb,ub,dim,fobj)

%% ��Ӭ����ز���
temp = 10;%��ӬѰ�Ų���
noP=N;
wMax=0.99;%������Ȩ��ֵ
wMin=0.1;%��С����Ȩ��ֵ
c1=2;%ѧϰ����c1
c2=2;%ѧϰ����c2

%% ��ʼ��ӬȺ��λ�ü���ز���
velx=rands(noP,dim);
vely=rands(noP,dim);
pBestScore=zeros(noP);
pBest=zeros(noP,dim);
gBest=zeros(1,dim);
fit=zeros(noP,1);
isFight=zeros(noP,1);
cg_curve=[];

%��ȫ�ַ�Χ�ڳ�ʼ����Ⱥ����
pos=initialization(noP,dim,ub,lb);
x_foa=1/sqrt(2)./pos;
y_foa=1/sqrt(2)./pos;
for i=1:noP
    pBestScore(i)=inf;
end
 gBestScore=inf;
it=1;    
FEs=0;
%����Խ�����
Flag4ub=pos(i,:)>ub;
Flag4lb=pos(i,:)<lb;
pos(i,:)=(pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
%% ��ӬѰ�ſ�ʼ
while  FEs < MaxFEs
    
    for i=1:size(pos,1)     
        %����ÿ�������Ŀ�꺯����Ӧ��ֵ�����������Ÿ�����Ϣ
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
    %Ѱ�����Ÿ�����������
    fit_best=min(fit);
    fit_worst=max(fit);
    for i=1:noP
        %���ݵ�ǰ���������Ÿ����������֮��Ĳ�ֵ���ж��´η��в���PSO�������Ի���Levy Fight
        if abs(fit(i)-fit_worst)>abs(fit(i)-fit_best)
            isFight(i)=1;
        end
    end
    %���¹���Ȩ�� w ʹ��������������Ӷ���С
    w=wMax-(wMax-wMin)*MaxFEs;
    %����Levy Fight �ķ��в���
    LF=LevyFight(noP,dim);
    for i=1:size(pos,1)
        for j=1:size(pos,2)   
            %�ϲ���Ⱥ����PSO�������Խ���λ�ø���
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
                %������Ⱥ����Levy Fight����λ�ø���
                x_foa(i,j)=gBest_xfoa(j)+temp*LF(i,j);
                y_foa(i,j)=gBest_yfoa(j)+temp*LF(i,j);
                pos(i,j)=1/sqrt(x_foa(i,j)^2+y_foa(i,j)^2);
            end
        end
    end
    %����Խ�����
    pos(i,pos(i,:)>ub)=ub;
    pos(i,pos(i,:)<lb)=lb;
    cg_curve(it)=gBestScore;
    it=it+1;
    bestPos=gBest;
end

end