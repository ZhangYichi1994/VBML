%% Generate the evolutionary game data
% Input: Adjacent Matrx--Adj
% Output: Evolutionary game data--Stra and Unity
% Stra--Data of strategies, Stra(i,j)-the strategy of agent i in step j
% Unity--Data of fitness, Unity(i,j)-the fitness of agent i in step j

function [Stra,Unity] = Game_withnoise(Adj,b,amp)
N = size(Adj,1);% N为网络节点数目

Stra=[];
stra=(rand(N,1)>0.5)*1;% 策略初始分布，等概率的选择合作或者非合作的行为
Unity=[];% 适应度
for T=1:4*N
    for i=1:N
        player1=i;%从群体选择一个个体i
        stra_player1=stra(player1);%得到个体i的策略
        score1=0;
        Neig=[];
        for j=1:N
%             if Adj(player1,j)==1
            if Adj(player1,j) > 0
                Neig=[Neig,j];%记录个体player1的邻居
                player2=j;%得到个体的邻居及其策略
                stra_player2=stra(player2);
                score1=score1+Payoff_noise(stra_player1,stra_player2,b,amp);%得到个体player1的收益  
            end            
        end
        unity(i)=score1;%计算各个个体的收益
    end
    Unity=[Unity,unity'];%第i列为i时刻所有个体的收益
    
    Stra=[Stra,stra];
    for i=1:N        
        player1=i;%从群体选择一个个体i
        stra_player1=stra(player1);%得到个体i的策略
        Neig=[];
        for j=1:N
%             if Adj(player1,j)==1
            if Adj(player1,j) > 0
                Neig=[Neig,j];%记录个体player1的邻居            
            end            
        end
        Neig_Size=size(Neig,2);%个体自身的度
        Neig_rand0=randi(Neig_Size);
        Neig_rand=Neig( Neig_rand0);%随机选择一个邻居         
        player11=Neig_rand;
        Neig1=[];
        for j=1:N
%             if Adj(player11,j)==1
            if Adj(player11,j) > 0
                Neig1=[Neig1,j];%记录个体player11的邻居            
            end            
        end   
        Neig_Size1=size(Neig1,2);%邻居的度
        stra_neig=stra(Neig_rand);%得到个体的策略
        if stra_neig~=stra_player1
            score1=unity(player1);
            score2=unity(Neig_rand);        
            dscore=score2-score1;
            fermi=dscore/b/max(Neig_Size,Neig_Size1);
            
            if(fermi>rand(1)) 
                stra_player1=stra_neig;
            end
            stra(i)=stra_player1;
        end    
    end
    
end
%% 以上得到了各个时刻，所有个体的策略矩阵Stra；收益矩阵Unity
end




