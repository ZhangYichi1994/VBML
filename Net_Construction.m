%% Network structure construction based on the strategies and unity data
% Input: Stra--strategies data; Unity--unity data
% Input:  Lamda1--Coefficient of L-inf; Lamda2--Coefficient of L2
% Input:  Index--0-5
% Input:  Length--the length of data, Length*SIZE
% Input:  b--coefficient of game model 
% Output: Adj--Reconstructed Adjacent Matrix; 
function Adj=Net_Construction(Stra,Unity,Index,Length,b)
    SIZE = size(Stra,1);
    Adj = zeros(SIZE,SIZE);
%      KNN = Commu_Maxt(SIZE,SIZE);        
    %% 计算U=Fai*x
    Length1 = floor(Length*SIZE);
    Unity_temp = Unity(:,1:Length1);
    Unity_use=reshape(Unity_temp', Length1*SIZE,1);% 得到 U
    Stra_temp = Stra(:,1:Length1);
    
    Stra_use=zeros(Length1*SIZE,SIZE*SIZE);    
    for k=1:SIZE
        TEMP=zeros(Length1,SIZE);
        player1=k;
        for t = 1:Length1 
            stra_player1=Stra_temp(player1,t);
            for i=1:SIZE
                player2=i;
                stra_player2=Stra_temp(player2,t);
                TEMP(t,i)=Payoff(stra_player1,stra_player2,b);
            end
        end
        aa=(k-1)*Length1;
        bb=(k-1)*SIZE;
        for row = 1:Length1
            for col = 1:SIZE                
                Stra_use(aa+row,bb+col)= TEMP(row,col); % 得到 Fai
            end
        end       
    end   

    %% lasso
    if (Index==0)
        for i = 1:SIZE
            Stra_use1 = Stra_use((i-1)*Length1+1:i*Length1,(i-1)*SIZE+1:i*SIZE);
            Unity_use1 = Unity_use((i-1)*Length1+1:i*Length1);
            [B,info] = lasso(Stra_use1,Unity_use1,'CV',10);
            Adj(:,i) = B(:,info.Index1SE);
        end       
    end
   
    %% VBML
    if(Index == 1)    
        for i = 1:SIZE
            Stra_use1 = Stra_use((i-1)*Length1+1:i*Length1,(i-1)*SIZE+1:i*SIZE);
            Unity_use1 = Unity_use((i-1)*Length1+1:i*Length1);
            [xhat, errbar,alpha0] = VBML(Stra_use1, Unity_use1, SIZE);
           Adj(:,i) = xhat;
        end           
    end     
end
