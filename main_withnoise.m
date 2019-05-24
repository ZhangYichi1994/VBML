clear all;clc;close all;
%% %%%%%%%%%%%%%%%%%% Main Program  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the Original Network Connection List Data, and Transfer the Connection List to Adjacent Matrix
% Input: Connection List
% Output: Adjacent Matrix
Cele1 = load('sf-25-1.txt');%Import the network data, Connection List
% Cele1 = load('sw-1.txt');%Import the network data, Connection List

SIZE = max(max(Cele1(:,1)),max(Cele1(:,2)));
Adj = zeros(SIZE,SIZE);% Get the Size of Adjacent Matrix
% Transfer the Connection List to Adjacent Matrix
for i=1:size(Cele1)
    a=Cele1(i,1);
    b=Cele1(i,2);
    Adj(a,b)=Cele1(i,3);
    Adj(b,a)=Cele1(i,3);   %zhang karate dataset
end
% diag(Adj)=1, self-interaction
for i=1:SIZE
     Adj(i,i)=1;
end

%% set parameters
b=1.2; % Define the game parameter
amp = 0.09;   % Define the amplifier of noise;
Length=1; % 数据量 


%% Get the Evolutionary Game Data, Including the Strategies and Unity
%     for amp = 0:0.05:1.2
%%  clustering G(or A), to explore the group porperty of X
% stra 为策略，列向量
% Unity 为收益，列向量   
[Stra,Unity] = Game_withnoise(Adj,b,amp);
%% Reconstruct the Network Structure Based on the Evolutionary Game Data

%% lasso 
Index = 0;
Adj_Re_BOMP = Net_Construction(Stra,Unity,Index,Length,b);
accuracy_lasso = EvaluationFcn(Adj_Re_BOMP,Adj);
adj=reshape(Adj,SIZE*SIZE,1);
adj_Re=reshape(Adj_Re_BOMP,SIZE*SIZE,1);
% prec_rec(adj_Re,adj,'holdFigure', 1); %ROC and PR
[AUROC_lasso, AUPR_lasso, prec3, tpr3, fpr3, thresh3]=prec_rec(adj_Re,adj,'holdFigure', 2); %ROC and PR

%% clusterVB
Index = 1;
Adj_Re_clusster = Net_Construction(Stra,Unity,Index,Length,b);
accuracy_clusterVB = EvaluationFcn(Adj_Re_clusster,Adj);
adj=reshape(Adj,SIZE*SIZE,1);
adj_Re=reshape(Adj_Re_clusster,SIZE*SIZE,1);
% prec_rec(adj_Re,adj,'holdFigure', 1); %ROC and PR
[AUROC_clusster ,AUPR_clusster, prec3, tpr3, fpr3, thresh3]=prec_rec(adj_Re,adj,'holdFigure', 2); %ROC and PR
