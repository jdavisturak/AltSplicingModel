%% Initial script for testing alt splicing model
% 6/25/2012


%% Test ODE
global S VnodeIndex VnodeMultiplier;
for( N =1:9)
    load(sprintf('AltSplicingMatrices%2d.mat',N))
    
    fprintf('%d Exons.  ', N);    tic
    
    X = zeros(NumNodes,1);
    X(1) = 1; % Start with everything in the first state
    
    [t,x] = ode23s(@Simulate,[0 .1],X);
    
    fprintf('Completely unspliced: %f\n', x(end,1)); toc
end

%% Next step:
% Include rate parameters (A vector: one for each of the N*(N+1)/2 moves)
%     global K;
% Step-wise splicing of individual regions:
%   Run the ODE for each phase, 1:N (N is # of introns in the gene).
%       After each phase, the final concentrations are the STARTING concentrations of the next
%       phase.

%% Assuming there is some vector K , Times
clear
StartingConcentration = [1 0];

global S VnodeIndex K VmoveIndex;
K = 2*ones(45,1);
Times = ones(10,1);

load('AltSplicing_1-9_structures.mat')
saveData = {};
for( N =1:7)
    S = S_cell{N};
    VnodeIndex = VnodeIndex_cell{N};
    VmoveIndex = VmoveIndex_cell{N};
    
    fprintf('Intron %d.  ', N);    tic
    
    % Add new empty nodes
    StartingConcentration = [StartingConcentration zeros(1,size(S,1)-length(StartingConcentration))];
    
    % Solve the ODE's.
    [t,x] = ode23s(@Simulate2,[0 Times(N)], StartingConcentration);
    
    fprintf('Completely unspliced: %f\n', x(end,1)); toc
    
    % Save the ending concentration
    StartingConcentration = x(end,:);
    saveData{N} = {t,x};  % save all the data too.    
end
FinalConcentration = StartingConcentration;

%% Figure out how much incomplete splicing still remains
% I.e., find out the weight of every internal node (non-leaf node).
nonLeafNodes = unique(VnodeIndex_cell{N});
totalUnspliced = sum(FinalConcentration(nonLeafNodes));

while totalUnspliced >= 0.01
    disp('doing more until all spliced');
    [t,x] = ode23s(@Simulate2,[0 Times(N)], FinalConcentration);
    
    FinalConcentration = x(end,:);
    totalUnspliced = sum(FinalConcentration(nonLeafNodes))

end


%% Test function versions
clear
load('AltSplicing_1-9_structures.mat'); % load spliceGraph
N=5;
[FinalConcentration,OutputTimes,OutputX] = SpliceDuringTxn(N,ones(N,1),0.1*ones(N*(N+1)/2,1),spliceGraph);
[a b c] = SpliceAllPostTxn(FinalConcentration,spliceGraph(N),1,0.0001);

%% Test speed
clear
load('AltSplicing_1-9_structures.mat'); % load spliceGraph
a=[]
for N=3:8;
    for(i =1:5)
        tic ;
        [FinalConcentration,OutputTimes,OutputX] = SpliceDuringTxn(N,5*ones(N,1),0.1*ones(N*(N+1)/2,1),spliceGraph);
        a=[a toc]
    end
end
times = reshape(a,5,numel(a)/5)

%


%%
a=ones(gene.exonEnds(end),1);
Ks = [10 1 10 1 1 10];
% Ks = [1 1 1 1 1 1];
Ke = a*30000;    Kp = a*60;  Ku = a*30;
[FinalConcentration,OutputTimes,OutputX , FinalConcentration2,OutputTimes2,OutputX2,tx,complete ]=runFullModel_max9(gene,Ke,Kp,Ku,Ks);
FinalConcentration2(complete)


%% 7/16/12
% Test CTS vs PTS ...
clear
load('AltSplicing_1-9_structures.mat'); % load spliceGraph
N=2;
times = [1 1];
TranscriptNames = spliceGraph(N).TranscriptNames;
whichCompleteSplicers = completeSplicers(spliceGraph(N))

% CTS + PTS
K =ones(N*(N+1)/2,1)
[FinalConcentration,OutputTimes,OutputX] = SpliceDuringTxn(N,times,K,spliceGraph);
[FinalCTS OutputTimes_PTS1 OutputX_PTS1] = SpliceAllPostTxn(FinalConcentration,spliceGraph(N),times(1),K);
FinalCTS(whichCompleteSplicers)

figure(1)
plot(OutputTimes{1},OutputX{1});ylim([0 1]);
hold on;
plot(OutputTimes{2}+times(1),OutputX{2});ylim([0 1]);
hold on;
plot(OutputTimes_PTS1{1}+times(1)+times(2),OutputX_PTS1{1});ylim([0 1]);xlim([0 10]);
hold off
legend(TranscriptNames);
CTS_times = [OutputTimes{1}; OutputTimes{2}+times(1);OutputTimes_PTS1{1}+times(1)+times(2)];

% PTS only: tune parameters to get the same % splicing
K2 = [1; .2793; 1;]; % .2793
[FinalPTS OutputTimes_PTS2 OutputX_PTS2] = SpliceAllPostTxn([1 0 0 0 0],spliceGraph(N),times(1),K2);
[FinalPTS(whichCompleteSplicers);FinalCTS(whichCompleteSplicers)]

figure(2)
plot(OutputTimes_PTS2{1},OutputX_PTS2{1})
legend(TranscriptNames);

%% Plot INTRON abundances
I1_index = [1 4];
E2_index1 = [1 2];
E2_index2 = [1 2 4 5];
I2_index = [1 2];
E3_index = [1:5];

regions = {'01i','02e','02i','03e'};
PTS_Regions = [sum(OutputX_PTS2{1}(:,I1_index),2) sum(OutputX_PTS2{1}(:,E2_index2),2) sum(OutputX_PTS2{1}(:,I2_index),2) sum(OutputX_PTS2{1}(:,E3_index),2)];
% PTS txn:
PTS_txn_t = [0 1 2]';
PTS_txn_x = [1 0 0 0; 1 1 1 0; 1 1 1 1];

% Try CTS where I only let the region count if it is added AFTER the splice period
CTS_v1_Regions_1 = [sum(OutputX{1}(:,1),2) zeros(size(OutputX{1},1),3)];
CTS_v1_Regions_2 = [sum(OutputX{2}(:,I1_index),2) sum(OutputX{2}(:,E2_index1),2) sum(OutputX{2}(:,I2_index),2) 0.*sum(OutputX{2}(:,E3_index),2)];
CTS_v1_Regions_3 = [sum(OutputX_PTS1{1}(:,I1_index),2) sum(OutputX_PTS1{1}(:,E2_index2),2) sum(OutputX_PTS1{1}(:,I2_index),2) sum(OutputX_PTS1{1}(:,E3_index),2)];
CTS_v1_Regions = [CTS_v1_Regions_1; CTS_v1_Regions_2;CTS_v1_Regions_3];

% CTS where I add the regions at the beginning
CTS_v2_Regions_1 = [sum(OutputX{1}(:,1),2) sum(OutputX{1}(:,E2_index1),2) zeros(size(OutputX{1},1),2)];
CTS_v2_Regions_2 = [sum(OutputX{2}(:,I1_index),2) sum(OutputX{2}(:,E2_index2),2) sum(OutputX{2}(:,I2_index),2) sum(OutputX{2}(:,E3_index),2)];
CTS_v2_Regions_3 = [sum(OutputX_PTS1{1}(:,I1_index),2) sum(OutputX_PTS1{1}(:,E2_index2),2) sum(OutputX_PTS1{1}(:,I2_index),2) sum(OutputX_PTS1{1}(:,E3_index),2)];

CTS_v2_Regions = [CTS_v2_Regions_1; CTS_v2_Regions_2;CTS_v2_Regions_3];

%CTS_v2_Regions_2 = [sum(OutputX{2}(:,I1_index),2) sum(OutputX{2}(:,E2_index2),2) sum(OutputX{2}(:,I2_index),2) sum(OutputX{2}(:,E3_index),2)];

figure(3)
subplot(1,3,1);plot([PTS_txn_t; OutputTimes_PTS2{1}+2],[PTS_txn_x;PTS_Regions]); legend(regions);
subplot(1,3,2);plot(CTS_times     ,  CTS_v1_Regions); legend(regions);
subplot(1,3,3);plot(CTS_times     ,  CTS_v2_Regions); legend(regions);


%% 7/17/12
% Show CTS has more high-fidelity splicing
clear
load('AltSplicing_1-9_structures.mat'); % load spliceGraph
N=5;
times = ones(1,N)  / 3;
TranscriptNames = spliceGraph(N).TranscriptNames;
whichCompleteSplicers = completeSplicers(spliceGraph(N))

% CTS + PTS
K = ones(N*(N+1)/2,1);
[FinalConcentration,OutputTimes,OutputX] = SpliceDuringTxn(N,times,K,spliceGraph);
[FinalCTS OutputTimes_PTS1 OutputX_PTS1] = SpliceAllPostTxn(FinalConcentration,spliceGraph(N),times(1),K);
FinalCTS(whichCompleteSplicers(1));

[FinalPTS OutputTimes_PTS2 OutputX_PTS2] = SpliceAllPostTxn([1 zeros(1,spliceGraph(N).NumNodes-1)],spliceGraph(N),times(1),K);
[FinalPTS(whichCompleteSplicers(1));
FinalCTS(whichCompleteSplicers(1))]



























