function [FinalConcentration,OutputTimes,OutputX] = SpliceAllPostTxn(FinalConcentration,spliceGraph,Time,Ks,Threshold)
% This function runs the ODE solver until all species have spliced.  

if (nargin < 5)
    Threshold = 0.001;
end


% Figure out how much incomplete splicing still remains
% I.e., find out the weight of every internal node (non-leaf node).
nonLeafNodes = unique(spliceGraph.VnodeIndex);
totalUnspliced = sum(FinalConcentration(nonLeafNodes));


global K SG;

K = Ks;
SG = spliceGraph;

OutputTimes = {};
OutputX = {};

% Keep allowing it to splice over time until virtually 100% is spliced
while totalUnspliced >= Threshold
    disp('doing more until all spliced');
    %[t,x] = ode23s(@Simulate3,[0 Time*10], FinalConcentration);
    
    % 5/15/13 Using Q-matrix
    Q = makeQMatrixFromS(SG.S,K(SG.VmoveIndex));        
    ExtraTime = 10^10;
    Q_T = SimulateMarkov(Q,ExtraTime);
    x = FinalConcentration * Q_T;
       
    FinalConcentration = x(end,:);
    totalUnspliced = sum(FinalConcentration(nonLeafNodes));
    OutputX = [OutputX; {x}];
    OutputTimes = [OutputTimes; {ExtraTime}];
end


end