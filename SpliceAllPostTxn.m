function [FinalConcentration,OutputTimes,OutputX, Q] = SpliceAllPostTxn(FinalConcentration,spliceGraph,Time,Ks,Threshold, verbose,Q)
% This function runs the ODE solver until all species have spliced.  

if (nargin < 5)
    Threshold = 0.001;
end

if (nargin < 6)
    verbose = true;
end


% Figure out how much incomplete splicing still remains
% I.e., find out the weight of every internal node (non-leaf node).
nonLeafNodes = unique(spliceGraph.VnodeIndex);
totalUnspliced = sum(FinalConcentration(nonLeafNodes));


global K SG;

K = Ks;
SG = spliceGraph;

if (nargin < 7)
     Q = makeQMatrixFromS(SG.S,K(SG.VmoveIndex));           
end


OutputTimes = {};
OutputX = {};

% Keep allowing it to splice over time until virtually 100% is spliced
while totalUnspliced >= Threshold
    if verbose
        disp('doing more until all spliced');
    end
    %[t,x] = ode23s(@Simulate3,[0 Time*10], FinalConcentration);
    
    % 5/15/13 Using Q-matrix
   ExtraTime = 10^10;
    Q_T = SimulateMarkov(Q,ExtraTime);
    x = FinalConcentration * Q_T;
       
    FinalConcentration = x(end,:);
    totalUnspliced = sum(FinalConcentration(nonLeafNodes));
    OutputX = [OutputX; {x}];
    OutputTimes = [OutputTimes; {ExtraTime}];
end


end