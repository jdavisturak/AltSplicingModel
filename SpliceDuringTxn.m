function [FinalConcentration,OutputTimes,OutputX, Q] = SpliceDuringTxn(NumIntrons,Times,Ks,spliceGraph, verbose)
% This function runs the ODE solver for the specified times,  which are calculated indepently.
% UPDATE: Now solve 'analyztically' using matrix exponentiation 

if NumIntrons > length(spliceGraph)
    fprintf('Sorry, not implemented for more than %d introns!',length(spliceGraph));
    return;
end


if (nargin < 5)
    verbose = true;
end


global K SG;

K = Ks;
StartingConcentration = [1 0];

OutputTimes = {};
OutputX = {};

for( N =1:NumIntrons)
    SG = spliceGraph(N);
    
    if verbose
        fprintf('Intron %d.  ', N);    tic
    end
    
    % Add new empty nodes
    StartingConcentration = [StartingConcentration zeros(1,size(SG.S,1)-length(StartingConcentration))];
    
    % New way! Analytical solution using matrix diagonalization.
    Q = makeQMatrixFromS(SG.S,K(SG.VmoveIndex));        
    Q_T = SimulateMarkov(Q,Times(N));
    t = Times(N);
    x = StartingConcentration * Q_T;
    

    %%% Solve the ODE's.
    %%[t,x] = ode23s(@Simulate3,[0 Times(N)], StartingConcentration);
    
    %     fprintf('Completely unspliced: %f\n', x(end,1)); toc
    
    % Save the ending concentration
    StartingConcentration = x(end,:);
    saveData{N} = {t,x};  % save all the data too.    
    
    OutputX = [OutputX; {x}];
    OutputTimes = [OutputTimes; {t}];
end
FinalConcentration = StartingConcentration;

end