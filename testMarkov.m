%% Test alt splice using Markov chain
% 5/12/2013

% Testing just one matrix at a time
load('AltSplicing_1-9_structures.mat');

global K SG;
N=8; 
SG = spliceGraph(N);
K = rand(length(SG.Moves),1);
startConc = [1 zeros(1,spliceGraph(N).NumNodes-1)];

% Method 1: analytical solution with Q matrix
Q = makeQMatrixFromS(spliceGraph(N).S,K(spliceGraph(N).VmoveIndex));    
% Method 2.1: analytical with symbolic: 
% tic
% Qs = makeQMatrixFromS_symbolic(spliceGraph(N).S, spliceGraph(N).VmoveIndex, length(SG.Moves));
% toc

T= 1;
tic
Q_T = SimulateMarkov(Q,T);
toc

% Method 2: ODE
tic
[t,x] = ode23s(@Simulate3,[0 T], startConc);
toc

% [startConc*Q_T; x(end,:); ]

% N = 6;
%Elapsed time is 0.025248 seconds.
%Elapsed time is 1.223179 seconds

% N = 7;
%Elapsed time is 0.355612 seconds.
%Elapsed time is 9.304419 seconds.

% For N = 8:
%Elapsed time is 10.114014 seconds.
% Elapsed time is 108.490873 seconds.

% For N = 9:
% 125 seconds
%


clear
load('AltSplicing_1-9_structures.mat'); % load spliceGraph
N=2;
global K SG;

SG = spliceGraph(N);
syms k1 k2 k3 t1 t2 t3
times = [ t1 t2 t3];
K = [k1 k2 k3]
TranscriptNames = SG.TranscriptNames;
whichCompleteSplicers = completeSplicers(SG)


startConc = [1 zeros(1,spliceGraph(N).NumNodes-1)];

% Method 1: analytical solution with Q matrix
Q = makeQMatrixFromS_symbolic(SG.S,SG.VmoveIndex,length(SG.Moves))    
Q_T = SimulateMarkov(Q,t1)



