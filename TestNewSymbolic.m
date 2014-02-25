% Script to obtain symbolic solution for 2-intron problem ...

load('AltSplicingRecruit_1-3_structures.mat')

Times = createSymData(4,'t');
Ks = createSymData(7,'k');

Q2 = makeQMatrixFromS_symbolic(spliceGraph(1).S(1:2,1),Ks(spliceGraph(1).VmoveIndex));
Theta2 = [1 0] * SimulateMarkov(Q2,Times(1))

Q3 = makeQMatrixFromS_symbolic(spliceGraph(1).S,Ks(spliceGraph(1).VmoveIndex));
Theta3 = [Theta2 0 0 0] * SimulateMarkov(Q3,Times(2))

Q4 = makeQMatrixFromS_symbolic(spliceGraph(2).S(1:10,1:15),Ks(spliceGraph(2).VmoveIndex));
Theta4 = [Theta3 0 0 0 0 0] * SimulateMarkov(Q4,Times(3))

Q5 = makeQMatrixFromS_symbolic(spliceGraph(2).S,Ks(spliceGraph(2).VmoveIndex));
%Theta5 = [Theta4 zeros(1,19)] * SimulateMarkov(Q5,Times(4))