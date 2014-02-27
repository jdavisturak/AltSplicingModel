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

%% Create a system of linear equations for hitting probabilities of all 16 non-spliced states in the 2-intron case. 
% H(x) = P(reach 'include' isoform | state i = x)
Ks = createSymData(7,'k');
A = sym(zeros(16));
B = sym(zeros(16,1));
syms k1 k2 k3 k4 k5 k6 k7

% k1 is recruit 5' on intron 1
% k2 is recruit 3' on intron 1
% k3 is constitutively splice intron 1
% k4 is recruit 5' on intron 2
% k5 is recruit 3' on intron 2
% k6 is splice alternative form
% k7 is constitutively splice intron 2


%  1.   H{1,2,4,5} = (k3+k7)/(k3 + k7 + k6)

A(1,1) = 1; 
B(1) = (k3+k7)/(k3 + k7 + k6);

%  2.   H{1,2,5} = (k3 + k4*H{1,2,4,5}) / (k3 + k4 + k6)
%       H{1,2,5} - k4*H{1,2,4,5}/ (k3 + k4 + k6) = k3/ (k3 + k4 + k6)

A(2,1:2) = [-k4/(k3 + k4 + k6) 1];
B(2) = k3/ (k3 + k4 + k6);

%  3.   H{1,4,5} - k4*H{1,2,4,5} / (k2 + k7 + k6 ) = k7/(k2 + k7 + k6)

A(3,[1 3]) = [-k4 / (k2 + k7 + k6) 1];
B(3) = k7/(k2 + k7 + k6);

%  4.   H{1,5} - (k2*H{1,2,5} + k4 *H{1,4,5}) / (k2 + k4 + k6) = 0

A(4,2:4) = [ -k2/(k2 + k4 + k6), -k4/(k2 + k4 + k6),  1];
B(4) = 0;

%  5.   H{1,2,4} - k5* (H{1,2,4,5}) / (k3 + k5) = k3/(k3 + k5)

A(5,[1 5]) = [-k5/(k3 + k5) 1];
B(5) = k3/(k3 + k5);

%  6.   H{2,4,5} - k1 * H{1,2,4,5} / (k1 + k7) = k7/(k1 + k7)

A(6,[1 6]) = [-k1/(k1 + k7) 1];
B(6) = k7/(k1 + k7);

%  7.   H{2,4} - (k1 * H{1,2,4} + k5 *H{2,4,5})/(k1 + k5)  = 0

A(7, [ 5 6 7]) = [ -[k1 , k5]/(k1 + k5) 1];
B(7) = 0;

%  8.   H{2,5} - (k1 * H{1,2,5} + k4 *H{2,4,5})/(k1 + k4)  = 0

A(8, [ 4 6 8]) = [ -[k1 , k4]/(k1 + k4) 1];
B(8) = 0;

%  9.   H{4,5} - (k1*H{1,4,5} + k2 *H{2.4.5})/ (k1 + k2 + k7) = k7/(k1 + k2 + k7)

A(9,[3 6 9]) = [ -[k1 , k2]/(k1 + k2 + k7) 1];
B(9) = k7/(k1 + k2 + k7);

% 10.   H{1,4} - (k2 * H{1,2,4} + k5 *H{1,4,5})/(k2 + k5)  = 0

A(10, [ 5 3 10]) = [ -[k2 , k5]/(k2 + k5) 1];
B(10) = 0;

% 11.   H{1,2} - (k4*H{1,2,4} + k5 *H{1,2,5})/ (k4 + k5 + k3) = k3/(k4 + k5 + k3)

A(11,[5 2 11]) = [ -[k4 , k5]/(k4 + k5 + k3) 1];
B(11) = k3/(k4 + k5 + k3);


% 12.   H{1) - (k2*H{1,2} + k4*H{1,4} + k5*H{1,5})/(k2 + k4 + k5) = 0

A(12,[11 10 4 12]) = [-[k2, k4, k5]/(k2+k4+k5)   1 ];
B(12) = 0;

% 13.   H{2) - (k1*H{1,2} + k4*H{2,4} + k5*H{2,5})/(k1 + k4 + k5) = 0

A(13,[11 7 8 13]) = [-[k1, k4, k5]/(k1+k4+k5)   1 ];
B(13) = 0;

% 14.   H{4) - (k1*H{1,4} + k2*H{2,4} + k5*H{4,5})/(k2 + k1 + k5) = 0

A(14,[10 7 9 14]) = [-[k1, k2, k5]/(k2+k1+k5)   1 ];
B(14) = 0;

% 15.   H{5) - (k1*H{1,5} + k2*H{2,5} + k4*H{4,5})/(k2 + k1 + k4) = 0

A(15,[4 8 9 15]) = [-[k1, k2, k4]/(k2+k1+k4)   1 ];
B(15) = 0;


% 16.   H{0) - (h1 + H{1} + k2*H{2} + k4*H{4} + k5*H{5})/(k1 + k2 + k4 + k5) = 0

A(16,12:16) = [-[k1 k2, k4, k5]/(k1+k2+k4+k5)   1 ];
B(16) = 0;



X=linsolve(A,B);


%% Now take the solutions of X and put them in the correct order to use with Q5 ...
%[spliceGraph(2).NodesTable]
%     't'
%     '1'
%     '2'
%     '1,2'
%     '1,2,3'
%     '2,4'
%     '1,2,4'
%     '1,2,3,4'
%     '4'
%     '1,4'
%     '1,4,5'
%     '1,2,4,5'
%     '1,2,3,4,5'
%     '1,2,3,4,5,7'
%     '1,2,4,5,6'
%     '1,2,4,5,7'
%     '1,4,5,6'
%     '1,4,5,7'
%     '1,2,5'
%     '1,2,3,5'
%     '1,2,5,6'
%     '2,4,5'
%     '2,4,5,7'
%     '4,5'
%     '4,5,7'
%     '2,5'
%     '1,5'
%     '1,5,6'
%     '5'

SolutionNodes = {'1,2,4,5','1,2,5','1,4,5','1,5','1,2,4','2,4,5','2,4','2,5','4,5','1,4','1,2','1','2','4','5','t'};
OneNodes  = {'1,2,3','1,2,3,4','1,2,3,4,5','1,2,3,4,5,7','1,2,4,5,7','1,4,5,7','1,2,3,5','2,4,5,7','4,5,7'};
ZeroNodes = {'1,2,4,5,6','1,4,5,6','1,2,5,6','1,5,6'};

AllNodes = {SolutionNodes{:} OneNodes{:} ZeroNodes{:} };
AllProbs = [X; ones(length(OneNodes),1); zeros(length(ZeroNodes),1)];

% Re-arrange

NewIndex = cellfun(@(node)find(strcmp(node,AllNodes)),spliceGraph(2).NodesTable);

FullProbs = AllProbs(NewIndex);

Theta5 = [Theta4 zeros(1,19)] * FullProbs;















