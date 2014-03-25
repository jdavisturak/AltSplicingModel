%% Create function for running Alt-splicing model of Cassette Exon
% Analytic solution for determining the hitting probability of including
% the exon given elongation times and kinetic rates of splicing events (co-
% and post-transcriptionally)


%% Load spliceGraph

%GenerateSmatrices
load('AltSplicingRecruit_1-3_structures.mat')

% spliceGraph is a struct of length 3 that contains information necessary
% to model genes of 1, 2 and 3 introns, respectively.
%
% Here we are interested in the set of events surrounding 2 introns


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

% This script simplifies the final phase(s) by calculating hitting
% probabilities of each node that exists at the beginning of phase V.

% These nodes have probability of 1 of including the middle exon (because
% they have already spliced at least one of the introns)
OneNodes  = {'1,2,3','1,2,3,4','1,2,3,4,5','1,2,3,4,5,7','1,2,4,5,7','1,4,5,7','1,2,3,5','2,4,5,7','4,5,7'};

% These nodes have probability of 0 of including the middle exon (they are
% all versions of the final form that is missing the internal exon, but they differ in the recruitment
% events that previously took place around the middle exon)
ZeroNodes = {'1,2,4,5,6','1,4,5,6','1,2,5,6','1,5,6'};

% These nodes have a chance of resulting in either outcome.
% These are written in the order that matrix A
SolutionNodes = {'1,2,4,5','1,2,5','1,4,5','1,5','1,2,4','2,4,5','2,4','2,5','4,5','1,4','1,2','1','2','4','5','t'};

% These next two lines re-order the output of this exercise into the order
% that spliceGraph(2) uses

AllNodes = {SolutionNodes{:} OneNodes{:} ZeroNodes{:} };
NewIndex = cellfun(@(node)find(strcmp(node,AllNodes)),spliceGraph(2).NodesTable);                 
% NewIndex = [16 12 13 11 17  7  5 18 14 10  3  1 19 20 26 21 27 22  2 23 28  6 24  9 25  8  4 29 15]'

%% Generate Q-matrices and simulation results (symbolic) for Phases II-IV:

syms k1 k2 k3 k4 k5 k6 k7
Ks = [ k1 k2 k3 k4 k5 k6 k7];
syms t2 t3 t4

% Phase II: First two states, first reaction are possible
Q2 = makeQMatrixFromS_symbolic(spliceGraph(1).S(1:2,1),Ks(spliceGraph(1).VmoveIndex));
Sim2 = SimulateMarkov(Q2,t2);

% Phase III: All states from 1-intron case are possible
Q3 = makeQMatrixFromS_symbolic(spliceGraph(1).S,Ks(spliceGraph(1).VmoveIndex));
Sim3 = SimulateMarkov(Q3,t3);

% Phase IV: First 10 states, 15 reactions are possible
Q4 = makeQMatrixFromS_symbolic(spliceGraph(2).S(1:10,1:15),Ks(spliceGraph(2).VmoveIndex));
Sim4 = SimulateMarkov(Q4,t4);

% Phase V: All states from 2-intron case are possible
Q5 = makeQMatrixFromS_symbolic(spliceGraph(2).S,Ks(spliceGraph(2).VmoveIndex));
% Sim5 = SimulateMarkov(Q5,t5); takes too long!


%% Turn the simulation results into functions:
F_Sim2 = matlabFunction(Sim2);
F_Sim3 = matlabFunction(Sim3);
F_Sim4 = matlabFunction(Sim4);
F_Q5   = matlabFunction(Q5);
F_Q4   = matlabFunction(Q4);
F_Q3   = matlabFunction(Q3);
F_Q2   = matlabFunction(Q2);

%% Special S-matrix substitution for Exon Definition:
%
% Exon definition is mathematically defined here as dependent recruitment rates.  
% In other words, in transcripts where reaction 2 (recruit 3'ss of intron 1) has occurred, reaction 4 (recruit 5'ss of intron 2) is faster, and vice-versa.
%

% Find the appropriate states ... 
node_has2 = find(cellfun(@(x)~isempty(x),regexp(spliceGraph(2).NodesTable,'2')));
node_has4 = find(cellfun(@(x)~isempty(x),regexp(spliceGraph(2).NodesTable,'4')));

% Find the fluxes that use 2/4:
which_flux_use2 = find(spliceGraph(2).VmoveIndex==2);
which_flux_use4 = find(spliceGraph(2).VmoveIndex==4);

% Now subset the fluxes (i.e. find flux where SOURCE node 'has2' or 'has4')
flux_use2_has4 = which_flux_use2(ismember(spliceGraph(2).VnodeIndex(which_flux_use2), node_has4));
flux_use4_has2 = which_flux_use4(ismember(spliceGraph(2).VnodeIndex(which_flux_use4), node_has2));

% RESULT: index into fluxes that can act non-linearly...
% Now make new versions of the Q4 an Q5 matrices, with 2 extra variables (enhancement rates of the new reactions)
% Here I am defining the enhancements as follows: 
%   -E1*k4 is the rate at which reaction 4 (recruit 5'ss of intron 2) occurs if 3'ss of intron1 is there already.
%   -E2*k2 is the rate at which reaction 2 (recruit 3'ss of intron 1) occurs if 5'ss of intron2 is there already.

syms E1 E2;
FluxVector = Ks(spliceGraph(2).VmoveIndex);
FluxVector(flux_use4_has2) = k4*E1;
FluxVector(flux_use2_has4) = k2*E2;


% Make new S/Q matrices for phases IV, V
Q4_2 = makeQMatrixFromS_symbolic(spliceGraph(2).S(1:10,1:15),FluxVector);
Q5_2 = makeQMatrixFromS_symbolic(spliceGraph(2).S,FluxVector);
F_Q5_2   = matlabFunction(Q5_2);
F_Q4_2   = matlabFunction(Q4_2);

    


%% Generate system of linear equations for solving hitting probabilities AFTER phase 5
A = sym(zeros(16));
B = sym(zeros(16,1));

% k1 is recruit 5' on intron 1
% k2 is recruit 3' on intron 1
% k3 is constitutively splice intron 1
% k4 is recruit 5' on intron 2
% k5 is recruit 3' on intron 2
% k6 is splice alternative form
% k7 is constitutively splice intron 2

% ==> Codes for names of nodes <==
% {1,2,4,5} means that raections 1,2,4,5 have occurred

%  1.   H{1,2,4,5} = (k3+k7)/(k3 + k7 + k6) %

A(1,1) = 1; 
B(1) = (k3+k7)/(k3 + k7 + k6);

%  2.   H{1,2,5} = (k3 + k4*H{1,2,4,5}) / (k3 + k4 + k6)
%       H{1,2,5} - k4*H{1,2,4,5}/ (k3 + k4 + k6) = k3/ (k3 + k4 + k6)

A(2,1:2) = [-k4/(k3 + k4 + k6) 1];
B(2) = k3/ (k3 + k4 + k6);

%  3.   H{1,4,5} - k2*H{1,2,4,5} / (k2 + k7 + k6 ) = k7/(k2 + k7 + k6)

A(3,[1 3]) = [-k2 / (k2 + k7 + k6) 1];
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

A(8, [ 2 6 8]) = [ -[k1 , k4]/(k1 + k4) 1];
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

% Result:
% A = [[                  1,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                       0,                       0,                       0,                       0, 0];
% [ -k4/(k3 + k4 + k6),                  1,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                       0,                       0,                       0,                       0, 0];
% [ -k2/(k2 + k6 + k7),                  0,                  1,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                       0,                       0,                       0,                       0, 0];
% [                  0, -k2/(k2 + k4 + k6), -k4/(k2 + k4 + k6),                  1,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                       0,                       0,                       0,                       0, 0];
% [      -k5/(k3 + k5),                  0,                  0,                  0,                  1,                  0,                  0,                  0,                  0,                  0,                  0,                       0,                       0,                       0,                       0, 0];
% [      -k1/(k1 + k7),                  0,                  0,                  0,                  0,                  1,                  0,                  0,                  0,                  0,                  0,                       0,                       0,                       0,                       0, 0];
% [                  0,                  0,                  0,                  0,      -k1/(k1 + k5),      -k5/(k1 + k5),                  1,                  0,                  0,                  0,                  0,                       0,                       0,                       0,                       0, 0];
% [                  0,      -k1/(k1 + k4),                  0,                  0,                  0,      -k4/(k1 + k4),                  0,                  1,                  0,                  0,                  0,                       0,                       0,                       0,                       0, 0];
% [                  0,                  0, -k1/(k1 + k2 + k7),                  0,                  0, -k2/(k1 + k2 + k7),                  0,                  0,                  1,                  0,                  0,                       0,                       0,                       0,                       0, 0];
% [                  0,                  0,      -k5/(k2 + k5),                  0,      -k2/(k2 + k5),                  0,                  0,                  0,                  0,                  1,                  0,                       0,                       0,                       0,                       0, 0];
% [                  0, -k5/(k3 + k4 + k5),                  0,                  0, -k4/(k3 + k4 + k5),                  0,                  0,                  0,                  0,                  0,                  1,                       0,                       0,                       0,                       0, 0];
% [                  0,                  0,                  0, -k5/(k2 + k4 + k5),                  0,                  0,                  0,                  0,                  0, -k4/(k2 + k4 + k5), -k2/(k2 + k4 + k5),                       1,                       0,                       0,                       0, 0];
% [                  0,                  0,                  0,                  0,                  0,                  0, -k4/(k1 + k4 + k5), -k5/(k1 + k4 + k5),                  0,                  0, -k1/(k1 + k4 + k5),                       0,                       1,                       0,                       0, 0];
% [                  0,                  0,                  0,                  0,                  0,                  0, -k2/(k1 + k2 + k5),                  0, -k5/(k1 + k2 + k5), -k1/(k1 + k2 + k5),                  0,                       0,                       0,                       1,                       0, 0];
% [                  0,                  0,                  0, -k1/(k1 + k2 + k4),                  0,                  0,                  0, -k2/(k1 + k2 + k4), -k4/(k1 + k2 + k4),                  0,                  0,                       0,                       0,                       0,                       1, 0];
% [                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  0, -k1/(k1 + k2 + k4 + k5), -k2/(k1 + k2 + k4 + k5), -k4/(k1 + k2 + k4 + k5), -k5/(k1 + k2 + k4 + k5), 1]]
%  
% 
% B = [(k3 + k7)/(k3 + k6 + k7);
%         k3/(k3 + k4 + k6);
%         k7/(k2 + k6 + k7);
%                         0;
%              k3/(k3 + k5);
%              k7/(k1 + k7);
%                         0;
%                         0;
%         k7/(k1 + k2 + k7);
%                         0;
%         k3/(k3 + k4 + k5);
%                         0;
%                         0;
%                         0;
%                         0;
%                         0;]

%% Substitute simpler components for the denominators

A2 = A;
A2=subs(A2,'k1+k2+k4+k5','G16');
A2=subs(A2,'k3 + k7 + k6','G1');
A2=subs(A2,'k3 + k4 + k6','G2');
A2=subs(A2,'k2 + k7 + k6','G3');
A2=subs(A2,'k2 + k4 + k6','G4');
A2=subs(A2,'k1 + k2 + k7','G9');
A2=subs(A2,'k4 + k5 + k3','G11');
A2=subs(A2,'k2+k4+k5','G12');
A2=subs(A2,'k1+k4+k5','G13');
A2=subs(A2,'k2+k1+k5','G14');
A2=subs(A2,'k2+k1+k4','G15');
A2=subs(A2,'k3 + k5','G5');
A2=subs(A2,'k1 + k7','G6');
A2=subs(A2,'k1 + k5','G7');
A2=subs(A2,'k1 + k4','G8');
A2=subs(A2,'k2 + k5','G10');

B2 = B;
B2=subs(B2,'k1+k2+k4+k5','G16');
B2=subs(B2,'k3 + k7 + k6','G1');
B2=subs(B2,'k3 + k4 + k6','G2');
B2=subs(B2,'k2 + k7 + k6','G3');
B2=subs(B2,'k2 + k4 + k6','G4');
B2=subs(B2,'k1 + k2 + k7','G9');
B2=subs(B2,'k4 + k5 + k3','G11');
B2=subs(B2,'k2+k4+k5','G12');
B2=subs(B2,'k1+k4+k5','G13');
B2=subs(B2,'k2+k1+k5','G14');
B2=subs(B2,'k2+k1+k4','G15');
B2=subs(B2,'k3 + k5','G5');
B2=subs(B2,'k1 + k7','G6');
B2=subs(B2,'k1 + k5','G7');
B2=subs(B2,'k1 + k4','G8');
B2=subs(B2,'k2 + k5','G10');

%% Solve the system of linear equations: A2*XX = B2
XX=linsolve(A2,B2);
F_XX = matlabFunction(XX); % turns the symbolic variable XX into a regular Matlab function

% Create a function that substitutes the parameters back into the function
F_XX_full = @(k1,k2,k3,k4,k5,k6,k7) F_XX(k3+k7+k6, k3+k4+k6, k2+k7+k6, k2+k4+k6, k3+k5, k1+k7, k1+k7, k1+k4, k1+k2+k7, k2+k5, k4+k5+k3, k2+k4+k5, k1+k4+k5, k2+k1+k5, k2+k1+k4,k1+k2+k4+k5,k1,k2,k3,k4,k5,k7);

%% The function F_XX_full(k1,k2,k3,k4,k5,k6,k7) returns the hitting probabilities of the 16 'SolutionStates'

save('CassetteExonModel.mat','F_XX_full','F_XX','NewIndex','OneNodes','ZeroNodes');

% Now make these functions into 'real' functions
makeFunction = @(handle, Name) regexprep(char(handle),'^@(\([\w,]*\))(.*)',sprintf('function output = %s$1\n\toutput=$2;\nend',Name));

f=fopen('F_Sim2.m','w');
fprintf(f,makeFunction(F_Sim2,'F_Sim2'));
fclose(f);

f=fopen('F_Sim3.m','w');
fprintf(f,makeFunction(F_Sim3,'F_Sim3'));
fclose(f);

f=fopen('F_Sim4.m','w');
fprintf(f,makeFunction(F_Sim4,'F_Sim4'));
fclose(f);

f=fopen('F_Q2.m','w');
fprintf(f,makeFunction(F_Q2,'F_Q2'));
fclose(f);

f=fopen('F_Q3.m','w');
fprintf(f,makeFunction(F_Q3,'F_Q3'));
fclose(f);

f=fopen('F_Q4.m','w');
fprintf(f,makeFunction(F_Q4,'F_Q4'));
fclose(f);

f=fopen('F_Q5.m','w');
fprintf(f,makeFunction(F_Q5,'F_Q5'));
fclose(f);

f=fopen('F_XX_full.m','w');
fprintf(f,makeFunction(F_XX_full,'F_XX_full'));
fclose(f);

f=fopen('F_XX.m','w');
fprintf(f,makeFunction(F_XX,'F_XX'));
fclose(f);

f=fopen('F_Q4_exonDefinition.m','w');
fprintf(f,makeFunction(F_Q4_2,'F_Q4_2'));
fclose(f);

f=fopen('F_Q5_exonDefinition.m','w');
fprintf(f,makeFunction(F_Q5_2,'F_Q5_2'));
fclose(f);


%% test
% 
% N=1000;KKs = [1;1.0001; 1.000056; 1; 1.001; 1.0005;1.0007]; T=[1/3 1/30 1/3 1];
% 
% disp('simulations')
% 
% tic
% for i=1:N
%    SimulateCassetteExon( KKs, T );
% end
% toc









