% Script to obtain symbolic solution for 2-intron problem ...

load('AltSplicingRecruit_1-3_structures.mat')

Times = createSymData(4,'t');
Ks = createSymData(7,'k');

Q2 = makeQMatrixFromS_symbolic(spliceGraph(1).S(1:2,1),Ks(spliceGraph(1).VmoveIndex));
Theta2 = [1 0] * SimulateMarkov(Q2,Times(1));

Q3 = makeQMatrixFromS_symbolic(spliceGraph(1).S,Ks(spliceGraph(1).VmoveIndex));
Theta3 = [Theta2 0 0 0] * SimulateMarkov(Q3,Times(2));

Q4 = makeQMatrixFromS_symbolic(spliceGraph(2).S(1:10,1:15),Ks(spliceGraph(2).VmoveIndex));
Theta4 = [Theta3 0 0 0 0 0] * SimulateMarkov(Q4,Times(3));

Q5 = makeQMatrixFromS_symbolic(spliceGraph(2).S,Ks(spliceGraph(2).VmoveIndex));
%Theta5 = [Theta4 zeros(1,19)] * SimulateMarkov(Q5,Times(4))

%% Create a system of linear equations for hitting probabilities of all 16 non-at-all-spliced states in the 2-intron case. 
% H(x) = P(reach 'include' isoform | state i = x)
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

% ==> Codes for names of 'nodes' <==
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

disp('Solving linear equations ... ')

X=linsolve(A,B);

%% Try making some simpler sub-components and then doing solving, substituting

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

XX=linsolve(A2,B2);
F_XX = matlabFunction(XX);

% This next function must substitute from k's into G's, then do the X
% thingy
F_XX_full = @(k1,k2,k3,k4,k5,k6,k7) F_XX(k3 + k7 + k6, k3 + k4 + k6, k2 + k7 + k6, k2 + k4 + k6, k3 + k5, k1 + k7, k1 + k7, k1 + k4, k1 + k2 + k7, k2 + k5, k4 + k5 + k3, k2+k4+k5, k1+k4+k5, k2+k1+k5, k2+k1+k4,k1+k2+k4+k5,k1,k2,k3,k4,k5,k7);


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

% The nodes we just created X for are here:
SolutionNodes = {'1,2,4,5','1,2,5','1,4,5','1,5','1,2,4','2,4,5','2,4','2,5','4,5','1,4','1,2','1','2','4','5','t'};

% These nodes have probability of 1 of including the middle exon
OneNodes  = {'1,2,3','1,2,3,4','1,2,3,4,5','1,2,3,4,5,7','1,2,4,5,7','1,4,5,7','1,2,3,5','2,4,5,7','4,5,7'};

% These nodes have probability of 0 of including the middle exon (they are
% all versions of the spliced-out form, but they differ in the recruitment
% events that previously took place around the middle exon)
ZeroNodes = {'1,2,4,5,6','1,4,5,6','1,2,5,6','1,5,6'};

AllNodes = {SolutionNodes{:} OneNodes{:} ZeroNodes{:} };
NewIndex = cellfun(@(node)find(strcmp(node,AllNodes)),spliceGraph(2).NodesTable);


AllProbs = [X; ones(length(OneNodes),1); zeros(length(ZeroNodes),1)];

% Re-arrange

FullProbs = AllProbs(NewIndex);




%% Calculate final theta (PSI):
disp('Calculating dot product ... ')

Theta5 = [Theta4 zeros(1,19)] * FullProbs;




%% Make a (much less complicated version) where all reactions are dictated by only 3 kinetic parameters: 5' recruit, 3' recruit, and splice reaction.
syms K5 K3 Ksplice

Theta4_2 = simplify(subs(subs(subs(subs(subs(subs(subs(Theta4,k1,K5),k2,K3),k3,Ksplice),k4,K5),k5,K3),k6,Ksplice),k7,Ksplice));
X2 = simplify(subs(subs(subs(subs(subs(subs(subs(X,k1,K5),k2,K3),k3,Ksplice),k4,K5),k5,K3),k6,Ksplice),k7,Ksplice));
AllProbs_2 = [X2; ones(length(OneNodes),1); zeros(length(ZeroNodes),1)];
Theta5_2 = [Theta4_2 zeros(1,19)] * AllProbs_2(NewIndex);
CHAR2 = char(Theta5_2);
length(CHAR2)  % 9281 characters!!!

%% Convert to matlab function:

F_Region2 = matlabFunction(Theta2);
F_Region3 = matlabFunction(SimulateMarkov(Q3,Times(2)));
F_Region4 = matlabFunction(SimulateMarkov(Q4,Times(3)));
F_X = matlabFunction(X);
%F_Region5 = matlabFunction(SimulateMarkov(Q5,Times(4)));

F_Theta5_2 = matlabFunction(Theta5_2);
F_Theta4_2 = matlabFunction(Theta4_2);
F_X2 = matlabFunction(X2);
F_X = matlabFunction(X);

X_simp=simplify(X)
CH_simp=char(X_simp)
length(CH_X)

F_AllProbs_2 = matlabFunction(AllProbs_2(NewIndex));


save('Analytic_Theta5.mat','Theta5','Theta4_2','Q5','X','X2', 'Theta5_2','CHAR2')
save('Analytic_Theta5_2.mat','Theta5_2')
save('Analytic_Function_Theta5_2.mat','F_Theta5_2')
save('Analytic_Full_Functions.mat','F_Region2','F_Region3','F_Region4','F_X','NewIndex','OneNodes','ZeroNodes');

k1=1;k2=1.0001; k3=1.000056;k4=k1;k5=k2;k6=k3;k7=k3;
T=[1/3 1/30 1/3];


Ks = [ k1 k2 k3 k4 k5 k6 k7];
N=300;
tic
for i=1:N
    F_Theta5_2(k1, k2, k3,T(1),T(2),T(3)); % .7211
end
toc % .012 seconds

%% Comparing methods that use all 7 parameters:
tic
for i=1:N
    X_ans = [F_X(k1,k2,k3,k4,k5,k6,k7); ones(length(OneNodes),1); zeros(length(ZeroNodes),1)];
    [[[F_Region2(k1,T(1)) zeros(1,3)] * F_Region3(k1, k2, k3, T(2)), zeros(1,5)] * F_Region4(k1,k2,k3,k4,T(3)) zeros(1,19)] * X_ans(NewIndex); % .7211
end
toc % 1.5 seconds

tic
for i=1:N
    T2 = [1 0] * SimulateMarkov(makeQMatrixFromS(spliceGraph(1).S(1:2,1),Ks(spliceGraph(1).VmoveIndex)),T(1));
    
    T3 = [T2 0 0 0] * SimulateMarkov(makeQMatrixFromS(spliceGraph(1).S,Ks(spliceGraph(1).VmoveIndex)),T(2));
    
    T4 = [T3 0 0 0 0 0] * SimulateMarkov(makeQMatrixFromS(spliceGraph(2).S(1:10,1:15),Ks(spliceGraph(2).VmoveIndex)),T(3));
    
    T5 = [T4 zeros(1,19)] * SimulateMarkov(makeQMatrixFromS(spliceGraph(2).S,Ks(spliceGraph(2).VmoveIndex)),T(3));
    X_ans = [F_XX_full(k1,k2,k3,k4,k5,k6,k7); ones(length(OneNodes),1); zeros(length(ZeroNodes),1)];
    T5*X_ans(NewIndex);
end
toc % 1.2 seconds


for i=1:N
    T2 = [1 0] * Sim2(k1,T(1));
    
    T3 = [T2 0 0 0] * Sim3(k1,k2,k3,T(2));
    
    T4 = [T3 0 0 0 0 0] * SimulateMarkov(makeQMatrixFromS(spliceGraph(2).S(1:10,1:15),Ks(spliceGraph(2).VmoveIndex)),T(3));
    
    T5 = [T4 zeros(1,19)] * SimulateMarkov(makeQMatrixFromS(spliceGraph(2).S,Ks(spliceGraph(2).VmoveIndex)),T(3));
    X_ans = [F_XX_full(k1,k2,k3,k4,k5,k6,k7); ones(length(OneNodes),1); zeros(length(ZeroNodes),1)];
    T5*X_ans(NewIndex);
end
toc % 1.2 seconds
















