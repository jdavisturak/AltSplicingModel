function [ output_args ] = SimulateCassetteExon_varRates_bareBones( KKs, T )
% SimulateCassetteExon_varRates_bareBones( KKs, T )
% This function simulates the model by solving the equations anew each time
% for all regions. (the non _bareBones version uses a precomputed solution
% for phases 2-4
%
% KKs: matrix of 5 *7 rate constants (7 constants, for each phase II through VI)
% T: vector of 4 phase durations (for phases II through V)


%prob2 = [1 0] *ouble(SimulateMarkov(makeQMatrixFromS_symbolic(spliceGraph(1).S(1:2,1),KKs(1,spliceGraph(1).VmoveIndex)),T(1)));
prob2 = [1 0] * SimulateMarkov(F_Q2(KKs(1,1)),T(1));

% prob3 = [prob2 0 0 0] * double(SimulateMarkov(makeQMatrixFromS_symbolic(spliceGraph(1).S,KKs(2,spliceGraph(1).VmoveIndex)),T(2)));
prob3 = [prob2 0 0 0] * SimulateMarkov(F_Q3(KKs(2,1),KKs(2,2),KKs(2,3)),T(2));

%prob4 = [prob3 0 0 0 0 0] * double(SimulateMarkov(makeQMatrixFromS_symbolic(spliceGraph(2).S(1:10,1:15),KKs(3,spliceGraph(2).VmoveIndex)),T(3)));
prob4 = [prob3 0 0 0 0 0] * SimulateMarkov(F_Q4(KKs(3,1),KKs(3,2),KKs(3,3),KKs(3,4)),T(3));

% using the pre-computed Q matrix here, for now... 
prob5 = [prob4 zeros(1,19)] * SimulateMarkov(F_Q5(KKs(4,1),KKs(4,2),KKs(4,3),KKs(4,4),KKs(4,5),KKs(4,6),KKs(4,7)),T(4));

prob6 = prob5  * SimulateMarkov(F_Q5(KKs(5,1),KKs(5,2),KKs(5,3),KKs(5,4),KKs(5,5),KKs(5,6),KKs(5,7)),1000*max(KKs(5,:)));

output_args = prob6(14);

end

