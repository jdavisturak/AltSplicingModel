function [ output_args ] = SimulateCassetteExon_varRates( KKs, T )
% KKs: matrix of 5 *7 rate constants (7 constants, for each phase II through VI)
% T: vector of 4 phase durations (for phases II through V)

% NewIndex = [16 12 13 11 17  7  5 18 14 10  3  1 19 20 26 21 27 22  2 23 28  6 24  9 25  8  4 29 15];

prob2 = [1 0] * F_Sim2(KKs(1,1),T(1));

prob3 = [prob2 0 0 0] * F_Sim3(KKs(2,1),KKs(2,2),KKs(2,3),T(2));

prob4 = [prob3 0 0 0 0 0] * F_Sim4(KKs(3,1),KKs(3,2),KKs(3,3),KKs(3,4),T(3));

prob5 = [prob4 zeros(1,19)] * SimulateMarkov(F_Q5(KKs(4,1),KKs(4,2),KKs(4,3),KKs(4,4),KKs(4,5),KKs(4,6),KKs(4,7)),T(4));

prob6 = prob5  * SimulateMarkov(F_Q5(KKs(5,1),KKs(5,2),KKs(5,3),KKs(5,4),KKs(5,5),KKs(5,6),KKs(5,7)),1000*max(KKs(5,:)));
output_args = prob6(14);

% X_ans = [F_XX_full(KKs(1),KKs(2),KKs(3),KKs(4),KKs(5),KKs(6),KKs(7)); ones(9,1); zeros(4,1)];

% output_args = prob5 * X_ans(NewIndex);

end

