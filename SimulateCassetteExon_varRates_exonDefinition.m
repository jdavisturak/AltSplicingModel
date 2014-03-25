function [ output_args ] = SimulateCassetteExon_varRates_exonDefinition( KKs, T, E1, E2)
% SimulateCassetteExon_varRates_bareBones( KKs, T, E1, E2 )
%
% KKs: matrix of 5 *7 rate constants (7 constants, for each phase II through VI)
% T: vector of 4 phase durations (for phases II through V)
% Exon definition:
%  E1: fold increase in kinetic rate of recruitment at 5'ss of intron 2, when spliceosome is already recruited to 3'ss of intron 1
%  E2: fold increase in kinetic rate of recruitment at 3'ss of intron 1,  when spliceosome is already recruited to 5'ss of intron 2

prob2 = [1 0] * SimulateMarkov(F_Q2(KKs(1,1)),T(1));

prob3 = [prob2 0 0 0] * SimulateMarkov(F_Q3(KKs(2,1),KKs(2,2),KKs(2,3)),T(2));

prob4 = [prob3 0 0 0 0 0] * SimulateMarkov(F_Q4_exonDefinition(E1,E2, KKs(3,1),KKs(3,2),KKs(3,3),KKs(3,4)),T(3));

prob5 = [prob4 zeros(1,19)] * SimulateMarkov(F_Q5_exonDefinition(E1,E2, KKs(4,1),KKs(4,2),KKs(4,3),KKs(4,4),KKs(4,5),KKs(4,6),KKs(4,7)),T(4));

prob6 = prob5  * SimulateMarkov(F_Q5_exonDefinition(E1,E2, KKs(5,1),KKs(5,2),KKs(5,3),KKs(5,4),KKs(5,5),KKs(5,6),KKs(5,7)),1000*max(KKs(5,:)));

output_args = prob6(14);

end
