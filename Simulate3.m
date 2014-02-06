function dxdt = Simulate3(t,X)

global K SG; % K is vector of K's for each move; SG is spliceGraph structure.

V = K(SG.VmoveIndex) .* X(SG.VnodeIndex);
dxdt = SG.S*V;
end