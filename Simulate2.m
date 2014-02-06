function dxdt = Simulate2(t,X)

global K S VnodeIndex VmoveIndex;

V = K(VmoveIndex) .* X(VnodeIndex);
dxdt = S*V;
end