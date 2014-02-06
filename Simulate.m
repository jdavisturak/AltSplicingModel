function dxdt = Simulate(t,X)

global S VnodeIndex;

%FluxByNode = (X.*VnodeMultiplier);
V = X(VnodeIndex);
dxdt = S*V;
end