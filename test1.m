function dxdt = test1(t,x)

S = [-1 -1 0 0;
       1 0 -1 0;
       0 1 0 -1;
       0 0 1 1];

V = x([1 1 2 3]);   % letting all k = 1
dxdt = S*V;
end