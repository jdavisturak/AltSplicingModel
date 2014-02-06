function dxdt = test2(t,x)

% Here there are 3 intermediate nodes

S = [-1 -1 -1 0 0 0 ;
       1 0 0 -1 0 0;
       0 1 0 0 -1 0 ;
       0 0 1 0 0 -1;
       0 0 0 1 1 1];

V = x([1 1 1 2 3 4 ]);   % letting all k = 1
dxdt = S*V;
end