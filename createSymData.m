function mySyms=createSymData(N)
% Create a string that we can evaluate to create tons of symbolic variables
A=arrayfun(@(x)sprintf('sym(''m%i'') ',x),1:N,'UniformOutput',false);
mySyms = eval(['[', A{:} ,']']);
end

