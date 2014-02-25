function mySyms=createSymData(N, start)
% Create a string that we can evaluate to create tons of symbolic variables

A=arrayfun(@(x)sprintf('sym(''%s%i'') ',start,x),1:N,'UniformOutput',false);
mySyms = eval(['[', A{:} ,']']);
end

