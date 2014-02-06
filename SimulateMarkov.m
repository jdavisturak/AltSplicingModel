function Q_T = SimulateMarkov(Q,T)

[V D] = eig(Q);

D2 = diag(exp(diag(D)*T));

Q_T = V*D2*inv(V); % This is equivalent to expm(Q*T)
end
