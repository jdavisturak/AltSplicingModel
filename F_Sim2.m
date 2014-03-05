function output = F_Sim2(k1,t2)
	output=reshape([exp(-k1.*t2),0.0,-exp(-k1.*t2)+1.0,1.0],[2,2]);
end