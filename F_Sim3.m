function output = F_Sim3(k1,k2,k3,t3)
	output=reshape([exp(-t3.*(k1+k2)),0.0,0.0,0.0,0.0,-exp(-t3.*(k1+k2))+exp(-k2.*t3),exp(-k2.*t3),0.0,0.0,0.0,-exp(-t3.*(k1+k2))+exp(-k1.*t3),0.0,exp(-k1.*t3),0.0,0.0,(exp(-t3.*(k1+k2)).*(k1+k2))./(k1+k2-k3)-(k1.*exp(-k1.*t3))./(k1-k3)-(k2.*exp(-k2.*t3))./(k2-k3)+(k1.*k2.*exp(-k3.*t3).*(k1+k2-k3.*2.0))./((k1-k3).*(k2-k3).*(k1+k2-k3)),-(k2.*exp(-k2.*t3))./(k2-k3)+(k2.*exp(-k3.*t3))./(k2-k3),-(k1.*exp(-k1.*t3))./(k1-k3)+(k1.*exp(-k3.*t3))./(k1-k3),exp(-k3.*t3),0.0,(k3.*exp(-k1.*t3))./(k1-k3)+(k3.*exp(-k2.*t3))./(k2-k3)-(k3.*exp(-t3.*(k1+k2)))./(k1+k2-k3)-(k1.*k2.*exp(-k3.*t3).*(k1+k2-k3.*2.0))./((k1-k3).*(k2-k3).*(k1+k2-k3))+1.0,-(k2.*exp(-k3.*t3))./(k2-k3)+(k3.*exp(-k2.*t3))./(k2-k3)+1.0,-(k1.*exp(-k3.*t3))./(k1-k3)+(k3.*exp(-k1.*t3))./(k1-k3)+1.0,-exp(-k3.*t3)+1.0,1.0],[5,5]);
end