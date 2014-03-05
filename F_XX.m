function output = F_XX(G1,G2,G3,G4,G5,G6,G7,G8,G9,G10,G11,G12,G13,G14,G15,G16,k1,k2,k3,k4,k5,k7)
	output=[(k3+k7)./G1;(G1.*k3+k3.*k4+k4.*k7)./(G1.*G2);(G1.*k7+k2.*k3+k2.*k7)./(G1.*G3);(G1.*G3.*k2.*k3+G1.*G2.*k4.*k7+G2.*k2.*k3.*k4+G3.*k2.*k3.*k4+G2.*k2.*k4.*k7+G3.*k2.*k4.*k7)./(G1.*G2.*G3.*G4);(G1.*k3+k3.*k5+k5.*k7)./(G1.*G5);(G1.*k7+k1.*k3+k1.*k7)./(G1.*G6);(G1.*G6.*k1.*k3+G1.*G5.*k5.*k7+G5.*k1.*k3.*k5+G6.*k1.*k3.*k5+G5.*k1.*k5.*k7+G6.*k1.*k5.*k7)./(G1.*G5.*G6.*G7);(G1.*G6.*k1.*k3+G1.*G2.*k4.*k7+G2.*k1.*k3.*k4+G2.*k1.*k4.*k7+G6.*k1.*k3.*k4+G6.*k1.*k4.*k7)./(G1.*G2.*G6.*G8);(G1.*G3.*G6.*k7+G1.*G3.*k2.*k7+G1.*G6.*k1.*k7+G3.*k1.*k2.*k3+G6.*k1.*k2.*k3+G3.*k1.*k2.*k7+G6.*k1.*k2.*k7)./(G1.*G3.*G6.*G9);(G1.*G3.*k2.*k3+G1.*G5.*k5.*k7+G3.*k2.*k3.*k5+G5.*k2.*k3.*k5+G3.*k2.*k5.*k7+G5.*k2.*k5.*k7)./(G1.*G3.*G5.*G10);(G1.*G2.*G5.*k3+G1.*G2.*k3.*k4+G1.*G5.*k3.*k5+G2.*k3.*k4.*k5+G5.*k3.*k4.*k5+G2.*k4.*k5.*k7+G5.*k4.*k5.*k7)./(G1.*G2.*G5.*G11);(G1.*G2.*G3.*G4.*G10.*k2.*k3.*k4+G1.*G2.*G3.*G4.*G11.*k2.*k3.*k4+G1.*G3.*G4.*G5.*G10.*k2.*k3.*k5+G1.*G2.*G4.*G5.*G11.*k4.*k5.*k7+G1.*G3.*G5.*G10.*G11.*k2.*k3.*k5+G1.*G2.*G5.*G10.*G11.*k4.*k5.*k7+G2.*G3.*G4.*G10.*k2.*k3.*k4.*k5+G2.*G3.*G4.*G11.*k2.*k3.*k4.*k5+G2.*G4.*G5.*G11.*k2.*k3.*k4.*k5+G3.*G4.*G5.*G10.*k2.*k3.*k4.*k5+G2.*G3.*G4.*G10.*k2.*k4.*k5.*k7+G2.*G3.*G4.*G11.*k2.*k4.*k5.*k7+G2.*G4.*G5.*G11.*k2.*k4.*k5.*k7+G3.*G4.*G5.*G10.*k2.*k4.*k5.*k7+G2.*G5.*G10.*G11.*k2.*k3.*k4.*k5+G3.*G5.*G10.*G11.*k2.*k3.*k4.*k5+G2.*G5.*G10.*G11.*k2.*k4.*k5.*k7+G3.*G5.*G10.*G11.*k2.*k4.*k5.*k7+G1.*G2.*G3.*G4.*G5.*G10.*k2.*k3)./(G1.*G2.*G3.*G4.*G5.*G10.*G11.*G12);(G1.*G2.*G6.*G7.*G8.*k1.*k3.*k4+G1.*G2.*G6.*G8.*G11.*k1.*k3.*k4+G1.*G5.*G6.*G7.*G8.*k1.*k3.*k5+G1.*G5.*G6.*G7.*G11.*k1.*k3.*k5+G1.*G2.*G5.*G7.*G11.*k4.*k5.*k7+G1.*G2.*G5.*G8.*G11.*k4.*k5.*k7+G2.*G6.*G7.*G8.*k1.*k3.*k4.*k5+G2.*G5.*G7.*G11.*k1.*k3.*k4.*k5+G2.*G5.*G8.*G11.*k1.*k3.*k4.*k5+G5.*G6.*G7.*G8.*k1.*k3.*k4.*k5+G2.*G6.*G7.*G8.*k1.*k4.*k5.*k7+G2.*G6.*G8.*G11.*k1.*k3.*k4.*k5+G2.*G5.*G7.*G11.*k1.*k4.*k5.*k7+G5.*G6.*G7.*G11.*k1.*k3.*k4.*k5+G2.*G5.*G8.*G11.*k1.*k4.*k5.*k7+G5.*G6.*G7.*G8.*k1.*k4.*k5.*k7+G2.*G6.*G8.*G11.*k1.*k4.*k5.*k7+G5.*G6.*G7.*G11.*k1.*k4.*k5.*k7+G1.*G2.*G5.*G6.*G7.*G8.*k1.*k3)./(G1.*G2.*G5.*G6.*G7.*G8.*G11.*G13);(G1.*G3.*G6.*G7.*G9.*k1.*k2.*k3+G1.*G3.*G6.*G9.*G10.*k1.*k2.*k3+G1.*G3.*G5.*G7.*G10.*k2.*k5.*k7+G1.*G5.*G6.*G7.*G9.*k1.*k5.*k7+G1.*G3.*G5.*G9.*G10.*k2.*k5.*k7+G1.*G5.*G6.*G7.*G10.*k1.*k5.*k7+G3.*G5.*G7.*G10.*k1.*k2.*k3.*k5+G3.*G6.*G7.*G9.*k1.*k2.*k3.*k5+G3.*G5.*G9.*G10.*k1.*k2.*k3.*k5+G5.*G6.*G7.*G9.*k1.*k2.*k3.*k5+G3.*G6.*G9.*G10.*k1.*k2.*k3.*k5+G5.*G6.*G7.*G10.*k1.*k2.*k3.*k5+G3.*G5.*G7.*G10.*k1.*k2.*k5.*k7+G3.*G6.*G7.*G9.*k1.*k2.*k5.*k7+G3.*G5.*G9.*G10.*k1.*k2.*k5.*k7+G5.*G6.*G7.*G9.*k1.*k2.*k5.*k7+G3.*G6.*G9.*G10.*k1.*k2.*k5.*k7+G5.*G6.*G7.*G10.*k1.*k2.*k5.*k7+G1.*G3.*G5.*G6.*G7.*G10.*k5.*k7)./(G1.*G3.*G5.*G6.*G7.*G9.*G10.*G14);(G1.*G3.*G4.*G6.*G9.*k1.*k2.*k3+G1.*G2.*G3.*G4.*G8.*k2.*k4.*k7+G1.*G2.*G3.*G4.*G9.*k2.*k4.*k7+G1.*G2.*G4.*G6.*G8.*k1.*k4.*k7+G1.*G3.*G6.*G8.*G9.*k1.*k2.*k3+G1.*G2.*G6.*G8.*G9.*k1.*k4.*k7+G2.*G3.*G4.*G8.*k1.*k2.*k3.*k4+G2.*G3.*G4.*G9.*k1.*k2.*k3.*k4+G2.*G4.*G6.*G8.*k1.*k2.*k3.*k4+G2.*G3.*G4.*G8.*k1.*k2.*k4.*k7+G2.*G3.*G4.*G9.*k1.*k2.*k4.*k7+G3.*G4.*G6.*G9.*k1.*k2.*k3.*k4+G2.*G4.*G6.*G8.*k1.*k2.*k4.*k7+G2.*G6.*G8.*G9.*k1.*k2.*k3.*k4+G3.*G4.*G6.*G9.*k1.*k2.*k4.*k7+G3.*G6.*G8.*G9.*k1.*k2.*k3.*k4+G2.*G6.*G8.*G9.*k1.*k2.*k4.*k7+G3.*G6.*G8.*G9.*k1.*k2.*k4.*k7+G1.*G2.*G3.*G4.*G6.*G8.*k4.*k7)./(G1.*G2.*G3.*G4.*G6.*G8.*G9.*G15);(G1.*G2.*G3.*G4.*G5.*G6.*G7.*G8.*G9.*G10.*G12.*G14.*G15.*k1.*k2.*k3+G1.*G2.*G3.*G4.*G5.*G6.*G7.*G8.*G9.*G10.*G13.*G14.*G15.*k1.*k2.*k3+G1.*G2.*G3.*G4.*G5.*G6.*G7.*G8.*G10.*G11.*G12.*G13.*G14.*k4.*k5.*k7+G1.*G2.*G3.*G4.*G5.*G6.*G7.*G8.*G10.*G11.*G12.*G13.*G15.*k4.*k5.*k7+G1.*G2.*G3.*G4.*G6.*G7.*G8.*G9.*G10.*G12.*G14.*G15.*k1.*k2.*k3.*k4+G1.*G2.*G3.*G4.*G6.*G7.*G8.*G9.*G11.*G12.*G13.*G15.*k1.*k2.*k3.*k4+G1.*G2.*G3.*G4.*G6.*G7.*G8.*G9.*G10.*G13.*G14.*G15.*k1.*k2.*k3.*k4+G1.*G2.*G3.*G4.*G6.*G7.*G8.*G9.*G11.*G13.*G14.*G15.*k1.*k2.*k3.*k4+G1.*G2.*G3.*G4.*G6.*G8.*G9.*G10.*G11.*G12.*G13.*G15.*k1.*k2.*k3.*k4+G1.*G2.*G3.*G4.*G6.*G8.*G9.*G10.*G11.*G12.*G14.*G15.*k1.*k2.*k3.*k4+G1.*G3.*G4.*G5.*G6.*G7.*G8.*G9.*G10.*G12.*G14.*G15.*k1.*k2.*k3.*k5+G1.*G3.*G4.*G5.*G6.*G7.*G8.*G9.*G10.*G13.*G14.*G15.*k1.*k2.*k3.*k5+G1.*G3.*G4.*G5.*G6.*G7.*G9.*G10.*G11.*G12.*G13.*G14.*k1.*k2.*k3.*k5+G1.*G2.*G3.*G4.*G5.*G7.*G8.*G10.*G11.*G12.*G13.*G14.*k2.*k4.*k5.*k7+G1.*G3.*G4.*G5.*G6.*G7.*G9.*G10.*G11.*G12.*G14.*G15.*k1.*k2.*k3.*k5+G1.*G2.*G3.*G4.*G5.*G7.*G8.*G10.*G11.*G12.*G13.*G15.*k2.*k4.*k5.*k7+G1.*G2.*G3.*G4.*G5.*G7.*G9.*G10.*G11.*G12.*G13.*G14.*k2.*k4.*k5.*k7+G1.*G2.*G4.*G5.*G6.*G7.*G8.*G9.*G11.*G12.*G13.*G15.*k1.*k4.*k5.*k7+G1.*G2.*G4.*G5.*G6.*G7.*G8.*G10.*G11.*G12.*G13.*G14.*k1.*k4.*k5.*k7+G1.*G3.*G5.*G6.*G7.*G8.*G9.*G10.*G11.*G12.*G13.*G14.*k1.*k2.*k3.*k5+G1.*G2.*G3.*G4.*G5.*G7.*G9.*G10.*G11.*G12.*G14.*G15.*k2.*k4.*k5.*k7+G1.*G2.*G3.*G4.*G5.*G8.*G9.*G10.*G11.*G12.*G13.*G15.*k2.*k4.*k5.*k7+G1.*G2.*G4.*G5.*G6.*G7.*G8.*G10.*G11.*G12.*G13.*G15.*k1.*k4.*k5.*k7+G1.*G2.*G3.*G4.*G5.*G8.*G9.*G10.*G11.*G12.*G14.*G15.*k2.*k4.*k5.*k7+G1.*G2.*G4.*G5.*G6.*G7.*G8.*G9.*G11.*G13.*G14.*G15.*k1.*k4.*k5.*k7+G1.*G3.*G5.*G6.*G7.*G8.*G9.*G10.*G11.*G13.*G14.*G15.*k1.*k2.*k3.*k5+G1.*G2.*G5.*G6.*G7.*G8.*G9.*G10.*G11.*G12.*G13.*G14.*k1.*k4.*k5.*k7+G1.*G2.*G5.*G6.*G7.*G8.*G9.*G10.*G11.*G13.*G14.*G15.*k1.*k4.*k5.*k7+G2.*G3.*G4.*G5.*G7.*G8.*G10.*G11.*G12.*G13.*G14.*k1.*k2.*k3.*k4.*k5+G2.*G3.*G4.*G5.*G7.*G8.*G10.*G11.*G12.*G13.*G15.*k1.*k2.*k3.*k4.*k5+G2.*G3.*G4.*G5.*G7.*G9.*G10.*G11.*G12.*G13.*G14.*k1.*k2.*k3.*k4.*k5+G2.*G3.*G4.*G6.*G7.*G8.*G9.*G10.*G12.*G14.*G15.*k1.*k2.*k3.*k4.*k5+G2.*G3.*G4.*G6.*G7.*G8.*G9.*G11.*G12.*G13.*G15.*k1.*k2.*k3.*k4.*k5+G2.*G3.*G4.*G6.*G7.*G8.*G9.*G10.*G13.*G14.*G15.*k1.*k2.*k3.*k4.*k5+G2.*G3.*G4.*G5.*G7.*G9.*G10.*G11.*G12.*G14.*G15.*k1.*k2.*k3.*k4.*k5+G2.*G3.*G4.*G5.*G8.*G9.*G10.*G11.*G12.*G13.*G15.*k1.*k2.*k3.*k4.*k5+G2.*G3.*G4.*G6.*G7.*G8.*G9.*G11.*G13.*G14.*G15.*k1.*k2.*k3.*k4.*k5+G2.*G4.*G5.*G6.*G7.*G8.*G9.*G11.*G12.*G13.*G15.*k1.*k2.*k3.*k4.*k5+G2.*G4.*G5.*G6.*G7.*G8.*G10.*G11.*G12.*G13.*G14.*k1.*k2.*k3.*k4.*k5+G2.*G3.*G4.*G5.*G7.*G8.*G10.*G11.*G12.*G13.*G14.*k1.*k2.*k4.*k5.*k7+G2.*G3.*G4.*G5.*G8.*G9.*G10.*G11.*G12.*G14.*G15.*k1.*k2.*k3.*k4.*k5+G2.*G3.*G4.*G6.*G8.*G9.*G10.*G11.*G12.*G13.*G15.*k1.*k2.*k3.*k4.*k5+G2.*G4.*G5.*G6.*G7.*G8.*G10.*G11.*G12.*G13.*G15.*k1.*k2.*k3.*k4.*k5+G3.*G4.*G5.*G6.*G7.*G8.*G9.*G10.*G12.*G14.*G15.*k1.*k2.*k3.*k4.*k5+G2.*G3.*G4.*G5.*G7.*G8.*G10.*G11.*G12.*G13.*G15.*k1.*k2.*k4.*k5.*k7+G2.*G3.*G4.*G5.*G7.*G9.*G10.*G11.*G12.*G13.*G14.*k1.*k2.*k4.*k5.*k7+G2.*G3.*G4.*G6.*G7.*G8.*G9.*G10.*G12.*G14.*G15.*k1.*k2.*k4.*k5.*k7+G2.*G3.*G4.*G6.*G7.*G8.*G9.*G11.*G12.*G13.*G15.*k1.*k2.*k4.*k5.*k7+G2.*G3.*G4.*G6.*G8.*G9.*G10.*G11.*G12.*G14.*G15.*k1.*k2.*k3.*k4.*k5+G2.*G4.*G5.*G6.*G7.*G8.*G9.*G11.*G13.*G14.*G15.*k1.*k2.*k3.*k4.*k5+G3.*G4.*G5.*G6.*G7.*G8.*G9.*G10.*G13.*G14.*G15.*k1.*k2.*k3.*k4.*k5+G3.*G4.*G5.*G6.*G7.*G9.*G10.*G11.*G12.*G13.*G14.*k1.*k2.*k3.*k4.*k5+G2.*G3.*G4.*G6.*G7.*G8.*G9.*G10.*G13.*G14.*G15.*k1.*k2.*k4.*k5.*k7+G2.*G3.*G4.*G5.*G7.*G9.*G10.*G11.*G12.*G14.*G15.*k1.*k2.*k4.*k5.*k7+G2.*G3.*G4.*G5.*G8.*G9.*G10.*G11.*G12.*G13.*G15.*k1.*k2.*k4.*k5.*k7+G2.*G3.*G4.*G6.*G7.*G8.*G9.*G11.*G13.*G14.*G15.*k1.*k2.*k4.*k5.*k7+G2.*G4.*G5.*G6.*G7.*G8.*G9.*G11.*G12.*G13.*G15.*k1.*k2.*k4.*k5.*k7+G2.*G4.*G5.*G6.*G7.*G8.*G10.*G11.*G12.*G13.*G14.*k1.*k2.*k4.*k5.*k7+G3.*G4.*G5.*G6.*G7.*G9.*G10.*G11.*G12.*G14.*G15.*k1.*k2.*k3.*k4.*k5+G2.*G3.*G4.*G5.*G8.*G9.*G10.*G11.*G12.*G14.*G15.*k1.*k2.*k4.*k5.*k7+G2.*G3.*G4.*G6.*G8.*G9.*G10.*G11.*G12.*G13.*G15.*k1.*k2.*k4.*k5.*k7+G2.*G4.*G5.*G6.*G7.*G8.*G10.*G11.*G12.*G13.*G15.*k1.*k2.*k4.*k5.*k7+G2.*G5.*G6.*G7.*G8.*G9.*G10.*G11.*G12.*G13.*G14.*k1.*k2.*k3.*k4.*k5+G3.*G4.*G5.*G6.*G7.*G8.*G9.*G10.*G12.*G14.*G15.*k1.*k2.*k4.*k5.*k7+G2.*G3.*G4.*G6.*G8.*G9.*G10.*G11.*G12.*G14.*G15.*k1.*k2.*k4.*k5.*k7+G2.*G4.*G5.*G6.*G7.*G8.*G9.*G11.*G13.*G14.*G15.*k1.*k2.*k4.*k5.*k7+G3.*G4.*G5.*G6.*G7.*G8.*G9.*G10.*G13.*G14.*G15.*k1.*k2.*k4.*k5.*k7+G3.*G4.*G5.*G6.*G7.*G9.*G10.*G11.*G12.*G13.*G14.*k1.*k2.*k4.*k5.*k7+G3.*G5.*G6.*G7.*G8.*G9.*G10.*G11.*G12.*G13.*G14.*k1.*k2.*k3.*k4.*k5+G2.*G5.*G6.*G7.*G8.*G9.*G10.*G11.*G13.*G14.*G15.*k1.*k2.*k3.*k4.*k5+G3.*G4.*G5.*G6.*G7.*G9.*G10.*G11.*G12.*G14.*G15.*k1.*k2.*k4.*k5.*k7+G2.*G5.*G6.*G7.*G8.*G9.*G10.*G11.*G12.*G13.*G14.*k1.*k2.*k4.*k5.*k7+G3.*G5.*G6.*G7.*G8.*G9.*G10.*G11.*G13.*G14.*G15.*k1.*k2.*k3.*k4.*k5+G3.*G5.*G6.*G7.*G8.*G9.*G10.*G11.*G12.*G13.*G14.*k1.*k2.*k4.*k5.*k7+G2.*G5.*G6.*G7.*G8.*G9.*G10.*G11.*G13.*G14.*G15.*k1.*k2.*k4.*k5.*k7+G3.*G5.*G6.*G7.*G8.*G9.*G10.*G11.*G13.*G14.*G15.*k1.*k2.*k4.*k5.*k7)./(G1.*G2.*G3.*G4.*G5.*G6.*G7.*G8.*G9.*G10.*G11.*G12.*G13.*G14.*G15.*G16)];
end