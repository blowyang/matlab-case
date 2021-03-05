function G_Rk=Get_G1(G0,G2,R1,R2,Rk)

DR=R1-R2;
dr=Rk-R2;
G_Rk=G2*power(10,dr/DR*log10(G0/G2));
%G_Rk=(G0-G2)/DR*dr+G2;
end