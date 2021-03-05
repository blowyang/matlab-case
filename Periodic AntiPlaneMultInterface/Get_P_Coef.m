function P_Coef=Get_P_Coef(lambda1,x1,x2)
%lambda1=R1/a
y=control(x1,x2);
P_Coef=y*(143/84*lambda1)^x2;
%P_Coef=y;
end