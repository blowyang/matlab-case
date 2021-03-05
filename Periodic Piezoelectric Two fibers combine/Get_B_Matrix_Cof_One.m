function B_Cof=Get_B_Matrix_Cof_One(N,M,R1,R2,x1,x2,a,j)

NP=80;
h=2*pi/NP;
zd=x1;
inte_B_Cof=0;
z0=zd+R1*exp(1i*0);
B_Cof=B_Matrix_Cof_One(N,M,R1,R2,x1,x2,a,z0);
inte_B_Cof=inte_B_Cof+h/6*B_Cof*exp(-j*1i*0);
zNP=zd+R1*exp(1i*2*pi);
B_Cof=B_Matrix_Cof_One(N,M,R1,R2,x1,x2,a,zNP);
inte_B_Cof=inte_B_Cof+h/6*B_Cof*exp(-j*1i*2*pi);

for np=1:(NP-1)
    zi=zd+R1*exp(1i*np*h);
    B_Cof=B_Matrix_Cof_One(N,M,R1,R2,x1,x2,a,zi);
    inte_B_Cof=inte_B_Cof+h/3*B_Cof*exp(-j*1i*np*h);
end

for np=1:NP
    zi=zd+R1*exp(1i*(np-0.5)*h);
    B_Cof=B_Matrix_Cof_One(N,M,R1,R2,x1,x2,a,zi);
    inte_B_Cof=inte_B_Cof+2*h/3*B_Cof*exp(-j*1i*(np-0.5)*h);
end
B_Cof=inte_B_Cof/2/pi;

end