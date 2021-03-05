function [L1,L2,L3,L4,L5]=Get_L_List(N,M,R1,R2,x1,x2,a,z)
L1=zeros(N,1);
L2=zeros(N,1);
L3=zeros(N,1);
L4=zeros(N,1);

L5=zeros(M,1);

for n=1:N
    L1(n,1)=((z-x1)/R1)^n;
    L2(n,1)=((z-x2)/R2)^n;
    L3(n,1)=(R1/(z-x1))^n;
    L4(n,1)=(R2/(z-x2))^n;
end
for m=1:M
%     L2(m,1)=Get_P(m,a,z);
%     L4(m,1)=Get_P_Diff(m,a,z);
%     L6(m,1)=Get_P_2Diff(m,a,z);    
    
     L5(m,1)=(z/a)^m;
%    L5(m,1)=Get_P(m,a,z);

end
end