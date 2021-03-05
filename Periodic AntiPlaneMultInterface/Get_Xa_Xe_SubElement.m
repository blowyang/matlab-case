function [Xa1,Xa2,Xd1,Xd2,Xe1,Xe2]=Get_Xa_Xe_SubElement(K,H_Matrix,B_Matrix,E_Matrix,N_Matrix)
P1=zeros(4*K,1);
P2=zeros(4*K,1);

for i=1:K
    P1(i,1)=1;
    P2(i+K,1)=1;
end

%Çódelt
Xe1=(H_Matrix'*H_Matrix)\(H_Matrix')*P1;
Xe2=(H_Matrix'*H_Matrix)\(H_Matrix')*P2;

Xd1=-N_Matrix*B_Matrix*Xe1;
Xd2=-N_Matrix*B_Matrix*Xe2;

Xa1=-E_Matrix*B_Matrix*Xe1;
Xa2=-E_Matrix*B_Matrix*Xe2;
end