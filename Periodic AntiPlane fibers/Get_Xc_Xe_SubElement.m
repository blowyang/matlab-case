function [Xcd1,Xcd2,Xe1,Xe2]=Get_Xc_Xe_SubElement(K,H_Matrix,B_Matrix,F_Matrix)
P1=zeros(4*K,1);
P2=zeros(4*K,1);

for i=1:K
    P1(i,1)=1;
    P2(i+K,1)=1;
end

%Çódelt
Xe1=(H_Matrix'*H_Matrix)\(H_Matrix')*P1;
Xe2=(H_Matrix'*H_Matrix)\(H_Matrix')*P2;


Xcd1=-F_Matrix*B_Matrix*Xe1;
Xcd2=-F_Matrix*B_Matrix*Xe2;
end