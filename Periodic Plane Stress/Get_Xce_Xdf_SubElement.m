function [Xce1,Xce2,Xce3,Xce4,Xdf1,Xdf2,Xdf3,Xdf4]=Get_Xce_Xdf_SubElement(K,H_Matrix,B_Matrix,F_Matrix)

P1=zeros(8*K,1);
P2=zeros(8*K,1);
P3=zeros(8*K,1);
P4=zeros(8*K,1);

for i=1:K
    P1(i,1)=1;
    P2(i+K,1)=1;
    P3(i+2*K,1)=1;
    P4(i+3*K,1)=1;
end
%

%Xdf0=-(H_Matrix'*H_Matrix)\(H_Matrix')*C_Matrix*F_Matrix*ELoad;
Xdf1=(H_Matrix'*H_Matrix)\(H_Matrix')*P1;
Xdf2=(H_Matrix'*H_Matrix)\(H_Matrix')*P2;
Xdf3=(H_Matrix'*H_Matrix)\(H_Matrix')*P3;
Xdf4=(H_Matrix'*H_Matrix)\(H_Matrix')*P4;
%Xce0=-F_Matrix*B_Matrix*Xdf0+F_Matrix*ELoad;
Xce1=-F_Matrix*B_Matrix*Xdf1;
Xce2=-F_Matrix*B_Matrix*Xdf2;
Xce3=-F_Matrix*B_Matrix*Xdf3;
Xce4=-F_Matrix*B_Matrix*Xdf4;
% N2=size(Xdf0,1)/2;
% Xd0=zeros(N2,1);
% Xf0=zeros(N2,1);
% Xd1=zeros(N2,1);
% Xf1=zeros(N2,1);
% Xd2=zeros(N2,1);
% Xf2=zeros(N2,1);
% Xd3=zeros(N2,1);
% Xf3=zeros(N2,1);
% Xd4=zeros(N2,1);
% Xf4=zeros(N2,1);
% 
% Xc0=zeros(N2,1);
% Xe0=zeros(N2,1);
% Xc1=zeros(N2,1);
% Xe1=zeros(N2,1);
% Xc2=zeros(N2,1);
% Xe2=zeros(N2,1);
% Xc3=zeros(N2,1);
% Xe3=zeros(N2,1);
% Xc4=zeros(N2,1);
% Xe4=zeros(N2,1);
% for i=1:N2
%     Xd0(i,1)=Xdf0(i,1);
%     Xf0(i,1)=Xdf0(i+N2,1);
%     Xd1(i,1)=Xdf1(i,1);
%     Xf1(i,1)=Xdf1(i+N2,1);
%     Xd2(i,1)=Xdf2(i,1);
%     Xf2(i,1)=Xdf2(i+N2,1);
%     Xd3(i,1)=Xdf3(i,1);
%     Xf3(i,1)=Xdf3(i+N2,1);
%     Xd4(i,1)=Xdf4(i,1);
%     Xf4(i,1)=Xdf4(i+N2,1);
%     
%     Xc0(i,1)=Xce0(i,1);
%     Xe0(i,1)=Xce0(i+N2,1);
%     Xc1(i,1)=Xce1(i,1);
%     Xe1(i,1)=Xce1(i+N2,1);
%     Xc2(i,1)=Xce2(i,1);
%     Xe2(i,1)=Xce2(i+N2,1);
%     Xc3(i,1)=Xce3(i,1);
%     Xe3(i,1)=Xce3(i+N2,1);
%     Xc4(i,1)=Xce4(i,1);
%     Xe4(i,1)=Xce4(i+N2,1);
% end
end