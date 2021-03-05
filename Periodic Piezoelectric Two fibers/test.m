S022=1;
S011=0;
S012=0;
[delta1,delta2,delta3,delta4,delta_cof]=Solve_Delta(R1,R2,x1,x2,a,Xce1,Xce2,Xce3,Xce4,...
    Xdf1,Xdf2,Xdf3,Xdf4,S022,S011,S012);
% delta1=1;
% delta2=1;
% delta3=1;
% delta4=1;
% 
Xce=delta1*Xce1+delta2*Xce2+delta3*Xce3+delta4*Xce4;
Xdf=delta1*Xdf1+delta2*Xdf2+delta3*Xdf3+delta4*Xdf4;


NP=100;
del=2*pi/(NP+1);
stress=zeros(NP,3);
for np=1:NP
    theta=np*del;
    z=x1+R1*exp(1i*theta);
    
    [Sm22,Sm11,Sm12]=Get_Matrix_Stress(R1,R2,x1,x2,a,Xce,Xdf,z);
    
    Smrr=(Sm22+Sm11)/2+(Sm11-Sm22)/2*cos(2*theta)+Sm12*sin(2*theta);
    Smss=(Sm22+Sm11)/2-(Sm11-Sm22)/2*cos(2*theta)-Sm12*sin(2*theta);
    stress(np,1)=theta;
    stress(np,2)=Smss;
    stress(np,3)=Smrr;
%  
%     S13_DB_Sub_CA(np,1)=S13_DB(np,1)-S13_CA(np,1);
%     S23_AB_Sub_CD(np,1)=S23_AB(np,1)-S23_CD(np,1);
end