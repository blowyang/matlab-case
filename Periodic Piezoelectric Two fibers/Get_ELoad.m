function ELoad=Get_ELoad(N,e31_1,E3infinity_1,R1,e31_2,E3infinity_2,R2)
ELoad=zeros(16*N,1);
ELoad(1,1)=e31_1*E3infinity_1*R1;
ELoad(8*N+1,1)=e31_2*E3infinity_2*R2;
end