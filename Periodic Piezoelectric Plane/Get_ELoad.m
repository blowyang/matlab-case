function ELoad=Get_ELoad(N,e31,E3infinity,R)
ELoad=zeros(8*N,1);
ELoad(1,1)=e31*E3infinity*R;
end