K=100;
Sm12_AB=zeros(K,1);
Sm12_DB=zeros(K,1);
for k=1:NP
    zk_AB=k*a/(K+1)-a*(1-1i)/2;
    zk_CD=conj(zk_AB);
    zk_DB=k*a*1i/(K+1)+a*(1-1i)/2;
    zk_CA=conj(-zk_DB);
    
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce1,Xdf1,zk_AB);
    Sm12_AB(k,1)=Sm12;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce1,Xdf1,zk_DB);
    
    Sm12_DB(k,1)=Sm12;

end