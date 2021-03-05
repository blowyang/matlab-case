NP=100;
del=2*pi/NP;
MisesMatrix=zeros();
for np=1:NP
    theta=del*np;
    Mises=((1+cos(2*theta))^2+(sin(2*theta))^2)^0.5;
    MisesMatrix(np,1)=Mises;
end