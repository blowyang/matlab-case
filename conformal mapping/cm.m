clear;
clc;
K=100;
del=2*pi/(K+1);
x=zeros(K,1);
y=zeros(K,1);
for i=1:K
    theta=del*i;
    sigma=exp(1i*theta);
    z=sigma+1/15*sigma^(-5);
    x(i,1)=real(z);
    y(i,1)=imag(z);
end
plot(x,y);