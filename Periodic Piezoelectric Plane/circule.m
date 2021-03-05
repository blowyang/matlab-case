clear;
clc;
% �������м���
core_number=10;            %��Ҫ���õĴ���������
parpool('local',core_number);
% % % % % % ��������������ʾ��
% Starting parallel pool (parpool) using the 'local' profile ...
% connected to 2 workers.
mises_matrix=zeros();
mises_L1=zeros();
%��ά�����������
% n

R1=1;
% z1=-1;
% z2=1;
z1=0;

n=1;
NR=100;
dr_1=R1/(NR+1);

for nr=1:NR
    R1_Show=dr_1*nr;
    COU=4*nr;
    del_ang=2*pi/(COU+1);

    for cou=1:COU
        theta=del_ang*cou;
        zM_circle_1=z1+R1_Show*exp(theta*1i);
        MisesMatrix=0;
        mises_matrix(n,1)=real(zM_circle_1);
        mises_matrix(n,2)=imag(zM_circle_1);
        mises_matrix(n,3)=MisesMatrix;
        n=n+1; 
%         zM_circle_2=z2+R2_Show*exp(theta*1i);
%         MisesMatrix=0;
%         mises_matrix(n,1)=real(zM_circle_2);
%         mises_matrix(n,2)=imag(zM_circle_2);
%         mises_matrix(n,3)=MisesMatrix;
%         n=n+1;
    end    
end

x = mises_matrix(:,1); 
y = mises_matrix(:,2); 
z = mises_matrix(:,3);

%title(['G_{m}:G_{f}=',num2str(Gm),':',num2str(Gf),'; e_{31}E_{0}^{\infty}:\sigma_{11}^{\infty}=',num2str(de/S011)]);%��ӱ���
shp=alphaShape(x,y,'HoleThreshold',0.5);
tri = alphaTriangulation(shp);

patch('faces',tri,'vertices',[x,y],... 
    'facevertexcdata',z,'edgecolor','none','facecolor','interp'); 


axis image
colorbar EastOutside
%% �رղ��м���
delete(gcp('nocreate'));
% % % % % �رպ���������ʾ��
% Parallel pool using the 'local' profile is shutting down.
