clear;
clc;
% 启动并行计算
core_number=10;            %想要调用的处理器个数
parpool('local',core_number);
% % % % % % 启动后有如下提示：
% Starting parallel pool (parpool) using the 'local' profile ...
% connected to 2 workers.
%在确定的模量比例下，和确定的界面相厚度的前提下，改变纤维的体积分数VF2,研究
%等效切边模量的变化规律
filename= 'mov2.gif';
N=20;
M=15;
K=60;
Gm=10;
Gf=1;
kf=1.8;
km=1.8;
R=1;

e31=-6.5;
%E3infinity=0.2;

S022=0;
S011=1;
S012=0;
MisesMatrixInf=((S022^2+S011^2+(S022-S011)^2+6*S012^2)/2)^0.5;
vf=0.3;
a=R*(pi/vf)^0.5;
[E_Matrix,F_Matrix]=Get_E_F_Matrix(N,Gm,Gf,kf,km);
B_Matrix=Get_B_Matrix(N,M,km,R,a);

[C_Matrix,D_Matrix]=Get_CD_Matrix(K,N,M,a,R,km,Gm);
H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,F_Matrix,B_Matrix);



mises_matrix=zeros();
mises_L1=zeros();
%纤维体积分数增量
% n
for i=1:101
    de=-1+(i-1)*2/100;
    E3infinity=de/e31;
    ELoad=Get_ELoad(N,e31,E3infinity,R);
    [Xce0,Xce1,Xce2,Xce3,Xce4,Xdf0,Xdf1,Xdf2,Xdf3,Xdf4]=Get_Xce_Xdf_SubElement(K,H_Matrix,B_Matrix,C_Matrix,F_Matrix,ELoad);

    [delta1,delta2,delta3,delta4,delta_cof]=Solve_Delta(a,R,Xce0,Xce1,Xce2,Xce3,Xce4,...
        Xdf0,Xdf1,Xdf2,Xdf3,Xdf4,S022,S011,S012);
    Xce=Xce0+delta1*Xce1+delta2*Xce2+delta3*Xce3+delta4*Xce4;
    Xdf=Xdf0+delta1*Xdf1+delta2*Xdf2+delta3*Xdf3+delta4*Xdf4;
    
    n=1;
    NP=40;
    del_1=0.75*a/(NP+1);
    COU=100;
    del_2=2*pi/(COU+1);

    for np=1:NP+2
        Ri=R+del_1*(np-1);

        for cou=1:COU+1
            zx=Ri*cos(del_2*cou);
            zy=Ri*sin(del_2*cou);
            if(abs(zx)<=a/2) && (abs(zy)<=a/2)
                zM_np=zx+zy*1i;
                MisesMatrix=Get_MisesMatrix(a,R,Xce,Xdf,zM_np)/MisesMatrixInf;
                mises_matrix(n,1)=zx;
                mises_matrix(n,2)=zy;
                mises_matrix(n,3)=MisesMatrix;
                n=n+1;
            end

        end
    end

    x = mises_matrix(:,1); 
    y = mises_matrix(:,2); 
    z = mises_matrix(:,3);
    fig = figure;
    % xlable('');
    % ylable('');
    %'G_{m}:G_{f}=',num2str(Gm),':',num2str(Gf)\sigma{0}^{2}+\beta_{0}^{2}=1'
    %title(['Gm=',num2str(Gm)]);
    title(['G_{m}:G_{f}=',num2str(Gm),':',num2str(Gf),'; e_{31}E_{0}^{\infty}:\sigma_{11}^{\infty}=',num2str(de/S011)]);%添加标题
    shp=alphaShape(x,y,'HoleThreshold',1.5);
    tri = alphaTriangulation(shp);

    patch('faces',tri,'vertices',[x,y],... 
        'facevertexcdata',z,'edgecolor','none','facecolor','interp'); 

    %hal=patch('faces',tri,'vertices',[x,y],'facevertexcdata',z,'edgecolor','none','facecolor','interp'); 
    axis image
    colorbar EastOutside 
    frame=getframe(fig);
    im=frame2im(frame);%制作gif文件，图像必须是index索引图像
    [I,map]=rgb2ind(im,256);
    k=i-0;
    if k==1
        imwrite(I,map,filename,'gif','Loopcount',inf,...
            'DelayTime',0.1);
    else
        imwrite(I,map,filename,'gif','WriteMode','append',...
            'DelayTime',0.1);
    end
end

% frame = getframe(fig); % 获取frame
% img = frame2im(frame); % 将frame变换成imwrite函数可以识别的格式
% imwrite(img,'a.png'); % 保存到工作目录下，名字为"a.png"
% for cou=1:COU
%     zM_np=R1*cos(del_2*cou)+R1*sin(del_2*cou)*1i;
%     MisesMatrix=Get_MisesMatrix_Divided_G0(Xd,Xe,a,R1,zM_np); 
%     mises_L1(cou,1)=del_2*cou;
%     mises_L1(cou,2)=MisesMatrix;
% end

%% 关闭并行计算
delete(gcp('nocreate'));
% % % % % 关闭后有如下提示：
% Parallel pool using the 'local' profile is shutting down.
