clc;
clear ;
close all;

num_image = 8;
dst_dir = 'GIF\';

filename= 'bird.gif'; %���gif�ļ�������
for i=1:num_image
    idx=sprintf('%03d',i);
    str=[dst_dir idx '.jpg'];
    Img=imread(str);
    figure(i)
    imshow(Img);
    frame=getframe(i);
    im=frame2im(frame);%����gif�ļ���ͼ�������index����ͼ��
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