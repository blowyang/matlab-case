%生成图像并保存至指定位置
count=1;
x = -10:0.1:10; % x赋值
y = x.^2; % y赋值
fig = figure; % 新建一个figure，并将图像句柄保存到fig
plot(x,y,'.') % 用"."的形式将x，y表现在上面生成的图像中
legend({'y=x^2'},'Location','northwest') % 在图像的左上角生成图例
frame = getframe(fig); % 获取frame
img = frame2im(frame); % 将frame变换成imwrite函数可以识别的格式
%imwrite(img,'a.png'); % 保存到工作目录下，名字为"a.png"
imwrite(img,['count=',num2str(count),'a.png']); % 保存到工作目录下，保存位置可带变量