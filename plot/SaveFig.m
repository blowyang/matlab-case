%����ͼ�񲢱�����ָ��λ��
count=1;
x = -10:0.1:10; % x��ֵ
y = x.^2; % y��ֵ
fig = figure; % �½�һ��figure������ͼ�������浽fig
plot(x,y,'.') % ��"."����ʽ��x��y�������������ɵ�ͼ����
legend({'y=x^2'},'Location','northwest') % ��ͼ������Ͻ�����ͼ��
frame = getframe(fig); % ��ȡframe
img = frame2im(frame); % ��frame�任��imwrite��������ʶ��ĸ�ʽ
%imwrite(img,'a.png'); % ���浽����Ŀ¼�£�����Ϊ"a.png"
imwrite(img,['count=',num2str(count),'a.png']); % ���浽����Ŀ¼�£�����λ�ÿɴ�����