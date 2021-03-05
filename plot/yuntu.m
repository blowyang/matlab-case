clear; 
clc;
M=xlsread('mises_matrix.xlsx');
x = M(:,1); y = M(:,2); z = M(:,3);

shp=alphaShape(x,y,'HoleThreshold',1.5);
tri = alphaTriangulation(shp);

patch('faces',tri,'vertices',[x,y],... % 'facecolor' == 'interp' ?
    'facevertexcdata',z,'edgecolor','none','facecolor','interp'); % 

% patch('faces',tri,'vertices',[x,y],... % 'facecolor' == 'flat' 
%     'facevertexcdata',mean(z(tri),2),'edgecolor','none','facecolor','flat'); % 

axis image
colorbar EastOutside