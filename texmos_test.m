clear all, close all; clc;

u0     = im2double(imread('tablecloth_and_desktop1.png'))*1;
[m,n] = size(u0);
[x,y] = meshgrid(1:n,1:m);
v0     = zeros(m,n);
a = 0.4;
v0(1:m/2,1:n/2)     = a*cos(2*pi*128/m*x(1:m/2,1:n/2)).*cos(2*pi*128/n*y(1:m/2,1:n/2));%sum(v0(:));
v0(m/2+1:end,1:n/2) = a*cos(2*pi*64/m*x(m/2+1:end,1:n/2));%sum(v0(:));
v0(1:m/2,n/2+1:end) = a*cos(2*pi*64*(x(1:m/2,n/2+1:end)/m+y(1:m/2,n/2+1:end)/n)) ;%sum(v0(:));
v0(m/2+1:end,m/2+1:end) = a*cos((2*pi*128)/m*y(m/2+1:end,1:n/2));

Im     = u0+v0;

sigma  = 3; 
lambda = 1e4;
Param.Reglambda = lambda;
Param.Sigma     = sigma;
Im              = im2double(Im);

tic; [uu,OutPut] = ImSmoothL0TVQP(Im, Param); t=toc;
figure(90); imshow(uu);
figure(91); imshow((Im-uu)+0.5);





%改进后


% clear all, close all; clc;
% 
% % 1. 加载图像并初始化
% filename = 'tablecloth_and_desktop1.png';
% if exist(filename, 'file')
%     u0 = im2double(imread(filename)); % 输入图像
% else
%     error('File %s not found.', filename);
% end
% 
% % 初始化参数
% [m,n] = size(u0);
% [x,y] = meshgrid(1:n,1:m);
% v0     = zeros(m,n);
% a = 0.4;
% 
% % 生成不同纹理模式
% v0(1:m/2,1:n/2)     = a*cos(2*pi*128/m*x(1:m/2,1:n/2)).*cos(2*pi*128/n*y(1:m/2,1:n/2));
% v0(m/2+1:end,1:n/2) = a*cos(2*pi*64/m*x(m/2+1:end,1:n/2));
% v0(1:m/2,n/2+1:end) = a*cos(2*pi*64*(x(1:m/2,n/2+1:end)/m+y(1:m/2,n/2+1:end)/n));
% v0(m/2+1:end,m/2+1:end) = a*cos((2*pi*128)/m*y(m/2+1:end,1:n/2));
% 
% Im     = u0+v0;
% 
% % 2. 设置算法参数
% sigma  = 3; 
% lambda = 1e4;
% beta = 0.1;  % 添加动态参数控制
% tau = 0.01;  % 新增多尺度参数控制
% 
% Param.Reglambda = lambda;
% Param.Sigma     = sigma;
% Param.Beta      = beta;
% Param.Tau       = tau; % 添加多尺度参数
% Im              = im2double(Im);
% 
% % 3. 调用算法
% tic; [uu,OutPut] = ImSmoothL0TVQP(Im, Param); t=toc; % 使用改进版函数
% figure(90); imshow(uu);
% title('Smoothed Image');
% imwrite(uu, 'smoothed_result.png'); % 保存结果
% 
% figure(91); imshow((Im-uu)+0.5);
% title('Residual Image');
% imwrite((Im-uu)+0.5, 'residual_result.png'); % 保存残差图像
% 
% disp(['Processing Time: ', num2str(t), ' seconds']);
