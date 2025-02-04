% function [u,OutPut] = ImSmoothL0TVQP(f, Param)
% 
% 
% 
% MaxIter = 5;
% SolRE   = 1e-4;
% kappa   = 2;            lambda  = 200;
% betamax = 1e5;          sigma   = 2;
% % beta    = 1;     
% %betamax = 1e5; 测试 1e5,1e8,1e10,1e15
% if nargin == 2
%     if isfield(Param,'PrimalStep'), kappa     = Param.kappa;      end
%     if isfield(Param,'MaxIter'),    MaxIter   = Param.MaxIter;    end
%     if isfield(Param,'SolRE'),      SolRE     = Param.SolRE;      end
%     if isfield(Param,'Reglambda'),  lambda    = Param.Reglambda;  end
%     if isfield(Param,'Sigma'),      sigma     = Param.Sigma;      end
%     if isfield(Param,'Beta'),       beta      = Param.Beta;       end
%     if isfield(Param,'BetaMax'),    betamax   = Param.BetaMax;    end
% end
% kouter  = 0;       k   = 0;
% [M,N,D]=size(f);
% sr  = 2 * ceil(3*sigma) + 1;
% psf = fspecial('gaussian', [sr,sr], sigma);
% H   = ClassBlurMatrix(psf,size(f));
% g   = f; 
% Hg  = H * g;        u   = Hg;
% Htg = H' * Hg;      HtH = H' * H;
% 
% R   = ClassImGrad('cir',size(g),1);%R = ClassFoImDiff('cir',size(g));
% RtR = R' * R;      Lap = RtR;
% 
% if ~exist('beta','var'), Ru  = abs(R*u);     beta = 1/mean(Ru(:));  end   
% 
% lamiter = lambda/8; 
% ratio   = floor(log(betamax/beta)/log(kappa))-1;
% ratio   = 1.5^(5/ratio);
% 
% while beta < betamax
%     kouter = kouter + 1;
%     A      = lamiter * HtH + beta*Lap;
%     for i = 1:MaxIter
%         k   = k + 1;
%         z = ImGradHardTh(R * u,1/beta);
%         %zz  = sqrt(sum(sum(z.^2,4),3));
%         %figure(1); imshow(1-zz/max(zz(:)));
%         %figure(2); imshow(u); pause(0.05)
%         
%         rhs = lamiter * Htg + beta * R' * z;
%         unew = A\rhs;
%        % re = norm(unew(:)-u(:))/norm(u(:));
%        
%         
%         u = unew; %u(u>1)=1; u(u<0) = 0;
%         re=norm(unew(:)-u(:))/norm(u(:));
%      %   e = H * u - f; norm(e(:))^2;
%         
%   %      fprintf('%3d-th iteration, beta: %1.2e, relative error: %1.2e\n',k,beta,re);
%         if re<SolRE, break; end
%     end
%     lamiter = min(lamiter*ratio, lambda); 
%     beta    = beta*kappa;
%     
%     
% end
%  
% 
% [lamiter kouter]
% %figure; imshow(u);
% OutPut.OuterIter = kouter;
% OutPut.TotalIter = k;
% 
% function [z,zz] = ImGradHardTh(z,lambda)
% zz  = sqrt(sum(sum(z.^2,4),3));
% % zz  = imerode(zz,ones(3,3));
% % zz  = imdilate(zz,ones(3,3));
% mu  = 2*sqrt(lambda);  %mu =max(mu,5e-3);
% idx = zz<mu;
% zz(idx) = 0;zz(idx) = 0;
% idx = repmat(idx,[1 1 size(z,3) size(z,4)]);
% z(idx)  = 0;
% end
% 
% function z = ImGradHardTh2(z,lambda)
% % zz  = sqrt(sum(z.^2,3));
% mu  = 2*sqrt(lambda);  %mu =max(mu,5e-3);
% % idx = zz<mu;
% % idx = repmat(idx,[1 1 size(z,3) 1]);
% % z(idx) = 0;
% 
% if D==1
%         idx = sqrt(sum(z.^2))<mu;
%     else
%         idx = sqrt(sum(z.^2,3))<mu;
% end
% idx = repmat(idx,[1 1 D 1]);
% z(idx) = 0;
% end
% 
% end





function [u,OutPut] = ImSmoothL0TVQP(f, Param)



MaxIter = 5;
SolRE   = 1e-4;
kappa   = 2;            lambda  = 200;
betamax = 1e5;          sigma   = 2;
beta    = 1;     

if nargin == 2
    if isfield(Param,'PrimalStep'), kappa     = Param.kappa;      end
    if isfield(Param,'MaxIter'),    MaxIter   = Param.MaxIter;    end
    if isfield(Param,'SolRE'),      SolRE     = Param.SolRE;      end
    if isfield(Param,'Reglambda'),  lambda    = Param.Reglambda;  end
    if isfield(Param,'Sigma'),      sigma     = Param.Sigma;      end
    if isfield(Param,'Beta'),       beta      = Param.Beta;       end
    if isfield(Param,'BetaMax'),    betamax   = Param.BetaMax;    end
end
kouter  = 0;       k   = 0;

sr  = 2 * ceil(3*sigma) + 1;
psf = fspecial('gaussian', [sr,sr], sigma);
H   = ClassBlurMatrix(psf,size(f));
g   = f; 
Hg  = H * g;        u   = Hg;
Htg = H' * Hg;      HtH = H' * H;

R   = ClassImGrad('cir',size(g),1);%R = ClassFoImDiff('cir',size(g));
RtR = R' * R;      Lap = RtR;

if ~exist('beta','var'), Ru  = abs(R*u);     beta = 1/mean(Ru(:));  end   

lamiter = lambda/8; 
ratio   = floor(log(betamax/beta)/log(kappa))-1;
ratio   = 1.5^(5/ratio);

while beta < betamax
    kouter = kouter + 1;
    A      = lamiter * HtH + beta*Lap;
    for i = 1:MaxIter
        k   = k + 1;
        [z,zz] = ImGradHardTh(R * u,1/beta);
        zz  = sqrt(sum(sum(z.^2,4),3));
        figure(1); imshow(1-zz/max(zz(:)));
        figure(2); imshow(u); pause(0.05)
        
        rhs = lamiter * Htg + beta * R' * z;
        unew = A\rhs;
        re = norm(unew(:)-u(:))/norm(u(:));
        
        u = unew; %u(u>1)=1; u(u<0) = 0;
        e = H * u - f; norm(e(:))^2;
        
       fprintf('%3d-th iteration, beta: %1.2e, relative error: %1.2e\n',k,norm(e(:))^2,lamiter);
        if re<SolRE, break; end
    end
    lamiter = min(lamiter*ratio, lambda); 
    beta    = beta*kappa;
    
    
end
 

[lamiter kouter]
figure; imshow(u);
OutPut.OuterIter = kouter;
OutPut.TotalIter = k;

function [z,zz] = ImGradHardTh(z,lambda)
zz  = sqrt(sum(sum(z.^2,4),3));
zz  = imerode(zz,ones(3,3));
zz  = imdilate(zz,ones(3,3));
mu  = 2*sqrt(lambda);  %mu =max(mu,5e-3);
idx = zz<mu;
zz(idx) = 0;zz(idx) = 0;
idx = repmat(idx,[1 1 size(z,3) size(z,4)]);
z(idx)  = 0;


function z = ImGradHardTh2(z,lambda)
zz  = sqrt(sum(z.^2,3));
mu  = 2*sqrt(lambda);  %mu =max(mu,5e-3);
idx = zz<mu;
idx = repmat(idx,[1 1 size(z,3) 1]);
z(idx) = 0;







%改进后

% function [u,OutPut] = ImSmoothL0TVQP_Extended(f, Param)
% 
% 
% MaxIter = 5;
% SolRE   = 1e-4;         % 收敛阈值
% kappa   = 2;            % 增长因子
% lambda  = 200;          % 正则化参数
% betamax = 1e5;          % 最大惩罚系数
% sigma   = 2;            % 高斯平滑参数
% tau     = 0.01;         % 多尺度平滑参数
% beta    = 1;            % 初始惩罚参数
% 
% 
% if nargin == 2
%     if isfield(Param,'PrimalStep'), kappa     = Param.kappa;      end
%     if isfield(Param,'MaxIter'),    MaxIter   = Param.MaxIter;    end
%     if isfield(Param,'SolRE'),      SolRE     = Param.SolRE;      end
%     if isfield(Param,'Reglambda'),  lambda    = Param.Reglambda;  end
%     if isfield(Param,'Sigma'),      sigma     = Param.Sigma;      end
%     if isfield(Param,'Beta'),       beta      = Param.Beta;       end
%     if isfield(Param,'BetaMax'),    betamax   = Param.BetaMax;    end
%     if isfield(Param,'Tau'),        tau       = Param.Tau;        end
% end
% 
% kouter  = 0;       k   = 0;
% 
% 
% sr  = 2 * ceil(3*sigma) + 1;
% psf = fspecial('gaussian', [sr,sr], sigma);
% H   = ClassBlurMatrix(psf,size(f));
% g   = f; 
% Hg  = H * g;        
% u   = Hg;           % 初始图像
% Htg = H' * Hg;      
% HtH = H' * H;
% 
% R   = ClassImGrad('cir',size(g),1);
% RtR = R' * R;      
% Lap = RtR;         % 一阶梯度正则化
% 
% if ~exist('beta','var'), Ru  = abs(R*u);     beta = 1/mean(Ru(:));  end   
% 
% lamiter = lambda/8; 
% ratio   = floor(log(betamax/beta)/log(kappa))-1;
% ratio   = 1.5^(5/ratio);
% 
% 
% while beta < betamax
%     kouter = kouter + 1;
%     A      = lamiter * HtH + beta * Lap + tau * del2(u); % 添加二阶正则化
%     for i = 1:MaxIter
%         k   = k + 1;
% 
%         % 更新梯度和硬阈值
%         [z,zz] = ImGradHardTh(R * u,1/beta);
% 
%         % 更新右侧项
%         rhs = lamiter * Htg + beta * R' * z + tau * del2(u); % 二阶梯度控制
%         unew = A \ rhs; % 解线性方程组
% 
%         % 判断收敛条件
%         re = norm(unew(:)-u(:))/norm(u(:));
%         u = unew; 
% 
%         if re < SolRE
%             fprintf('Iteration %d converged with relative error %.2e\n', k, re);
%             break; 
%         end
%     end
% 
%     % 更新参数
%     lamiter = min(lamiter * ratio, lambda); 
%     beta    = beta * kappa;
%     tau     = max(tau / 1.2, 1e-3); % 动态调整tau
% 
%     % 输出进度
%     fprintf('Outer Iter %d, beta = %.2e, tau = %.2e, re = %.2e\n', kouter, beta, tau, re);
% end
% 
% 
% OutPut.OuterIter = kouter;
% OutPut.TotalIter = k;
% 
% 
% function [z,zz] = ImGradHardTh(z,lambda)
% zz  = sqrt(sum(z.^2,3)); 
% mu  = 2*sqrt(lambda);  
% idx = zz < mu;
% idx = repmat(idx,[1 1 size(z,3)]);
% z(idx) = 0;
