% Copyright May. 25, 2019, Dr.WEN You-Wei
% email: wenyouwei@gmail.com

classdef ClassBlurMatrix
    properties
        psf
        eigBM
        boundarycond  = 'cir'
        tgv        = 0
    end
    
    methods
        function ob = ClassBlurMatrix(psf, imsize, flag)
            if nargin < 3, flag = 'cir';   end
            ob.boundarycond = flag;
            if length(imsize)>2 && imsize(3)>1 && size(psf,3) == 1
                psf  = repmat(psf,1,1,imsize(3));
            end
            ob.psf   = psf;
            ob.eigBM = zeros(imsize);
            
            s        = fix(([size(psf,1) size(psf,2)]+1)/2);
            switch ob.boundarycond
                case {'cir','circular'}
                    ob.eigBM(1:size(psf,1), 1:size(psf,2),:) = psf;
                    ob.eigBM = circshift(ob.eigBM, 1-s);
                    ob.eigBM = fft2(ob.eigBM);
                case {'refl','symmetric'}
                    psfa = psf(end:-1:1,end:-1:1,:);
                    if norm(psfa(:)-psf(:),'fro')>1e-7
                        error('The blur kernel should be symetric!');
                    end
                    h1 = psf(s(1):end, s(2):end,:);
                    h2 = h1(2:end,:,:); h2(s(1),:,:) = 0;
                    h3 = h1 + h2;
                    h4 = h3(:,2:end,:); h4(:,s(2),:) = 0;
                    h5 = h4 + h3;
                    
                    ob.eigBM(1:s(1),1:s(2),:) = h5;
                    e1 = zeros(imsize); e1(1,1,:) = 1;
                    for k = 1:imsize(3)
                        ob.eigBM(:,:,k) = dct2(ob.eigBM(:,:,k))./dct2(e1(:,:,k));
                    end
                    %ob.ForTF         = @(x)dct2(x);
                    %ob.BackTF        = @(x)idct2(x);
                otherwise
                    error('Wrong. Flag should be cir or refl or imfilter');
            end
        end
        
        function ob = ctranspose(A) %% written by wenyouwei@gmail.com
            ob       = A;
            ob.eigBM = conj(A.eigBM);
        end
        function y = mldivide(A,x)
            if ~isa(A,'ClassBlurMatrix')
                error('A must be the class of ClassBlurMatrix');
            end
            B = A; B.eigBM = 1./A.eigBM;
            y = B * x; 
        end
        
        function y = mrdivide(x,A)
            B = A; B.eigBM = 1./A.eigBM;
            y = B * x; 
        end
        
        function ob = abs(A) %% written by wenyouwei@gmail.com
            ob       = A;
            ob.eigBM = abs(A.eigBM);
        end
        
        function ob = inv(A) %% written by wenyouwei@gmail.com
            ob       = A;
            ob.eigBM = 1./(A.eigBM);
        end
        
        function ob = plus(a,b)
            if isa(a,'ClassBlurMatrix')
                ob = a;
                if isa(b,'ClassBlurMatrix')
                    ob.eigBM = a.eigBM + b.eigBM;
                else
                    ob.eigBM = a.eigBM + b;
                end
            else
                ob = b;
                ob.eigBM = a + b.eigBM;
            end
        end
        function ob = minus(a,b)
            if isa(a,'ClassBlurMatrix')
                ob = a;
                if isa(b,'ClassBlurMatrix')
                    ob.eigBM = a.eigBM - b.eigBM;
                else
                    ob.eigBM = a.eigBM - b;
                end
            else
                ob       = b;
                ob.eigBM = a - b.eigBM;
            end
        end
        
        
        function y = mtimes(a,x)%% written by wenyouwei@gmail.com
            if ~isa(a, 'ClassBlurMatrix')
                if numel(a) ~= 1
                    error('Wrong, 1th input: a scalar or a class');
                else
                    y = x;
                    y.eigBM = a * x.eigBM;
                end
                return;
            end
            if isa(x,'ClassBlurMatrix')
                y = a;    y.eigBM = y.eigBM .* x.eigBM;
                return;
            end
            
            if numel(x) == 1
                y = a;    y.eigBM = y.eigBM  * x;
                return;
            end
            switch a.boundarycond
                case {'cir','circular'}
                    y = ifft2(a.eigBM  .* fft2(x));
                    y = real(y); 
                case {'refl','symmetric'}
                    y = x; 
                    for j = 1:size(x,3)
                        y(:,:,j) = idct2(a.eigBM(:,:,j) .* dct2(x(:,:,j)));
                    end
            end
        end
        
    end
    
end







% Copyright May. 25, 2019, Dr.WEN You-Wei
% email: wenyouwei@gmail.com

%{
classdef ClassBlurMatrix
    % ClassBlurMatrix 类用于表示具有不同边界条件的模糊矩阵。它允许处理与模糊矩阵相关的操作。
    
    properties
        psf                 % PSF（点扩散函数，表示模糊核）
        eigBM               % 储存模糊矩阵的特征（通常为频域表示）
        boundarycond = 'cir' % 边界条件，默认为'cir'（圆形边界）
        tgv = 0             % 变量，用于控制TGV（全变差正则化，total generalized variation）
    end
    
    methods
        % 构造函数，初始化 ClassBlurMatrix 对象
        function ob = ClassBlurMatrix(psf, imsize, flag)
            if nargin < 3, flag = 'cir'; end   % 如果没有传入 flag 参数，默认使用 'cir'
            ob.boundarycond = flag;           % 设置边界条件
            
            % 如果输入图像是三维的（例如彩色图像），而PSF是单通道的
            if length(imsize) > 2 && imsize(3) > 1 && size(psf, 3) == 1
                psf = repmat(psf, 1, 1, imsize(3));  % 将PSF扩展到每个通道
            end
            
            ob.psf = psf;             % 存储模糊核
            ob.eigBM = zeros(imsize); % 初始化模糊矩阵（特征矩阵），大小与图像相同
            
            s = fix(([size(psf, 1), size(psf, 2)] + 1) / 2); % 计算PSF的中心点
            
            switch ob.boundarycond
                case {'cir', 'circular'} % 圆形边界条件
                    ob.eigBM(1:size(psf, 1), 1:size(psf, 2), :) = psf; % 将PSF赋给模糊矩阵
                    ob.eigBM = circshift(ob.eigBM, 1 - s);  % 将PSF在矩阵中循环位移
                    ob.eigBM = fft2(ob.eigBM);  % 对模糊矩阵进行二维傅里叶变换
                case {'refl', 'symmetric'} % 对称边界条件
                    psfa = psf(end:-1:1, end:-1:1, :);  % 反转PSF以满足对称性
                    if norm(psfa(:) - psf(:), 'fro') > 1e-7
                        error('The blur kernel should be symmetric!');  % 如果模糊核不对称，抛出错误
                    end
                    
                    % 计算图像分块，以适应对称边界条件
                    h1 = psf(s(1):end, s(2):end, :);  
                    h2 = h1(2:end, :, :); h2(s(1), :, :) = 0;
                    h3 = h1 + h2;
                    h4 = h3(:, 2:end, :); h4(:, s(2), :) = 0;
                    h5 = h4 + h3;
                    
                    ob.eigBM(1:s(1), 1:s(2), :) = h5;  % 将结果赋值给模糊矩阵
                    e1 = zeros(imsize); e1(1, 1, :) = 1;
                    for k = 1:imsize(3)
                        ob.eigBM(:, :, k) = dct2(ob.eigBM(:, :, k)) ./ dct2(e1(:, :, k)); % 对模糊矩阵进行离散余弦变换（DCT）
                    end
                otherwise
                    error('Wrong. Flag should be cir or refl or imfilter');  % 错误处理：不支持的边界条件
            end
        end
        
        % 共轭转置操作
        function ob = ctranspose(A)
            ob = A;
            ob.eigBM = conj(A.eigBM);  % 对模糊矩阵进行共轭转置
        end
        
        % 求解线性方程组 A * x = b，A 为 ClassBlurMatrix 类型
        function y = mldivide(A, x)
            if ~isa(A, 'ClassBlurMatrix')
                error('A must be the class of ClassBlurMatrix');  % 确保 A 是 ClassBlurMatrix 类型
            end
            B = A; B.eigBM = 1 ./ A.eigBM;  % 将 A 的模糊矩阵元素反转
            y = B * x;  % 计算解 B * x
        end
        
        % 右除操作
        function y = mrdivide(x, A)
            B = A; B.eigBM = 1 ./ A.eigBM;  % 将 A 的模糊矩阵元素反转
            y = B * x;  % 计算解 x / A
        end
        
        % 求模（绝对值）
        function ob = abs(A)
            ob = A;
            ob.eigBM = abs(A.eigBM);  % 对模糊矩阵进行元素级别的绝对值操作
        end
        
        % 求逆操作
        function ob = inv(A)
            ob = A;
            ob.eigBM = 1 ./ A.eigBM;  % 对模糊矩阵元素进行求逆操作
        end
        
        % 加法操作
        function ob = plus(a, b)
            if isa(a, 'ClassBlurMatrix')
                ob = a;
                if isa(b, 'ClassBlurMatrix')
                    ob.eigBM = a.eigBM + b.eigBM;  % 两个 ClassBlurMatrix 对象相加
                else
                    ob.eigBM = a.eigBM + b;  % 一个 ClassBlurMatrix 对象和标量相加
                end
            else
                ob = b;
                ob.eigBM = a + b.eigBM;  % 标量和 ClassBlurMatrix 对象相加
            end
        end
        
        % 减法操作
        function ob = minus(a, b)
            if isa(a, 'ClassBlurMatrix')
                ob = a;
                if isa(b, 'ClassBlurMatrix')
                    ob.eigBM = a.eigBM - b.eigBM;  % 两个 ClassBlurMatrix 对象相减
                else
                    ob.eigBM = a.eigBM - b;  % ClassBlurMatrix 对象减去标量
                end
            else
                ob = b;
                ob.eigBM = a - b.eigBM;  % 标量减去 ClassBlurMatrix 对象
            end
        end
        
        % 乘法操作
        function y = mtimes(a, x)
            if ~isa(a, 'ClassBlurMatrix')
                if numel(a) ~= 1
                    error('Wrong, 1th input: a scalar or a class');
                else
                    y = x;
                    y.eigBM = a * x.eigBM;  % 标量与 ClassBlurMatrix 对象相乘
                end
                return;
            end
            if isa(x, 'ClassBlurMatrix')
                y = a;    y.eigBM = y.eigBM .* x.eigBM;  % 两个 ClassBlurMatrix 对象逐元素相乘
                return;
            end
            
            if numel(x) == 1
                y = a;    y.eigBM = y.eigBM * x;  % ClassBlurMatrix 对象与标量相乘
                return;
            end
            switch a.boundarycond
                case {'cir', 'circular'}  % 圆形边界条件
                    y = ifft2(a.eigBM .* fft2(x));  % 使用傅里叶变换进行卷积操作
                    y = real(y);  % 获取实部
                case {'refl', 'symmetric'}  % 对称边界条件
                    y = x; 
                    for j = 1:size(x, 3)
                        y(:, :, j) = idct2(a.eigBM(:, :, j) .* dct2(x(:, :, j)));  % 对每个通道进行 DCT 操作
                    end
            end
        end
    end
end
%}