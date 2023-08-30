classdef Tomography
    %TOMOGRAPHYIMAGE Compute tomography image from object 3d's slices
    %{ 
    %====================================================
    AUTHOR      Mingling Luo
    CONTACT     ml-lrem@outlook.com
    DATE        2023
    %====================================================
    USAGE       [obj] = TomographyImag e(slices)
    
    properties  slices       - feild*layers array         - Slices of 3d object 
                params       - struct                     - optimize, projector, material params
                errorRates   - 1*MaxIter    array         - iter error rate
                images       - img*angles   array         - Tomography images in all angle for projection
        
    Methods     Configure                       - Configure all params
                Compute                         - Compute tomography images
    %====================================================
    %}
    properties
        slices;                                                             % 源切片数据
        optimizeParams;
        projectorParams;
        materialParams;
        errorRates;
        images;                                                             % 最终输出的断层扫描数据，也即各角度的投影图像
    end
    
    methods
        function obj = Tomography(slices)
            obj.slices = slices;
        end

        function obj = Configure(obj, optimizeParams, projectorParams, materialParams)
            %CONFIGURE Configure projector\material\optimize params
            obj.optimizeParams = optimizeParams;
            obj.projectorParams = projectorParams;
            obj.materialParams = materialParams;
            % default values: optimize
            if ~isfield(optimizeParams, 'filter')
                obj.optimizeParams.filter = true;                           % 启用滤波器
            end
            if ~isfield(optimizeParams, 'maxIter')
                obj.optimizeParams.maxIter = 50;                            % 优化迭代次数
            end
            if ~isfield(optimizeParams, 'sigmoid')
                obj.optimizeParams.omega = 100;                             % omega参数是sigmoid函数的连接权值，该值越大函数的阶梯性越明显，意味着阈值分割的锐度越大
            end                                                             
            if ~isfield(optimizeParams, 'learningRate')
                obj.optimizeParams.learningRate = 0.005;                   % 迭代梯度下降法时的步长参数，一般取0.001~0.1
            end
            if ~isfield(optimizeParams, 'beta')
                obj.optimizeParams.beta = 0;                                % ??
            end
            if ~isfield(optimizeParams, 'theta')
                obj.optimizeParams.theta = 0;                               % ??
            end
            if ~isfield(optimizeParams, 'rho')
                obj.optimizeParams.rho = 0;                                 % ??
            end
            % default values: project
            if ~isfield(projectorParams, 'deltaAngle')   
                projectorParams.deltaAngle = 1;                                             % 旋转精度，决定了成型的水平/横向分辨率
            end
            obj.projectorParams.angles = linspace(0, 179, ceil(180/projectorParams.deltaAngle));
            if ~isfield(projectorParams, 'bit8')   
                obj.projectorParams.bit8 = false;                           % 投影仪8bit
            end
            if ~isfield(projectorParams, 'equalize8bit')   
                obj.projectorParams.equalize8bit = false;                   % ??
            end
            if ~isfield(projectorParams, 'equalize8bit')   
                obj.projectorParams.equalize8bit = false;                   % ??
            end
            if ~isfield(projectorParams, 'mask')
                obj.projectorParams.mask = false;                           % 投影遮罩
            end
            % default values: matrial
            if ~isfield(materialParams, 'threshold')
                obj.materialParams.threshold = NaN;                         % 决定树脂是否被固化的分割阈值
            end
        end
        
        function obj = Compute(obj,show_switch)
            %COMPUTE Compute tomography image using projector\material\optimize params
            angles = obj.projectorParams.angles;
            threshold = obj.materialParams.threshold;
            omega = obj.optimizeParams.omega;
            beta = obj.optimizeParams.beta;
            learningRate = obj.optimizeParams.learningRate;
            theta = obj.optimizeParams.theta;

            % STEP1: Radon变换并执行频域滤波   
            rSlices = CALFun.Radon(obj.slices, obj.projectorParams.angles);
            if obj.optimizeParams.filter
                rSlices = CALFun.Filter(rSlices, 'ram-lak');
                rSlices = max(rSlices, 0);                                  % 将滤波后出现的负数置零,与零比较取最大值
            end
            
            % STEP2: 反投影迭代
            rSlicesDeltaPre = single(zeros(size(rSlices)));
            obj.errorRates = zeros(obj.optimizeParams.maxIter,1);
            for iter = 1:obj.optimizeParams.maxIter
                 % 调整位深度和对比度
                rSlices = single(uint8(rSlices/max(rSlices(:))*255))/255;
                rSlices = imadjustn(rSlices);

                % 反投影，用于计算误差和梯度
                slicesBack = CALFun.IRadon(rSlices, angles);                   
                slicesBack = slicesBack/max(slicesBack(:));                
                slicesBackThreshold = CALFun.Sigmoid(slicesBack, threshold, omega);   % 使用sigmoid函数实现阈值分类，threshold为分界点
                % 与源切片的差值: Delta. 该值就是梯度下降法中当次迭代的梯度
                slicesDelta = slicesBackThreshold - obj.slices;
                obj.errorRates(iter) = CALFun.ErrorRate(obj.slices, slicesBack);     % 【为什么要用slicesBack, 而不用slicesBackThreshold? 前者才能真实反映树脂各点吸收的光的能量，也才能反映材料的固化情况】

                % 松弛/动量梯度下降: gradientMomentum-动量梯度(指数加权平均)-斜率; 
                rSlicesDelta = CALFun.Radon(slicesDelta, angles);
                gradientMomentum = ((1-beta)*rSlicesDelta + beta*rSlicesDeltaPre)/(1-beta^iter);
                rSlices = rSlices - gradientMomentum*learningRate;                    % w' = w - k*learningRate
                rSlices = rSlices.*(double(rSlices >= 0) + theta*double(rSlices < 0)); % 施加正向约束             
                rSlicesDeltaPre = rSlicesDelta;

                % 迭代停止
                if show_switch
                    Show.DynamicPlot(obj.errorRates,'iter','errorRsate');
                end
                if obj.errorRates(iter) <= 1e-5
                    break;
                end
            end
            slcesOptimized = slicesBackThreshold;
            rSlicesOptimized = rSlices;
            
            if show_switch
                Show.Slices(obj.slices,'vol');
                Show.Slices(slcesOptimized,'vol');
            end

            % STEP3: 转换为断层扫描图像=各角度投影仪图像
            obj.images = permute(rSlicesOptimized,[3,1,2]);
        end
        
        function ShowImages(obj, show_switch)
            if show_switch
                Show.Img3D(obj.images, 'Tomography images / Projection images');
            end
        end
    end
end

