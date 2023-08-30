classdef CALFun
    %CALFUNCTION A collection of other functions that need to be used in CAL
    %{ 
    %====================================================
    AUTHOR      Mingling Luo
    CONTACT     ml-lrem@outlook.com
    DATE        2023
    %====================================================
    USAGE       result = CALFun.MyFun()

    Methods     Radon            - Configure all params
                IRadon           - Compute tomography images
                Filter           - Ram-Lak filter
                Sigmod           - Sigmod Function
                ErrorRate        - Compute iter error rate
    %====================================================
    %}
    methods (Static=true)
        function result = Radon(slices,angles)
            %CALRADON compute slices's radon;   foward projection
            layers = size(slices,3);                                                % 切片层数
            slicesSize = size(slices,1);                                            % 一层切片的大小 = XY方向上的体素数量
            result = zeros(slicesSize,length(angles),layers);
            for z = 1:layers
                [resultLayer,~] = radon(slices(:,:,z),angles);
                % 沿中心切除多余数据，保证Radon变换后的投影数据和源数据尺寸相同，便于后续处理。
                % 之所以能切除多于部分，是因为在获取切片时，对切片进行了拓展
                % 在切除时，需要保证切除数量一定为整数，避免出现不合理的浮点
                if mod(size(resultLayer,1),2) == mod(slicesSize,2)
                    d = size(resultLayer,1)/2 - slicesSize/2;
                    rxCount = d + 1;
                    ryCount = size(resultLayer,1) - d;
                else
                    d = size(resultLayer,1)/2 - (slicesSize+1)/2;
                    rxCount = d + 1;
                    ryCount = size(resultLayer,1) - d - 1;
                end
                resultLayer = resultLayer((rxCount:ryCount),:);
                result(:,:,z) = resultLayer;
            end
            result = single(result);
            % 并行算法执行Radon，能大幅加快运算速度
        end

        function result = IRadon(rSlices,angles)
            %CALRADON compute slices's iradon;  back projection
            layers = size(rSlices,3);                                       % 切片层数
            slicesSize = size(rSlices,1);                                   % 一层切片的大小 = XY方向上的体素数量
            
            result = zeros(slicesSize,slicesSize,layers);
            for z = 1:layers
                [resultLayer,~] = iradon(rSlices(:,:,z),angles,'none',slicesSize);  % 不滤波
                result(:,:,z) = resultLayer;
            end
            % 将切片时填充的区域置零：这些区域是进行切片时为了在变换时不损失精度而故意填充的，在进行逆radon变换时，需要将这些区域置0
            [Y,X] = meshgrid(linspace(-1,1,size(result,1)),linspace(-1,1,size(result,2)));
            R = sqrt(X.^2 + Y.^2);                                          % R大于1的部分为填充域
            R = repmat(R,[1,1,size(result,3)]);                             % 重复R坐标系下的副本
            result= result.*(R<=1);
            result = single(result);
            % 并行算法执行Radon，能大幅加快运算速度
        end
        
        % 理解/修改/优化
        function [p_out,H] = Filter(p_in, filter)
            p = p_in;
            d = 1;      % 截止频率
            % Design the filter
            len = size(p,1);
            H = designFilter(filter, len, d);
            if strcmpi(filter, 'none')
                return;
            end          
            p_out = zeros(length(H),size(p,2),size(p,3));     
            for nz=1:size(p,3)
                p_z = p(:,:,nz);    % 2维
                p_z(length(H),1)=0;     % 将p_z的行数拓展到H，保证傅里叶变换后的频域范围 
                p_z = fft(p_z);         % 180个角度，每个角度对应一个一维傅里叶变换，更加合理的写法应该是：pz_fft = fft(pz) or fpz = fft(pz)   
                p_z = bsxfun(@times, p_z, H);       % 数组乘法，傅里叶变换后乘以滤波器，频域处理
                p_z = ifft(p_z,'symmetric');
                p_out(:,:,nz) = p_z;
            end 
            p_out(len+1:end,:,:) = [];  % 后面的都是无效的频域范围，只取主要部分 
            p_out = single(p_out);
        end
        
        function result = Sigmoid(x,mu,omega)
            %CALSIGMOID sigmoid function, Limit output to [0,1]
            result = 1./(1+exp(-(x-mu)*omega));
        end
        
        % 理解/修改/优化
        function result = ErrorRate(target,recon)
            [X,Y,~] = meshgrid(linspace(-1,1,size(target,1)),linspace(-1,1,size(target,1)),linspace(-1,1,size(target,3)));
            R = sqrt(X.^2 + Y.^2);
            circle_mask = logical(R.*(R<=1));
            gel_inds = circle_mask & target==1;                             % 目标对象的圆域内的固化区域和圆域内非固化区域
            void_inds = circle_mask & ~target;          

            num_gel_void = sum(gel_inds(:)) + sum(void_inds(:));
            min_gel_dose = min(recon(gel_inds),[],'all');                   % 重建的目标中，固化区域的计算最小值
            void_doses = recon(void_inds);                                  % 重建的目标中，非固化区域的所有值
            n_pix_overlap = sum(void_doses>=min_gel_dose);                  % 非固化区域大于等于固化区域最小值的所有数量 = 溢出数量，意味着在实际打印时固化和非固化会产生异常
            result = n_pix_overlap/num_gel_void;                            % 异常数量占总数的比例，当该比例下降到0时，所有重建质量达标
        end
        
        function varargout = STLRead(file)
            % STLREAD imports geometry from an STL file into MATLAB.
                if ~exist(file,'file')
                    error(['File ''%s'' not found. If the file is not on MATLAB''s path' ...
                           ', be sure to specify the full path to the file.'], file);
                end
                
                fid = fopen(file,'r');    
                if ~isempty(ferror(fid))
                    error(lasterror); %#ok
                end
                
                M = fread(fid,inf,'uint8=>uint8');
                fclose(fid);
                
                [f,v,n] = stlbinary(M);
                %if( isbinary(M) ) % This may not be a reliable test
                %    [f,v,n] = stlbinary(M);
                %else
                %    [f,v,n] = stlascii(M);
                %end
                
                varargout = cell(1,nargout);
                switch nargout        
                    case 2
                        varargout{1} = f;
                        varargout{2} = v;
                    case 3
                        varargout{1} = f;
                        varargout{2} = v;
                        varargout{3} = n;
                    otherwise
                        varargout{1} = struct('faces',f,'vertices',v);
                end
            
            end
    end
end

function filt = designFilter(filter, len, d)
% Returns the Fourier Transform of the filter which will be
% used to filter the projections
%
% INPUT ARGS:   filter - either the string specifying the filter
%               len    - the length of the projections
%               d      - the fraction of frequencies below the nyquist
%                        which we want to pass
%
% OUTPUT ARGS:  filt   - the filter to use on the projections


order = max(64,2^nextpow2(2*len));

if strcmpi(filter, 'none')
    filt = ones(1, order);
    return;
end

% First create a bandlimited ramp filter (Eqn. 61 Chapter 3, Kak and
% Slaney) - go up to the next highest power of 2.

n = 0:(order/2); % 'order' is always even. 
filtImpResp = zeros(1,(order/2)+1); % 'filtImpResp' is the bandlimited ramp's impulse response (values for even n are 0)
filtImpResp(1) = 1/4; % Set the DC term 
filtImpResp(2:2:end) = -1./((pi*n(2:2:end)).^2); % Set the values for odd n
filtImpResp = [filtImpResp filtImpResp(end-1:-1:2)]; 
filt = 2*real(fft(filtImpResp)); 
filt = filt(1:(order/2)+1);

w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist

switch filter
    case 'ram-lak'
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
    case 'cosine'
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
    case 'hamming'
        filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
    case 'hann'
        filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
    case 'superlinear'
        filt(2:end) = filt(2:end).^3;

    otherwise
        error(message('images:iradon:invalidFilter'))
end

filt(w>pi*d) = 0;                      % Crop the frequency response
filt = [filt' ; filt(end-1:-1:2)'];    % Symmetry of the filter
%----------------------------------------------------------------------
end

function [F,V,N] = stlbinary(M)

    F = [];
    V = [];
    N = [];
    
    if length(M) < 84
        error('MATLAB:stlread:incorrectFormat', ...
              'Incomplete header information in binary STL file.');
    end
    
    % Bytes 81-84 are an unsigned 32-bit integer specifying the number of faces
    % that follow.
    numFaces = typecast(M(81:84),'uint32');
    %numFaces = double(numFaces);
    if numFaces == 0
        warning('MATLAB:stlread:nodata','No data in STL file.');
        return
    end
    
    T = M(85:end);
    F = NaN(numFaces,3);
    V = NaN(3*numFaces,3);
    N = NaN(numFaces,3);
    
    numRead = 0;
    while numRead < numFaces
        % Each facet is 50 bytes
        %  - Three single precision values specifying the face normal vector
        %  - Three single precision values specifying the first vertex (XYZ)
        %  - Three single precision values specifying the second vertex (XYZ)
        %  - Three single precision values specifying the third vertex (XYZ)
        %  - Two unused bytes
        i1    = 50 * numRead + 1;
        i2    = i1 + 50 - 1;
        facet = T(i1:i2)';
        
        n  = typecast(facet(1:12),'single');
        v1 = typecast(facet(13:24),'single');
        v2 = typecast(facet(25:36),'single');
        v3 = typecast(facet(37:48),'single');
        
        n = double(n);
        v = double([v1; v2; v3]);
        
        % Figure out where to fit these new vertices, and the face, in the
        % larger F and V collections.        
        fInd  = numRead + 1;        
        vInd1 = 3 * (fInd - 1) + 1;
        vInd2 = vInd1 + 3 - 1;
        
        V(vInd1:vInd2,:) = v;
        F(fInd,:)        = vInd1:vInd2;
        N(fInd,:)        = n;
        
        numRead = numRead + 1;
    end
    
end

function [F,V,N] = stlascii(M)
    warning('MATLAB:stlread:ascii','ASCII STL files currently not supported.');
    F = [];
    V = [];
    N = [];
end

% TODO: Change the testing criteria! Some binary STL files still begin with
% 'solid'.
function tf = isbinary(A)
% ISBINARY uses the first line of an STL file to identify its format.
    if isempty(A) || length(A) < 5
        error('MATLAB:stlread:incorrectFormat', ...
              'File does not appear to be an ASCII or binary STL file.');
    end    
    if strcmpi('solid',char(A(1:5)'))
        tf = false; % ASCII
    else
        tf = true;  % Binary
    end
end
