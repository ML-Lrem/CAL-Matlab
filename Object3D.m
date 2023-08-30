classdef Object3D
    %OBJECT3D 3D object Class, support {'.STL',...}
    %{ 
    %====================================================
    AUTHOR      Mingling Luo
    CONTACT     ml-lrem@outlook.com
    DATE        2023
    %====================================================
    USAGE       [obj] = Object3D(file_name)
    
    properties  file_name   - string          - 3D object's read path
                resolution  - numerical       - 3D object's slice resolution [units:um]
                fv          - n*[Dim] array   - 3D object's source data, include face & vertices
                slices      - n*[Dim] array   - 3D object's slices
                
                xyz_min/max - 2*[Dim] array   - XYZ direction's edge: Min(xyz),Max(xyz)
                xyz_size    - 1*[Dim] array   - X,Y,Z direction's range size
                xyz_count   - 1*[Dim] array   - The count of voxels in the XYZ direction at
                                                the specified resolution
    
    Methods     Slice()             - Generate 3D object's slices
                ShowSlices()        - Show 3D object's slices
    %====================================================
    %}
    properties
        file_name;
        slice_layers;
        fv;
        slices;
        dim = 3;
        xyz_min;
        xyz_max;
        xyz_size;
    end
    
    methods
        function obj = Object3D(file_name)
            %Object3D Create 3D object and import data from file name
            fprintf("Read STL ... cost: ");
            time_start = tic;
            %================
            obj.file_name = file_name;
            obj.fv = CALFun.STLRead(file_name);       
            obj.xyz_min = zeros(1,obj.dim);
            obj.xyz_max = zeros(1,obj.dim);
            obj.xyz_size = zeros(1,obj.dim);
            for i = 1:obj.dim
                obj.xyz_min(i) = min(obj.fv.vertices(:,i));
                obj.xyz_max(i) = max(obj.fv.vertices(:,i));
            end
            for i = 1:obj.dim
                obj.xyz_size(i) = obj.xyz_max(i)-obj.xyz_min(i);      % Max - Min
            end
            %================
            time_stlread = toc(time_start);
            fprintf("%0.2f s\n",time_stlread);
        end

        function obj = Slice(obj,resolution,mode)
            %Slice Slice the 3D object to slices
            fprintf("Slice STL ... cost: ");
            time_start = tic;
            %================
            if mode == "layer"
                obj.slice_layers = resolution;
            elseif mode == "absolute"
                obj.slice_layers = obj.xyz_size(3)/resolution;
            end
            xyz_count = zeros(1,obj.dim);       % xyz方向的体素数量
            for i = 1:obj.dim
                xyz_count(i) = ceil(obj.slice_layers*(obj.xyz_size(i)/obj.xyz_size(3)));
                % 调整count的奇偶性, 保证之后使用padarray时XY方向上的体素数量是一致的
                xyz_count(i) = xyz_count(i) - mod(xyz_count(i),2) + 1;          % 确保xyz_count为奇数    
            end
            gridX = linspace(obj.xyz_min(1), obj.xyz_max(1),xyz_count(1));      % X方向上的所有体素顶点位置
            gridY = linspace(obj.xyz_min(2), obj.xyz_max(2),xyz_count(2));      % Y方向上的所有体素顶点位置
            gridZ = linspace(obj.xyz_min(3), obj.xyz_max(3),xyz_count(3));      % Z方向上的所有体素顶点位置
            obj.slices = VOXELISE(gridX,gridY,gridZ,obj.fv);                    % 体素化
            r_count = ceil(sqrt(xyz_count(1)^2 + xyz_count(2)^2));              % 将切片拓展为一个最小圆面可包裹的范围
            r_count = r_count - mod(r_count,2) + 1;
            pad_range = [0.5*(r_count-xyz_count(1)), 0.5*(r_count-xyz_count(2))];
            obj.slices = padarray(obj.slices, pad_range, 0, 'both');            % 根据最小圆，填充XY方向上的体素
            % Empty the 3D object source data that no longer need
            obj.fv = [];
            %================
            time_voxel = toc(time_start);
            fprintf("%0.2f s\n",time_voxel);
        end

        function ShowSlices(obj,show_switch)
            %SHOWSLICES Show 3D object's slices in 3D viewer
            if show_switch
                Show.Slices(obj.slices,"vol");
            end
        end
    end
end
