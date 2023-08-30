%{
% Compute the 3D target geometry on 3D object
%============================================
% AUTHOR  Mingling Luo
% CONTACT ml-lrem@outlook.com
%============================================
%}

%% Define
clc;
clear;
close all;

%% Global Configure
showOn = true;
showOff = false;

%% Import & Slice STL
voxelResolution = 0.2;     % 体素精度, 按物理尺寸计(mm)，自适应切片层数【PS: 实际打印精度还受到投影仪焦平面像素宽度限制】

targetSTL = Object3D("stl\head.stl");      % stl单位默认为mm，如非如此需要提前进行单位转换
targetSTL = targetSTL.Slice(voxelResolution,'absolute');
targetSTL.ShowSlices(showOff);

%% Compute Tomography Image
projectorParams.deltaAngle = 1;         % 水平方向投影精度，单位：角度[0,180];
materialParams.threshold = 0.8;
optimizeParams.filter = true;
optimizeParams.learningRate = 0.001;
optimizeParams.maxIter = 150;

tomography = Tomography(targetSTL.slices);
tomography = tomography.Configure(optimizeParams,projectorParams,materialParams);
tomography = tomography.Compute(showOn);
projectorImages = tomography.images;
tomography.ShowImages(showOn);


