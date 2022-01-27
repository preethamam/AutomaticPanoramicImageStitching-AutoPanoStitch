%//%*************************************************************************%
%//%*                 Automatic Panorama Image Stitcher                     *%
%//%*           Stitches multiple images using feature points               *%
%//%*                                                                       *%
%//%*                    Name: Preetham Manjunatha                          *%
%//%*               GitHub: https://github.com/preethamam	                *%
%//%*                   Repo Name: AutoPanoStitch                           *%
%//%*                    Written Date: 07/28/2021                           *%
%%***************************************************************************%
%* Citation 1: Automatic Panoramic Image Stitching using Invariant Features.*% 
%* M. Brown and D. Lowe. International Journal of Computer Vision. 74(1),   *%
%* pages 59-73, 2007                                                        *%
%*                                                                          *%
%* Citation 2: Recognising Panoramas. M. Brown and D. G. Lowe.              *%
%* International Conference on Computer Vision (ICCV2003). pages 1218-1225, *%
%* Nice, France, 2003.                                                      *%
%****************************************************************************%

%% Start
%--------------------------------------------------------------------------
clear; close all; clc;
clcwaitbarz = findall(0,'type','figure','tag','TMWWaitbar');
delete(clcwaitbarz);
Start = tic;
warning('off','all')

%% Inputs
%--------------------------------------------------------------------------
if ismac
    % Code to run on Mac platform
elseif isunix
    % Code to run on Linux platform
    folderPath      = '../../../Data/Generic';
elseif ispc
    % Code to run on Windows platform
    folderPath      = '..\..\..\Data\Generic';
else
    disp('Platform not supported')
end

folderName      = '';
fileExt         = {'.jpg','.jpeg','.png', '.tiff', '.JPG','.JPEG','.PNG', '.TIFF'};
input.showPlot  = 0;

%% Inputs 2
% Feature matching
input.detector = 'SIFT'; % 'HARRIS' | 'SIFT' | 'FAST' | 'SURF' | 'BRISK' | 'ORB' | 'KAZE'
input.Matchingmethod = 'Approximate'; %'Exhaustive' (default) | 'Approximate'
input.Matchingthreshold = 1.5; %10.0 or 1.0 (default) | percent value in the range (0, 100] | depends on binary and non-binary features
input.Ratiothreshold = 0.6; % ratio in the range (0,1]

% Image matching (RANSAC)
input.Inliersconfidence = 99.9;
input.maxIter = 2000;
input.Transformationtype = 'projective'; %'rigid' | 'similarity' | 'affine' | 'projective'
input.MaxDistance = 1.50; %1.5;

% Image blending and panorama
input.sigmaN = 10;
input.sigmag = 0.1;
input.resizeImage = 1;
input.resizeStitchedImage = 0;
input.blending = 'multiband';       % 'multiband' | 'linear' | 'none'
input.bands = 2;

% Post-processing
input.canvas_color = 'black';
input.showCropBoundingBox = 0;
input.blackRange = 0;
input.whiteRange = 250;

%% Get image filenames and store image names
imgSetVector = imageSet(fullfile(folderPath,folderName),'recursive');
datasetName  = cat(1,{imgSetVector.Description});

%% Extract SIFT features
folderLen = length(imgSetVector);

for myImg = 8 %1:folderLen

    fprintf('Current folder: %s\n', imgSetVector(myImg).Description);
    
    %% Get feature matrices and keypoints    
    tic
    [keypoints, allDescriptors, images, imageinfo, imageFocals, numImg] = featureMatching(input, imgSetVector, myImg);
    toc
    
    %% Find matches        
    tic
    [allMatches, numMatches, initialTforms] = imageMatching(input, numImg, images, keypoints, allDescriptors);
    toc        
    
    %% Bundle adjustment
    tic
    [finalPanoramaTforms, concomps] = bundleAdjustmentLM(input, images, keypoints, allMatches, numMatches, initialTforms);
    toc

    %% Render panoramas
    tic
    displayPanorama(finalPanoramaTforms, input, images, imageFocals, concomps, myImg, datasetName);
    toc
end

%% End parameters
%--------------------------------------------------------------------------
clcwaitbarz = findall(0,'type','figure','tag','TMWWaitbar');
delete(clcwaitbarz);
statusFclose = fclose('all');
if(statusFclose == 0)
    disp('All files are closed.')
end
Runtime = toc(Start);
disp(Runtime);
currtime = datetime('now');
display(currtime)
