%//%*************************************************************************%
%//%*                 Automatic Panorama Image Stitcher                     *%
%//%*           Stitches multiple images using feature points               *%
%//%*                                                                       *%
%//%*                    Name: Dr. Preetham Manjunatha                      *%
%//%*               GitHub: https://github.com/preethamam	                *%
%//%*                   Repo Name: AutoPanoStitch                           *%
%//%*                    Written Date: 04/01/2022                           *%
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

%% Get inputs
%--------------------------------------------------------------------------
% Inputs file
%--------------------------------------------------------------------------
Inputs;

%% Parallel workers start
%--------------------------------------------------------------------------
% Parallel pools
%--------------------------------------------------------------------------
if(isempty(gcp('nocreate')))
    if strcmp(input.poolType,'numcores')
        parpool(input.numCores);
    else
        parpool('Threads')
    end
end

Start = tic;
warning('off','all')


%% Get image filenames and store image names
%--------------------------------------------------------------------------
% Image sets
%--------------------------------------------------------------------------
imgSetVector = imageSet(fullfile(folderPath,folderName),'recursive');
datasetName  = cat(1,{imgSetVector.Description})';
foldersLen = length(imgSetVector);

%% Panorama stitcher
%--------------------------------------------------------------------------
% Stitches panoramas
%--------------------------------------------------------------------------

for myImg = 22 %1:foldersLen

    fprintf('Current folder: %s\n', imgSetVector(myImg).Description);
    
    %% Load images
    tic
    imageFiles = loadImages(imgSetVector, myImg);
    fprintf('Loading images: %f seconds\n', toc);

    %% Get feature matrices and keypoints    
    tic
    [keypoints, allDescriptors, images, imageinfo, imageFocals, numImg] = featureMatching(input, imageFiles);
    fprintf('Feature matching: %f seconds\n', toc);
    
    %% Find matches        
    tic
    [allMatches, numMatches, initialTforms] = imageMatching(input, numImg, images, keypoints, allDescriptors);
    fprintf('Image matching: %f seconds\n', toc);               

    %% Bundle adjustment
    tic
    [finalPanoramaTforms, concomps] = bundleAdjustmentLM(input, images, keypoints, allMatches, numMatches, initialTforms);
    fprintf('Final alignment (Bundle adjustment): %f seconds\n', toc);

    %% Render panoramas
    tic
    [allPanoramas, croppedPanoramas] = displayPanorama(finalPanoramaTforms, input, images, imageFocals, concomps, myImg, datasetName);
    fprintf('Rendering panorama : %f seconds\n', toc);

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
fprintf('Total runtime : %f seconds\n', Runtime);
currtime = datetime('now');
display(currtime)
