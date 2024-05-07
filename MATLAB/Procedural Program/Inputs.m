%% Inputs
%--------------------------------------------------------------------------
% I/O

% Folder path
if ismac
    % Code to run on Mac platform
    folderPath      = '../../../../../../../Team Work/Team CrackSTITCH/Datasets/Generic';
elseif isunix
    % Code to run on Linux platform
    folderPath      = '../../../../../../../Team Work/Team CrackSTITCH/Datasets/Generic';
elseif ispc
    % Code to run on Windows platform
    folderPath      = '..\..\..\..\..\..\..\Team Work\Team CrackSTITCH\Datasets\Generic';
else
    disp('Platform not supported')
end

% Folder name that consists of the images set
folderName      = '';

%% Inputs 2
%--------------------------------------------------------------------------
% Parallel workers
input.numCores = 32;            % Number of cores for parallel processing

%% Inputs 3
% Warping
input.warpType = 'spherical';   % 'spherical' | 'cylindrical' | 'planar' (projective)

% Focal length
input.fx = 1600;       % focal length of camera in pixels
input.fy = 1600;       % focal length of camera in pixels

% Distortion coefficients [k1, k2, k3, p1, p2]
input.DC = [0, 0, 0, 0, 0];

% Feature matching
input.detector = 'SIFT';                % 'HARRIS' | 'SIFT' | 'FAST' | 'SURF' | 'BRISK' | 'ORB' | 'KAZE'
input.Matchingmethod = 'Approximate';   %'Exhaustive' (default) | 'Approximate'
input.Matchingthreshold = 1.5;  %       10.0 or 1.0 (default) | percent value in the range (0, 100] | depends on binary and non-binary features
input.Ratiothreshold = 0.6;             % ratio in the range (0,1]
input.NumLayersInOctave = 3;            % Number of layers in each octave -- SIFT only
input.ContrastThreshold = 0.00133;       % Contrast threshold for selecting the strongest features, 
                                        % specified as a non-negative scalar in the range [0,1]. 
                                        % The threshold is used to filter out weak features in 
                                        % low-contrast regions of the image. -- SIFT only
input.EdgeThreshold = 6;                % Edge threshold, specified as a non-negative scalar greater than or equal to 1. 
                                        % The threshold is used to filter out unstable edge-like features  -- SIFT only  

% Image matching (RANSAC)
input.Inliersconfidence = 99.9;         % Inlier confidence [0,100]
input.maxIter = 2000;                   % RANSAC maximum iterations
input.Transformationtype = 'affine'; % 'rigid' | 'similarity' | 'affine' | 'projective'
input.MaxDistance = 1.50;               % Maximum distance (pixels) 1.5

% Image blending and panorama
input.sigmaN = 10;                  % Standard deviations of the normalised intensity error
input.sigmag = 0.1;                 % Standard deviations of the gain
input.resizeImage = 1;              % Resize input images
input.resizeStitchedImage = 0;      % Resize stitched image
input.blending = 'multiband';       % 'multiband' | 'linear' | 'none'
input.bands = 6;                    % bands
input.MBBsigma = 5;                 % Multi-band Gaussian sigma
input.filtSize = [5,5];             % Gaussian kernel Filter size

% Post-processing
input.canvas_color = 'black';       % Panorama canvas color 'black' | 'white'
input.showCropBoundingBox = 1;      % Display cropping bounding box 0 | 1
input.blackRange = 0;               % Minimum dark pixel value to crop panaroma
input.whiteRange = 250;             % Minimum bright pixel value to crop panaroma
input.showPlot  = 0;                % Display keypoints plot (parfor suppresses this flag, so no use)
