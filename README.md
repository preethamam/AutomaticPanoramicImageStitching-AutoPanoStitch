# AutoPanoStitch
Automatic Panorama Stitching software in MATLAB. Spherical, cyclindrical and planar projections stitching is supported in this version and can recognize multiple panoramas.

# Stitched images 1:
| Type | Images |
| --- | --- |
| Stitched image | ![pano_full](https://user-images.githubusercontent.com/28588878/151394796-907b2a27-2054-412a-aa6c-aa5120294df5.jpg) |
| Crop box | ![pano_bbox](https://user-images.githubusercontent.com/28588878/151394950-cb1c0009-2ed4-4b2b-94dc-66dc18695445.jpg) |
| Cropped image | ![pano_crop](https://user-images.githubusercontent.com/28588878/151394973-c05b9c2c-c3b2-416a-8270-77afd79f484c.jpg) |

# Stitched images 2:
| Type | Images |
| --- | --- |
| Stitched image | ![result_26](https://user-images.githubusercontent.com/28588878/143264138-cbb7b009-569b-426e-81f2-14d8eacad415.jpg) |
| Cropped image | ![cropped](https://user-images.githubusercontent.com/28588878/143264182-d472d40c-8b24-4728-8304-42a7cfbbfed8.jpg) |

# Requirements
MATLAB <br />
MATLAB Computer Vision Toolbox <br />
MATLAB Image Processing Toolbox <br />
MATLAB Parallel Computing Toolbox <br />

# Run command
Please use the `Main_AutoPanoStitch.m` to run the program. Change the `folderPath      = '../../../Data/Generic';` to your desired folder path. Also, change the `folderName      = '';` to a valid name.

Change the hyper parameters accordingly if needed. But it is not required though.
```
%% Inputs 2
%--------------------------------------------------------------------------
% Parallel workers
input.numCores = 4;            % Number of cores for parallel processing

%% Inputs 3
% Warping
input.warpType = 'spherical';   % 'spherical' | 'cylindrical' | 'planar' (projective)
input.focalLength = 1200;       % focal length of camera in pixels
input.k1 = 0;                   % Radial distortion coefficient
input.k2 = 0;                   % Radial distortion coefficient
input.k3 = 0;                   % Radial distortion coefficient

% Feature matching
input.detector = 'SIFT';                % 'HARRIS' | 'SIFT' | 'FAST' | 'SURF' | 'BRISK' | 'ORB' | 'KAZE'
input.Matchingmethod = 'Approximate';   %'Exhaustive' (default) | 'Approximate'
input.Matchingthreshold = 1.5;  %       10.0 or 1.0 (default) | percent value in the range (0, 100] | depends on binary and non-binary features
input.Ratiothreshold = 0.6;             % ratio in the range (0,1]

% Image matching (RANSAC)
input.Inliersconfidence = 99.9;         % Inlier confidence [0,100]
input.maxIter = 2000;                   % RANSAC maximum iterations
input.Transformationtype = 'rigid'; % 'rigid' | 'similarity' | 'affine' | 'projective'
input.MaxDistance = 1.50;               % Maximum distance (pixels)

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
```

# Note
Currently, spherical, cyclindrical and planar projections stitching is supported in this version and can recognize multiple panoramas. This work is in progress, further improvements such as inlcusion of a full view 360 x 180 degree panoramas stitching (everything visible from a point), automatic panorama straightening, runtime speed optimization and Graphical User Interface (GUI) are under development. Depending on how your images are captured, panaroma being a `spherical`, `cylindrical` or `planar` should be selected judicially using the `input.warpType` and `input.Transformationtype`. Generally, `spherical` with `rigid` transformation should work well. Your patience will be appreciated.

# Known issues
1. Bundle adjustment using the `lsqnonlin` in MATLAB is computationally slow.
2. Gain compensation on large number of images is computationally slow.
3. Due to the rigid transformation and no straightening, some of the panaromas will be skewed.

# Adaptation of open source 
Some of the MATLAB functions are adapted from the [Kevin Luo's GitHub Repo](https://github.com/kluo8128/cs231_project) and heavily improved.

# Licensing conditions
This software is being made available for research purposes only. Please cite the relevant citations as provided in the main file.

# MATLAB Central
[![View Automatic panorama stitcher (AutoPanoStitch) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/105850-automatic-panorama-stitcher-autopanostitch)
