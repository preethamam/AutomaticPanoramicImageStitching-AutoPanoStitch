# AutoPanoStitch
[![View Automatic panorama stitcher (AutoPanoStitch) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/105850-automatic-panorama-stitcher-autopanostitch) [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=preethamam/AutomaticPanoramicImageStitching-AutoPanoStitch)

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
MATLAB Optimization Toolbox

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
input.warpType = 'cylindrical';   % 'spherical' | 'cylindrical' | 'planar' (projective)
input.focalLength = 1200;       % focal length of camera in pixels
input.k1 = 0;                   % Radial distortion coefficient
input.k2 = 0;                   % Radial distortion coefficient
input.k3 = 0;                   % Radial distortion coefficient

% Feature matching
input.detector = 'SIFT';                % 'HARRIS' | 'SIFT' | 'FAST' | 'SURF' | 'BRISK' | 'ORB' | 'KAZE'
input.Matchingmethod = 'Approximate';   %'Exhaustive' (default) | 'Approximate'
input.Matchingthreshold = 1.5;  %       10.0 or 1.0 (default) | percent value in the range (0, 100] | depends on binary and non-binary features
input.Ratiothreshold = 0.6;             % ratio in the range (0,1]
input.NumLayersInOctave = 3;            % Number of layers in each octave -- SIFT only
input.ContrastThreshold = 0.00133;      % Contrast threshold for selecting the strongest features, 
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
```

# Note
Depending on how your images are captured, panaroma being a `spherical`, `cylindrical` or `planar` should be selected judicially using the `input.warpType` and `input.Transformationtype`. Generally, `spherical` or `cylindrical` projections with `affine` or `rigid` transformation should work well in most of the cases. However, some panoramas specifically looks good in `projective` transformation, e.g. flatbed scanner or whiteboard (`affine` works well too) images.

Currently, spherical, cyclindrical and planar projections stitching is supported in this version and can recognize multiple panoramas. This work is in progress, further improvements such as inlcusion of a full view 360 x 180 degree panoramas stitching (everything visible from a point), automatic panorama straightening, runtime speed optimization and Graphical User Interface (GUI) are under development. Your patience will be appreciated.

# Known issues
1. Bundle adjustment using the `lsqnonlin` in MATLAB is computationally slow.
2. Gain compensation on large number of images is computationally slow.
3. Due to the rigid transformation and no straightening, some of the panoramas will be skewed.
4. Gain compensation values for some images are same.

# Adaptation of open source 
Bundle adjustment functions in MATLAB are adapted from the [Kevin Luo's GitHub Repo](https://github.com/kluo8128/cs231_project) and heavily improved.

# Licensing conditions
The original implementation of the automatic panaroma stitching by Dr. Brown was written in C++ and is `LICENSED under The University of British Columbia`. This is being programmed and made available to public for academic / research purposes only. Please cite the relevant citations as provided in the main file.

Please rate and provide feedback for the further improvements.
