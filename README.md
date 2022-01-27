# AutoPanoStitch
Automatic Panorama Stitching software in MATLAB.

# Sample 1: White background
| Type | Images |
| --- | --- |
| Stitched image | ![out_6_NSH](https://user-images.githubusercontent.com/28588878/143262116-10a768b1-d791-4758-86a9-2d2b906e8644.jpg) |
| Cropped image | ![cropped1](https://user-images.githubusercontent.com/28588878/143262859-213860cb-1e2f-4986-9c2d-e1bec2d368a2.jpg) |

# Sample 2: Black background
| Type | Images |
| --- | --- |
| Stitched image | ![result_26](https://user-images.githubusercontent.com/28588878/143264138-cbb7b009-569b-426e-81f2-14d8eacad415.jpg) |
| Cropped image | ![cropped](https://user-images.githubusercontent.com/28588878/143264182-d472d40c-8b24-4728-8304-42a7cfbbfed8.jpg) |

# Requirements
MATLAB <br />
MATLAB Computer Vision Toolbox <br />
MATLAB Image Processing Toolbox <br />
MATLAB Parallel Computing Toolbox <br />
[VLFeat](https://www.vlfeat.org/install-matlab.html)

# Run command
Please use the `Main_AutoPanoStitch.m` to run the program. Change the `folderPath      = '../../../Data/Generic';` to your desired folder path. Also, change the `folderName      = '';` to a valid name.

# Note
Currently, only planar projections stitching is supported in this version and can recognize multiple panoramas. This work is in progrees further improvements such as inlcusion of spherical, cylindrical projections, full view 360 x 180 degree panoramas stitching (everything visible from a point), automatic panorama straightening, runtime speed optimization and Graphical User Interface (GUI) are under development. Your patience will be appreciated.

# Known issues
1. Bundle adjustment using the `lsqnonlin` in MATLAB is computationally slow.
2. Gain compensation on large number of images is computationally slow.

# Adaptation of open source 
Some of the MATLAB functions are adapted from the [Kevin Luo's GitHub Repo](https://github.com/kluo8128/cs231_project) and heavily improved.

# Licensing conditions
This software is being made available for research purposes only. Please cite the relevant citations as provided in the main file.

# MATLAB Central
