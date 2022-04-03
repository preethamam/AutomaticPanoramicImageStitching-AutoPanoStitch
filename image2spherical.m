function imageSpherical = image2spherical(image, f, k1, k2, k3)

%%***********************************************************************%
%*                   Image to spherical projection                      *%
%*              Projects normal image to a spherical warp               *%
%*                                                                      *%
%* Code author: Preetham Manjunatha                                     *%
%* Github link: https://github.com/preethamam
%* Date: 08/02/2021                                                     *%
%************************************************************************%
%
%************************************************************************%
%
% Usage: imageSpherical = image2spherical(image, f, k1, k2, k3)
% Inputs: image - input image
%         f  - focal length in pixels (typically varies from 200 t0 800)
%              highly depends on the camera.
%         k1 - Radial distortion coefficient.
%         k2 - Radial distortion coefficient.
%         k3 - Radial distortion coefficient.
%
% Outputs: imageSpherical - Warpped image to spherical coordinates

% Get image size
[ydim, xdim, bypixs] = size(image);

% Initialize an array
imageSpherical = uint8(zeros(ydim, xdim, bypixs));

% Get the center of image
xc = round(xdim/2);
yc = round(ydim/2);

% Create X and Y coordinates grid
[X,Y] = meshgrid(1:xdim, 1:ydim);

% Perform the cylindrical projection
theta = (X - xc)/f;
phi   = (Y - yc)/f;

% Spherical coordinates to Cartesian
xcap = sin(theta) .* cos(phi);
ycap = sin(phi);
zcap = cos(theta) .* cos(phi);

xn = xcap ./ zcap;
yn = ycap ./ zcap;
r = xn.^2 + yn.^2;

% Lens distortion addition 
xd = xn .* (1 + k1 * r.^2 + k2 * r.^4 + k3 * r.^6);
yd = yn .* (1 + k1 * r.^2 + k2 * r.^4 + k3 * r.^6);

% Convert to floor
ximg = floor(f * xd + xc);
yimg = floor(f * yd + yc);

% Find boundary of the cylindrical projection
mask = ximg > 0 & ximg <= xdim & yimg > 0 & yimg <= ydim;

% Reshape the mask
if bypixs == 3
    mask = repmat(mask,1,1,3);
end

% Get projections
imageSpherical(mask) = image(mask);

end

