%%***********************************************************************%
%*                   Automatic panorama stitching                       *%
%*                        Panorama render                               *%
%*                                                                      *%
%* Code author: Preetham Manjunatha                                     *%
%* Github link: https://github.com/preethamam                           *%
%* Date: 01/27/2022                                                     *%
%************************************************************************%

function [panorama, gainpanorama, noBlendCompensationPanorama, gainImages, ...
            gainRGB, warpedImages, xCorrect, yCorrect] = ...
            renderPanorama(input, images, focalLengths, tforms, ccs, cc)
n = length(tforms);

% If the total area is too large, rescale the images to be smaller
if input.resizeStitchedImage
    [height, width] = getPanoramaSize(images, tforms, ccs, cc);
    area = height * width;
    maxArea = 3e6;
    if area > maxArea
        f = sqrt(area / maxArea);
        S_inv = inv(diag([f; f; 1]));
        for i = 1:n
            tf = tforms(i).T * S_inv; 
            tf(1:2,3) = 0;
            tf(3,3) = 1;

            tforms(i).T = single(tf);
        end
    end
end

xlim = zeros(n,2);
ylim = zeros(n,2);
hMax = 0;
wMax = 0;
indices = find(ccs == cc);

if numel(indices) == 1
    panorama = images{indices};
    gainpanorama = images{indices};
    noBlendCompensationPanorama = [];
    gainImages = [];
    gainRGB = []; 
    warpedImages = [];
    xCorrect = []; 
    yCorrect = [];

    warning('Single image detected. No panorama created!')
    
else
    for index = 1:length(indices)
        i = indices(index);
        h = size(images{i}, 1);
        w = size(images{i}, 2);
        hMax = max(h, hMax);
        wMax = max(w, wMax);
        [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1, w], [1, h]);
    end
    
    % Find the minimum and maximum output limits
    xMin = min(xlim(indices,1));
    xMax = max(xlim(indices,2));
    yMin = min(ylim(indices,1));
    yMax = max(ylim(indices,2));
    
    % Width and height of panorama
    width  = round(xMax - xMin);
    height = round(yMax - yMin);
    
    % Create a 2-D spatial reference object defining the size of the panorama.
    xLimits = [xMin, xMax];
    yLimits = [yMin, yMax];
    panoramaView = imref2d([height, width], xLimits, yLimits);
    
    k = length(indices);
    warpedImages = cell(k, 1);
    warpedWeights = cell(k, 1);
    
    xCorrect = zeros(5,k);
    yCorrect = zeros(5,k);
    
    % Initial no blend panaroma
    for index = 1:k
        i = indices(index);
        
        image = images{i};
        [imRows, imCols, imChannel] = size(image);
        
        % Get warped image corners 
        imCorners = [1,1; 1,imCols; imRows, imCols; imRows,1; 1,1];
        [x,y] = transformPointsForward(tforms(i), imCorners(:,2), imCorners(:,1));
        xCorrect(:,i) = x - panoramaView.XWorldLimits(1); 
        yCorrect(:,i) = y - panoramaView.YWorldLimits(1);
        
        warpedImage = imwarp(image, tforms(i), 'OutputView', panoramaView, 'FillValues', 0);
        
        switch input.blending
            case 'none'
                warpedImages{index} = warpedImage;
                
            case {'linear', 'multiband'}
                warpedImages{index} = warpedImage;
                weight = getWeight(size(image));
                warpedWeight = imwarp(weight, tforms(i), 'OutputView', panoramaView);
                warpedWeights{index} = warpedWeight;  
            otherwise
        end
      
    end
    
    switch input.blending
        case 'none'
            % Panorama creation
            warpedImagesMask = cellfun(@(x) repmat(imfill(imbinarize(rgb2gray(255 * x)), 'holes'), 1, 1, size(warpedImages{1},3)), ...
                                        warpedImages, 'UniformOutput',false);
            panorama = zeros(size(warpedImages{1}), 'uint8');
            for i = 1:length(warpedImages)
                panorama(warpedImagesMask{i}) = warpedImages{i}(warpedImagesMask{i});
            end
            
            % Gain compensation null values
            gainpanorama = []; 
            gainImages = [];  
            gainRGB = []; 
            noBlendCompensationPanorama = [];
    
        case 'linear'
            tic
            % Gain compensation
            [gainpanorama, gainImages, gainRGB] = gainCompensation(input, warpedImages);
            fprintf('Gain compensation: %f seconds\n', toc);

            % Apply linear blending
            tic
            panorama = linearBlending(input, gainImages, warpedWeights);
            fprintf('Linear blending: %f seconds\n', toc);

            % No compensation or blend panorama creation
            warpedImagesMask = cellfun(@(x) repmat(imfill(imbinarize(rgb2gray(255 * x)), 'holes'), 1, 1, size(warpedImages{1},3)), ...
                                        warpedImages, 'UniformOutput',false);
            noBlendCompensationPanorama = zeros(size(warpedImages{1}), 'uint8');
            for i = 1:length(warpedImages)
                noBlendCompensationPanorama(warpedImagesMask{i}) = warpedImages{i}(warpedImagesMask{i});
            end
    
        case 'multiband'
            tic
            % Gain compensation
            [gainpanorama, gainImages, gainRGB] = gainCompensation(input, warpedImages);
            fprintf('Gain compensation: %f seconds\n', toc);
    
            % Apply multi-band blending
            tic
            panorama = multiBandBlending(input, gainImages, warpedWeights);
            fprintf('Multi-band blending: %f seconds\n', toc);

            % No compensation or blend panorama creation
            warpedImagesMask = cellfun(@(x) repmat(imfill(imbinarize(rgb2gray(255 * x)), 'holes'), 1, 1, size(warpedImages{1},3)), ...
                                        warpedImages, 'UniformOutput',false);
            noBlendCompensationPanorama = zeros(size(warpedImages{1}), 'uint8');
            for i = 1:length(warpedImages)
                noBlendCompensationPanorama(warpedImagesMask{i}) = warpedImages{i}(warpedImagesMask{i});
            end

        otherwise
            error('Please select blending mode.')
    end
end
end

%--------------------------------------------------------------------------------------------------------
% Auxillary functions
%--------------------------------------------------------------------------------------------------------
% [weight] = getWeight(size)
%
% Given the size of an RGB image, returns a weight matrix with the same
% dimensions, w(x, y) = w(x) * w(y), where w(x) and w(y) vary linearly from
% 1 at the center of the image to 0 at the edges.
function [weight] = getWeight(size)

    h = size(1);
    w = size(2);
    c = size(3);

    wx = ones(1, w);
    wx(1:ceil(w/2)) = linspace(0, 1, ceil(w/2));
    wx(floor(w/2 + 1):w) = linspace(1, 0, w - floor(w/2));
    wx = repmat(wx, h, 1, c);
    wy = ones(h, 1);
    wy(1:ceil(h/2)) = linspace(0, 1, ceil(h/2));
    wy(floor(h/2 + 1):h) = linspace(1, 0, h - floor(h/2));
    wy = repmat(wy, 1, w, c);
    weight = wx .* wy;
end

% [height, width] = getPanoramaSize(images, tforms, ccs, cc)
%
% Returns the size of the panorama from applying projective transformations
% on the images in the connected component with index cc.
%
% Credits: Adapted from online MATLAB example "Feature Based Panoramic
% Image Stitching" at
% http://www.mathworks.com/examples/matlab-computer-vision/mw/vision_product-FeatureBasedPanoramicImageStitchingExample-feature-based-panoramic-image-stitching
function [height, width] = getPanoramaSize(images, tforms, ccs, cc)
n = length(tforms);
xlim = zeros(n,2);
ylim = zeros(n,2);
hMax = 0;
wMax = 0;

indices = find(ccs == cc);

for index = 1:length(indices)
    i = indices(index);
    h = size(images{i}, 1);
    w = size(images{i}, 2);
    hMax = max(h, hMax);
    wMax = max(w, wMax);
    [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1, w], [1, h]);
end

% Find the minimum and maximum output limits
xMin = min(xlim(indices,1));
xMax = max(xlim(indices,2));
yMin = min(ylim(indices,1));
yMax = max(ylim(indices,2));

% Width and height of panorama
width = round(xMax - xMin);
height = round(yMax - yMin);
end
