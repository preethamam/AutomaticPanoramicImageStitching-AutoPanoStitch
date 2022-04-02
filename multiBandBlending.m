%%***********************************************************************%
%*                   Automatic panorama stitching                       *%
%*                        Multiband blending                            *%
%*                                                                      *%
%* Code author: Preetham Manjunatha                                     *%
%* Github link: https://github.com/preethamam                           *%
%* Date: 01/27/2022                                                     *%
%************************************************************************%
% [im] = multiBandBlending(input, warpedImages, warpedWeights)
%
% Given a set of images with the same size, applies multi-band blending
% with a m-level Laplacian pyramid. Returns the resulting blended image.
%
% Credits: Adapted from online code by Hao Jiang, Boston College at
% http://www.phototalks.idv.tw/academic/?p=1275
% Modified and improved by Dr. Preetham Manjunatha
%
%

function [im] = multiBandBlending(input, warpedImages, warpedWeights) % assumes all images same size
m = input.bands;
k = length(warpedImages);
gauss = cell(k, m);
weights = cell(k, m);
laplace = cell(k, m);
blend = cell(1, m);

% Gaussian filter
h = fspecial('gaussian', input.filtSize, input.MBBsigma);

% Create Wmax weights/images
Wmax = cell(1,length(warpedWeights));

wsum = cellfun(@(x) sum(x,3), warpedWeights, 'UniformOutput',false);
wsum2 = cat(3,wsum{:});
wsum3 = sum(wsum2,3);
zeroInd = wsum3 == 0;

[~, index] = max(wsum2,[],3,'omitnan');
maxInd = index .* imcomplement(zeroInd);

for i = 1:length(warpedWeights)    
    Wmaxind = maxInd == i;
    WmaxTemp    = zeros(size(warpedWeights{i}));

    % Reshape the mask
    WmaxMask = repmat(Wmaxind,1,1,3);
    
    % Fill ones
    WmaxTemp(WmaxMask) = 1;

    % Wmax final
    Wmax{i} = WmaxTemp;
end

% Gaussian pyramid
for i = 1:k
    gauss{i,1} = double(warpedImages{i}) ./ 255;
    weights{i,1} = Wmax{i};
    for j = 2:m
        gauss{i,j} = imresize(imfilter(gauss{i,j-1}, h), 0.5);
        weights{i,j} = imresize(imfilter(weights{i,j-1}, h), 0.5, 'bilinear');
    end
end

% Laplacian pyramid
for i = 1:k
   for j = 1:m-1
       h = size(gauss{i,j}, 1);
       w = size(gauss{i,j}, 2);
       laplace{i,j} = gauss{i,j} - imresize(gauss{i,j+1}, [h, w]);
   end
   laplace{i,m} = gauss{i,m};
end

% Multi-band blending
for j = 1:m
    denominator = zeros(size(laplace{1,j}));
    blend{j} = zeros(size(laplace{1,j}));
    for i = 1:k
        blend{j} = blend{j} + laplace{i,j} .* weights{i,j};
        denominator = denominator + weights{i,j};
    end
    denominator(denominator == 0) = Inf;
    blend{j} = blend{j} ./ denominator;
end

% Laplacian pyramid reconstruction
im = blend{m};
for i = 1:m-1
    j = m - i;
    h = size(blend{j}, 1);
    w = size(blend{j}, 2);

    im = blend{j} + imresize(im, [h, w]);    
end
im = uint8(255 .* im);

end
