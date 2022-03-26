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
function [im] = multiBandBlending(input, warpedImages, warpedWeights) % assumes all images same size
m = input.bands;
k = length(warpedImages);
gauss = cell(k, m);
weights = cell(k, m);
laplace = cell(k, m);
blend = cell(1, m);

% Gaussian pyramid
for i = 1:k
    gauss{i,1} = double(warpedImages{i}) ./ 255;
    weights{i,1} = warpedWeights{i};
    for j = 2:m
        gauss{i,j} = imresize(gauss{i,j-1}, 0.5);
        weights{i,j} = imresize(weights{i,j-1}, 0.5, 'bilinear');
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
