%%***********************************************************************%
%*                   Automatic panorama stitching                       *%
%*                        Linear blending                               *%
%*                                                                      *%
%* Code author: Preetham Manjunatha                                     *%
%* Github link: https://github.com/preethamam                           *%
%* Date: 01/27/2022                                                     *%
%************************************************************************%
% [im] = linearBlending(input, warpedImages, warpedWeights)
%
% Given a set of images with the same size, applies linear blending
% Returns the resulting blended image.
%
function [im] = linearBlending(input, warpedImages, warpedWeights) % assumes all images same size

% Linear blending
numerator = zeros(size(warpedImages{1}));
denominator = zeros(size(warpedImages{1}));

for i = 1:length(warpedImages)
    numerator = numerator + (double(warpedImages{i}) .* warpedWeights{i});
    denominator = denominator + warpedWeights{i};
end

im = uint8(numerator./denominator);

end
