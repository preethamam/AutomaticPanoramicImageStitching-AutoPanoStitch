% [im] = linearBlending(input, warpedImages, warpedWeights)
%
% Given a set of images with the same size, applies multi-band blending
% with a 2-level Laplacian pyramid. Returns the resulting blended image.
%
% Credits: Adapted from online code by Hao Jiang, Boston College at
% http://www.phototalks.idv.tw/academic/?p=1275
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
