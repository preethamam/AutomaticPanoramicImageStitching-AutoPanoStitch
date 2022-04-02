function [allMatches, numMatches, tforms] = imageMatching(input, n, images, keypoints, allDescriptors)

%%***********************************************************************%
%*                   Automatic panorama stitching                       *%
%*                        Image matching                                *%
%*                                                                      *%
%* Code author: Preetham Manjunatha                                     *%
%* Github link: https://github.com/preethamam                           *%
%* Date: 01/27/2022                                                     *%
%************************************************************************%

% Initialize
allMatches = cell(n);
numMatches = zeros(n);
switch input.Transformationtype
    case 'rigid' 
        tforms(n,n) = rigid2d(eye(3));
    case 'similarity' 
        tforms(n,n) = affine2d(eye(3));
    case 'affine' 
        tforms(n,n) = affine2d(eye(3));
    case 'projective'
        tforms(n,n) = projective2d(eye(3));        
end
iterations = size(allMatches);
 
if strcmp(input.warpType,'spherical') || strcmp(input.warpType,'cylindrical')
    nfmin = 2;
elseif (strcmp(input.warpType,'planar') && strcmp(input.Transformationtype,'rigid')) || ...
        strcmp(input.Transformationtype,'similarity')
    nfmin = 2;
elseif strcmp(input.Transformationtype, 'affine')
    nfmin = 3;
else
    nfmin = 4;
end


parfor i = 1:prod(iterations) 
    % IND2SUB converts from a "linear" index into individual
    % subscripts
    [ii,jj] = ind2sub(iterations, i);
    
    if (ii~=jj)
        matches = getMatches(input, allDescriptors{ii}, allDescriptors{jj});
        nf = size(matches, 2);

        % Image matching
        % Filter matches using RANSAC (model maps keypt i to keypt j)
        if nf >= nfmin
            [inliers, model, status] = refineMatch(input, keypoints{ii}, keypoints{jj}, matches, ...
                                                   images{ii}, images{jj});
            ni = length(inliers);

            % Verify image matches using probabilistic model
            if strcmp(input.warpType,'spherical') || strcmp(input.warpType,'cylindrical')
                if ni > 2 % accept as correct image match
                    allMatches{i} = matches(:,inliers);
                    numMatches(i) = ni;
                    tforms(i) = model;
                end

            elseif strcmp(input.warpType,'planar') && strcmp(input.Transformationtype,'rigid') || ...
                   strcmp(input.Transformationtype, 'similarity')
                if ni > 2 % accept as correct image match
                    allMatches{i} = matches(:,inliers);
                    numMatches(i) = ni;
                    tforms(i) = model;
                end
            elseif strcmp(input.Transformationtype, 'affine')
                 if ni > 3 % accept as correct image match
                    allMatches{i} = matches(:,inliers);
                    numMatches(i) = ni;
                    tforms(i) = model;
                end
            else
                if ni > 8 + 0.3 * nf %5.9 + 0.22 * nf % accept as correct image match
                    allMatches{i} = matches(:,inliers);
                    numMatches(i) = ni;
                    tforms(i)     = model;
                end
            end

        end
    end
end

% Clean up the lower triangular part
iii = ones(size(allMatches));
idx = find(tril(iii,-1));
allMatches(idx) = {[]};
numMatches = triu(numMatches,0);

end

%--------------------------------------------------------------------------------------------------------
% Auxillary functions
%--------------------------------------------------------------------------------------------------------
% [matches] = getMatches(features1,features2)
function matches = getMatches(input,features1,features2)

    matches = matchFeatures(features1,features2, 'Method', input.Matchingmethod, ...
                            'MatchThreshold', input.Matchingthreshold, ...
                            'MaxRatio', input.Ratiothreshold);
    
    matches = double(matches)';
end

% [inliers, model] = refineMatch(P1, P2, matches)
%
% Returns the set of inliers and corresponding homography that maps matched
% points from P1 to P2.
%
% Adapted from Problem 3 of PS3
function [inliers, model, status] = refineMatch(input, P1, P2, matches, image1, image2)

matchedPts_1 = P1(:,matches(1,:))';
matchedPts_2 = P2(:,matches(2,:))';


[model, inliers, status] = estimateGeometricTransform2D(matchedPts_2, matchedPts_1, input.Transformationtype, ...
                        'Confidence', input.Inliersconfidence, 'MaxNumTrials', input.maxIter, ...
                        'MaxDistance', input.MaxDistance);

inliers = find(inliers);

% Plot of the feature matches
if(input.showPlot == 1)
    figure(1);
    subplot(1,3,1)
    montage({image1, image2});
    title('Original Images')
    
    subplot(1,3,2)
    showMatchedFeatures(image1, image2, matchedPts_1, matchedPts_2, 'montage')
    title('Putative Matches')
    
    subplot(1,3,3)
    showMatchedFeatures(image1, image2, matchedPts_1(inliers,:), matchedPts_2(inliers,:), 'montage')
    title('Inlier Matches')
end

end

