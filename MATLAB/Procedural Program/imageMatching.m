function [allMatches, numMatches, tforms] = imageMatching(input, n, keypoints, matchesAll, images)

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
    matSize = size(allMatches);

    % Use symmetry and run for upper triangular matrix
    linearIdxs = reshape(1:numel(allMatches), size(allMatches));
    IuptriIdx = nonzeros(triu(linearIdxs,1));
    IlowtriIdx = nonzeros(triu(linearIdxs',1));

    % Initialize
    allMatches_temp = cell(1,length(IuptriIdx));
    numMatches_temp = zeros(1,length(IuptriIdx));

    % Initialize transformation matrix
    switch input.Transformationtype
        case 'rigid' 
            tforms(n,n) = rigidtform2d(eye(3));
            tforms_ij_ji(2,length(IuptriIdx)) = rigidtform2d(eye(3));
        case 'similarity' 
            tforms(n,n) = simtform2d(eye(3));
            tforms_ij_ji(2,length(IuptriIdx)) = simtform2d(eye(3));
        case 'affine' 
            tforms(n,n) = affinetform2d(eye(3));
            tforms_ij_ji(2,length(IuptriIdx)) = affinetform2d(eye(3));
        case 'projective'
            tforms(n,n) = projtform2d(eye(3));
            tforms_ij_ji(2,length(IuptriIdx)) = projtform2d(eye(3));
        case 'translation'
            tforms(n,n) = transltform2d(eye(3));
            tforms_ij_ji(2,length(IuptriIdx)) = transltform2d(eye(3));
        otherwise
            error('Requires a transformation type!')
    end

    % Minimum number of features 
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

    % Match images
    parfor i = 1:length(IuptriIdx) 

        % IND2SUB converts from a "linear" index into individual
        % subscripts
        [ii,jj] = ind2sub(matSize, IuptriIdx(i));        

        % Keypoints matches
        matches = matchesAll{ii,jj};

        % Number of features
        nf = size(matches, 2);

        % Image matching
        % Filter matches using RANSAC (model maps keypt i to keypt j)
        if nf >= nfmin
            [inliers, model, status] = refineMatch(input, keypoints{ii}, keypoints{jj}, matches, ...
                                                   images{ii}, images{jj});

            % Number of inliers
            ni = length(inliers);

            % Verify image matches using probabilistic model
            if strcmp(input.warpType,'spherical') || strcmp(input.warpType,'cylindrical')
                if ni > 5 + 0.15 * nf %2 % accept as correct image match
                    allMatches_temp{i} = matches(:,inliers);
                    numMatches_temp(i) = ni;
                    model_temp = [model; model.invert];
                    tforms_ij_ji(:,i) = model_temp;
                end
            elseif strcmp(input.warpType,'planar') && strcmp(input.Transformationtype,'rigid') || ...
                   strcmp(input.Transformationtype, 'similarity')
                if ni > 5 + 0.025 * nf %2 % accept as correct image match
                    allMatches_temp{i} = matches(:,inliers);
                    numMatches_temp(i) = ni;
                    model_temp = [model; model.invert];
                    tforms_ij_ji(:,i) = model_temp;
                end
            elseif strcmp(input.Transformationtype, 'affine')
                 if ni > 5 + 0.15 * nf %3 % accept as correct image match
                    allMatches_temp{i} = matches(:,inliers);
                    numMatches_temp(i) = ni;
                    model_temp = [model; model.invert];
                    tforms_ij_ji(:,i) = model_temp;
                end
            else
                if ni > 8 + 0.3 * nf %5.9 + 0.22 * nf % accept as correct image match
                    allMatches_temp{i} = matches(:,inliers);
                    numMatches_temp(i) = ni;
                    model_temp = [model; model.invert];
                    tforms_ij_ji(:,i) = model_temp;
                end
            end
        end
    end

    % Populate all matches symmetric matrix
    allMatches(IuptriIdx) = allMatches_temp;

    % Populate number of matches symmetric matrix
    numMatches(IuptriIdx) = numMatches_temp;

    % Populate transformation matrix
    tforms(IuptriIdx) = tforms_ij_ji(1,:);     
    tforms(IlowtriIdx) = tforms_ij_ji(2,:);
end

%--------------------------------------------------------------------------------------------------------
% Auxillary functions
%--------------------------------------------------------------------------------------------------------
% [inliers, model] = refineMatch(P1, P2, matches)
%
% Returns the set of inliers and corresponding homography that maps matched
% points from P1 to P2.
function [inliers, model, status] = refineMatch(input, P1, P2, matches, image1, image2)
    % Matched points
    matchedPts_1 = P1(:,matches(1,:))';
    matchedPts_2 = P2(:,matches(2,:))';

    % Estimate the transformation matrix
    [model, inliers, status] = estgeotform2d(matchedPts_2, matchedPts_1, input.Transformationtype, ...
                            'Confidence', input.Inliersconfidence, 'MaxNumTrials', input.maxIter, ...
                            'MaxDistance', input.MaxDistance);

    % Find inliers
    inliers = find(inliers);

    % Plot of the feature matches
    if(input.showKeypointsPlot == 1)
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

