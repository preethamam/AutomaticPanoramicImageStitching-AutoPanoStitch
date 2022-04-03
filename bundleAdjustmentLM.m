function [finalPanoramaTforms, concomps] = bundleAdjustmentLM(input, images, keypoints, allMatches, numMatches, initialTforms)

%%***********************************************************************%
%*                   Automatic panorama stitching                       *%
%*                        Bundle adjustment                             *%
%*                                                                      *%
%* Code author: Preetham Manjunatha                                     *%
%* Github link: https://github.com/preethamam                           *%
%* Date: 01/27/2022                                                     *%
%************************************************************************%

% Find connected components of image matches
numMatchesG = graph(numMatches,'upper');
[concomps, ccBinSizes] = conncomp(numMatchesG);
panaromaCCs = find(ccBinSizes>=1);
ccnum = numel(panaromaCCs);
[tree] = getMST(numMatches);

finalPanoramaTforms = cell(1,ccnum);

% Find panoramas 
options = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt', ...
    'FunctionTolerance', 1e-6, 'StepTolerance', 1e-6, 'Display','off', 'UseParallel', true);

for cc = 1:ccnum
    indices = find(concomps == cc);
    k = length(indices);
        
    % Find the center image of the panorama that minimizes total area      
    areas = zeros(k, 1);
    for index = 1:k
        i = indices(index);
        finalTforms = getTforms(input, tree, i, initialTforms);

        [height, width] = getPanoramaSize(images, finalTforms, concomps, cc);

        areas(index) = height * width;
    end
    [~, index] = min(areas);
    center = indices(index);
    
    finalTforms = getTforms(input, tree, center, initialTforms);

    if k < 2 % Skip if only one image in connected component
        finalPanoramaTforms{cc} = finalTforms;
        continue;
    end
    
    % Get an ordering of the images in the panorama
    ordering = getOrdering(indices, tree);
       
    % Minimize the sum of squared projection error over all matches
    Hs_initial = zeros(9, k);
    Hs_LMfinal = zeros(9, k);
    for index = 1:k
        i = indices(index);
        H = finalTforms(i).T';
        Hs_initial(:,index) = reshape(H, [], 1);
    end

    for c = 2:k
        subordering = ordering(1:c);
        Phi = zeros(9, c);
        for j = 1:c
            Phi(:,j) = Hs_initial(:,subordering(j));
        end
        
        Phi = reshape(Phi(1:8,:), [], 1); % last entry of each H is 1
        
        fxn = @(Phi)projectionError(Phi, indices(subordering), ...
            keypoints, allMatches, numMatches);
        
        [Phi,resnorm,residual,exitflag,LMoutput] = lsqnonlin(fxn, Phi, [], [], options);
        
        Phi = [reshape(Phi, [], c); ones(1, c)];
        for j = 1:c
            Hs_LMfinal(:,subordering(j)) = Phi(:,j);
        end
    end

    % Update transforrmations
    refinedTforms = finalTforms;
    for index = 1:k
        i = indices(index);
        H = reshape(Hs_LMfinal(:,index), [], 3);
        if strcmp(input.warpType,'spherical') || strcmp(input.warpType,'cylindrical')
            tf = H';
            tf(1:2,3) = 0;
            refinedTforms(i).T = single(tf);
        elseif strcmp(input.warpType,'planar') && (strcmp(input.Transformationtype,'rigid') ...
                || strcmp(input.Transformationtype,'similarity') || ...
                   strcmp(input.Transformationtype,'affine'))
            tf = H';
            tf(1:2,3) = 0;
            refinedTforms(i).T = single(tf);
        else
            refinedTforms(i).T = H';
        end
    end   

    finalPanoramaTforms{cc} = refinedTforms;
end

end

%--------------------------------------------------------------------------------------------------------
% Auxillary functions
%--------------------------------------------------------------------------------------------------------

%--------------------------------------------------------------------------------------------------------
% [tree] = getMST(G)
%--------------------------------------------------------------------------------------------------------
%
% Returns the adjacency matrix of the maximum spanning tree of the
% undirected graph with weighted adjacency matrix G.
function [tree] = getMST(G)
n = size(G, 1);
ccs = (1:n)'; % list of component number of each vertex
components = cell(n, 1);
for i = 1:n
    components{i} = i;
end
tree = zeros(n);
numEdges = 0;
[values, indices] = sort(G(:), 'descend');
for k = 1:length(indices)
    if values(k) > 0
        i = mod(indices(k) - 1, n) + 1;
        j = ceil(indices(k) / n);
        if ccs(i) ~= ccs(j)
            tree(i,j) = values(k);
            tree(j,i) = values(k);
            components{ccs(i)} = [components{ccs(i)}; components{ccs(j)}];
            ccs(components{ccs(j)}) = ccs(i);
            numEdges = numEdges + 1;
        end
    end
    if numEdges == n - 1
        break;
    end
end
end

%--------------------------------------------------------------------------------------------------------
% [ordering] = getOrdering(indices, tree)
%--------------------------------------------------------------------------------------------------------
%
% Given the adjacency matrix for a weighted forest (tree), finds an
% ordering on the specified indices (a single tree in the forest) that
% greedily maximizes the cumulative weight, i.e., starts with vertices with
% the highest edge weight, then expands outward along tree, adding the
% vertex sharing the highest edge weight in the fringe each time, until all
% vertices in the tree have been added. Returns a permutation of 1 to k
% corresponding to the indices of the ordering, where k is the number of
% vertices in the tree.
function [ordering] = getOrdering(indices, tree)
k = length(indices);
ordering = zeros(k, 1);
visited = zeros(k, 1);
subtree = tree(indices,indices);
edges = getEdges(subtree);
[~, index] = max(edges(3,:));
index_i = edges(1,index);
index_j = edges(2,index);
ordering(1) = index_i;
ordering(2) = index_j;
visited(index_i) = 1;
visited(index_j) = 1;
c = 2;

fringe = [];
for index = 1:k
    if subtree(index,index_j) > 0 && ~visited(index)
        fringe = [fringe, [index; index_j; subtree(index,index_j)]];
    end
end
while c < k
    for index = 1:k
        if subtree(index,index_i) > 0 && ~visited(index)
            fringe = [fringe, [index; index_i; subtree(index,index_i)]];
        end
    end
    [~, index] = max(fringe(3,:));
    index_i = fringe(1,index);
    fringe(:,index) = [];
    c = c + 1;
    ordering(c) = index_i;
    visited(index_i) = 1;
end
end

%--------------------------------------------------------------------------------------------------------
% [edges] = getEdges(G)
%--------------------------------------------------------------------------------------------------------
%
% Returns a list of weighted edges (i, j, w) of the undirected graph with
% adjacency matrix G.
function [edges] = getEdges(G)
n = size(G, 1);
edges = zeros(3, n * (n - 1) / 2);
c = 0;
for i = 1:n
    for j = i + 1:n
        if G(i,j) > 0
            c = c + 1;
            edges(:,c) = [i; j; G(i,j)];
        end
    end
end
edges = edges(:,1:c);
end

%--------------------------------------------------------------------------------------------------------
% [height, width] = getPanoramaSize(images, tforms, ccs, cc)
%--------------------------------------------------------------------------------------------------------
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

%--------------------------------------------------------------------------------------------------------
% [tforms] = getTforms(G, i, keypoints, allMatches)
%--------------------------------------------------------------------------------------------------------
%
% Computes and returns the projective transformations for all images in the
% connected component of image i in the tree with adjacency matrix G, given
% the set of keypoints and matching indices. All tforms are calculated with
% respect to image i.
function [tforms] = getTforms(input, G, i, Tforms)
    n = size(G, 1);
    visited = zeros(n, 1);

    switch input.Transformationtype
        case 'rigid' 
            tforms(n) = rigid2d(eye(3));
        case 'similarity' 
            tforms(n) = affine2d(eye(3));
        case 'affine' 
            tforms(n) = affine2d(eye(3));
        case 'projective'
            tforms(n) = projective2d(eye(3));        
    end
    
    [tforms] = updateTforms(G, i, visited, tforms, Tforms);
end

%--------------------------------------------------------------------------------------------------------
% [tforms] = updateTforms(G, i, visited, tforms, Tforms)
%--------------------------------------------------------------------------------------------------------
%
% Updates and returns the projective transformations for each image j that
% shares an edge with image i in the tree with adjacency matrix G, given
% the corresponding keypoints and matching indices. Recursively updates
% the tforms of the neighbors of each image j.
function [tforms] = updateTforms(G, i, visited, tforms, Tforms)
    n = size(G, 1);
    visited(i) = 1;
    for j = 1:n
        if G(i,j) > 0 && ~visited(j)            
            if i < j
                tform = Tforms(j,i);
                tform = invert(tform); % j->i
            else
                tform = Tforms(i,j);
            end
            tform.T = tform.T * tforms(i).T;
            tforms(j).T = tform.T ./ tform.T(3,3);
            [tforms] = updateTforms(G, j, visited, tforms, Tforms);
        end
    end
end

%--------------------------------------------------------------------------------------------------------
% [error] = projectionError(Phi, indices, keypoints, allMatches, numMatches)
%--------------------------------------------------------------------------------------------------------
%
% Returns the projection error, given the homography values in the vector
% Phi, the keypoints, matching indices, and number of matches. Only
% considers the specified indices when computing the error.
function [error] = projectionError(Phi, indices, keypoints, allMatches, numMatches)
error = 0;
m = length(indices);

% extract homography transformations from Phi
Hs = [reshape(Phi, 8, []); ones(1,m)];

for index_i = 1:m
    i = indices(index_i);
    keypt_i = keypoints{i};
    H_i = reshape(Hs(:,index_i), [], 3); % i->c
    for index_j = index_i + 1:m
        j = indices(index_j);
        if numMatches(i,j) > 0 || numMatches(j,i) > 0
            keypt_j = keypoints{j};
            if i < j
                matches = allMatches{i,j};
            else
                matches = allMatches{j,i};
                matches = [matches(2,:); matches(1,:)];
            end
            
            %matches = allMatches{i,j};
            H_j = reshape(Hs(:,index_j), [], 3); % j->c
            H = H_i \ H_j; % j->i = (c->i) * (j->c)
            H_inv = H_j \ H_i; % i->j = (c->j) * (i->c)
            
            for k = 1:size(matches,2)
                u_i = [keypt_i(1,matches(1,k)); keypt_i(2,matches(1,k)); 1];
                u_j = [keypt_j(1,matches(2,k)); keypt_j(2,matches(2,k)); 1];
                
                % Matching error from image j to i
                p_ij = H * u_j;
                r_ij = (u_i(1:2) ./ u_i(3)) - (p_ij(1:2) ./ p_ij(3));
                error = error + norm(r_ij) ^ 2;
                
                % Matching error from image i to j
                p_ji = H_inv * u_i;
                r_ji = (u_j(1:2) ./ u_j(3)) - (p_ji(1:2) ./ p_ji(3));
                error = error + norm(r_ji) ^ 2;
            end
        end
    end
end
end
