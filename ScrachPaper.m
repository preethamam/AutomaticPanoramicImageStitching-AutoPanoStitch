clc; close all; clear;

Iijc = Iij;
iii = ones(size(Iij));
idx = find(tril(iii,-1));

Iijp = Iij';
Iijc(idx) = Iijp(idx);

%     GMR(1,1) = ((Nij(1,2) * Iij{1,2}(1)^2 + Nij(2,1) * Iji{1,2}(1)^2) / sigmaN^2) + ((Nij(1,1) + Nij(1,2)) / sigmag^2);
%     GMR(1,2) = - (Nij(1,2) * Iij{1,2}(1) * Iji{2,1}(1) + Nij(2,1) * Iij{2,1}(1) * Iji{2,1}(1)) / sigmaN^2;
%     GMR(2,1) = GMR(1,2);
%     GMR(2,2) = ((Nij(2,1) * Iij{2,1}(1)^2 + Nij(1,2) * Iji{2,1}(1)^2) / sigmaN^2) + ((Nij(2,1) + Nij(2,2)) / sigmag^2);
% 
%     GMG(1,1) = ((Nij(1,2) * Iij{1,2}(2)^2 + Nij(2,1) * Iji{1,2}(2)^2) / sigmaN^2) + ((Nij(1,1) + Nij(1,2)) / sigmag^2);
%     GMG(1,2) = - (Nij(1,2) * Iij{1,2}(2) * Iji{2,1}(2) + Nij(2,1) * Iij{2,1}(2) * Iji{2,1}(2)) / sigmaN^2;
%     GMG(2,1) = GMG(1,2);
%     GMG(2,2) = ((Nij(2,1) * Iij{2,1}(2)^2 + Nij(1,2) * Iji{2,1}(2)^2) / sigmaN^2) + ((Nij(2,1) + Nij(2,2)) / sigmag^2);
% 
%     GMB(1,1) = ((Nij(1,2) * Iij{1,2}(3)^2 + Nij(2,1) * Iji{1,2}(3)^2) / sigmaN^2) + ((Nij(1,1) + Nij(1,2)) / sigmag^2);
%     GMB(1,2) = - (Nij(1,2) * Iij{1,2}(3) * Iji{2,1}(3) + Nij(2,1) * Iij{2,1}(2) * Iji{2,1}(3)) / sigmaN^2;
%     GMB(2,1) = GMG(1,2);
%     GMB(2,2) = ((Nij(2,1) * Iij{2,1}(3)^2 + Nij(1,2) * Iji{2,1}(3)^2) / sigmaN^2) + ((Nij(2,1) + Nij(2,2)) / sigmag^2);
% 
%     bM (1,1) = (Nij(1,1) + Nij(1,2))/sigmag^2;
%     bM (2,1) = (Nij(2,2) + Nij(2,1))/sigmag^2;
%     
%     gR = GMR \ bM
%     gG = GMG \ bM
%     gB = GMB \ bM

%{
A=[1 2 3 4; 5 6 7 8; 9 10 11 12;13 14 15 16]
ii=ones(size(A))
idx=rot90(tril(rot90(ii)),-1);
A(~idx)=nan
%}

%{
aa = cell(3);
for i = 1:3
    for j = i:3
        aa{i,j} = [i,j]; 
    end
end

n = 5;
Ibarijvalue = zeros(n);
matSize = size(Ibarijvalue);
IuppeIdx = nonzeros(triu(reshape(1:numel(Ibarijvalue), size(Ibarijvalue))));
a = cell(length(IuppeIdx),1);
s = cell(n);

parfor i = 1:length(IuppeIdx)
    [ii,jj] = ind2sub(matSize, IuppeIdx(i));
    a{i} = [ii,jj];
end

s(triu(true(n))) = a;
%}

%{
function [gainpanorama, gainImages] = gainCompensation(input, warpedImages)
    
    % Initialze
    n = length(warpedImages);
    gainImages = cell(1,n);

    vistedpair = zeros(n^2,2);
    panoramasize = size(warpedImages{1});

    cnt = 1;
    for i = 1:n
        % Overlay the warpedImage onto the panorama.
        maski = imbinarize(rgb2gray(255 * warpedImages{i}));
        
        for j =  1:n
            exists = ismember([i j], vistedpair, 'rows');
            if (i~=j && ~ exists(1))

                Ibarij = zeros(panoramasize,'uint8');
                Ibarji = zeros(panoramasize,'uint8');

                vistedpair(cnt,:) = [j,i];
                
                maskj = imbinarize(rgb2gray(255 * warpedImages{j}));
                Nij   = maski & maskj;
                Nij = imfill(Nij, 'holes');

                Nijidx = repmat(Nij, 1, 1, 3);

                Iij = warpedImages{i};
                Iji = warpedImages{j};

                Ibarij(Nijidx) = Iij(Nijidx);
                Ibarji(Nijidx) = Iji(Nijidx);
                
                Ibarij_double = double(Ibarij);
                Ibarji_double = double(Ibarji);
                
                Ibarijvalue = sum(sum(Ibarij_double)) ./ sum(sum(Nij));
                Ibarjivalue = sum(sum(Ibarji_double)) ./ sum(sum(Nij)); 

                gi = Ibarjivalue ./ Ibarijvalue;
                gj = (Ibarijvalue ./ Ibarjivalue) .* gi;

                % Get gain compensated images
                if sum(sum(Nij)) > 0                     
                    if sum(gi(:,:,1:3) > 0.8) == 3 && sum(gj(:,:,1:3) > 0.8) == 3
                        if sum(gi(:,:,1:3) > 1) == 3
                            gi = ones(1,1,3);
                        end
                        if sum(gj(:,:,1:3) > 1) == 3
                            gj = ones(1,1,3);
                        end

                        gainImages{i} = uint8(double(warpedImages{i}) .* gi);
                        gainImages{j} = uint8(double(warpedImages{j}) .* gj);
                    end
                end
                cnt = cnt + 1;
            end
        end
    end

    % Apply no blending        
    iterations = size(gainImages{1},1:2);
    panorama = zeros(prod(iterations),3);
    
    tic
    for im = 1:length(gainImages)
        nextImage = gainImages{im};
        parfor i = 1:prod(iterations) 
            % IND2SUB converts from a "linear" index into individual
            % subscripts
            [ii,jj] = ind2sub(iterations, i);          
            if sum(nextImage(ii,jj,:)) ~= 0
                panorama(i,:) = nextImage(ii,jj,:);
            end
        end
    end
    gainpanorama = uint8(reshape(panorama,[iterations,3]));     
end
%}

%{
gi(1,1,1) = 0.81;
gi(1,1,2) = 0.81;
gi(1,1,3) = 0.81;

gj(1,1,1) = 0.81;
gj(1,1,2) = 0.81;
gj(1,1,3) = 0.81;

sum(gi(:,:,1:3) > 0.8 & gi(:,:,1:3) < 1)

sum(gj(:,:,1:3) > 0.8 & gj(:,:,1:3) < 1)
%}

%{
[img,cmap] = imread('lena.png');
img = img(:,:,1);
[ny,nx] = size(img);
% since there are an even number of columns and rows...
y1d = -ny/2+1/2:-1/2+ny/2;
x1d = -nx/2+1/2:-1/2+nx/2;
[X2D,Y2D] = meshgrid(x1d,y1d);
figure()
imagesc(x1d,y1d,img);
[theta2D,radius2D] = cart2pol(X2D,Y2D);
%}

%{
[theta,phi] = ndgrid(0:.1:2*pi, 0:.1:2*pi);
fx = 1./sin(theta).^2.*cos(phi).^4 + cos(theta).^2.*cos(phi).^2;
[X,Y,Z] = sph2cart(theta,phi,15);
surf(X,Y,Z);
%}

%{
warpedImageMask = cell(k, 1);
warpedImageNoBlend = zeros([panoramaView.ImageSize,3],'uint8');
% Overlay the warpedImage onto the panorama.
warpedImageMask{i} = imbinarize(rgb2gray(255 * warpedImage));

if (i~=1)
    warpedImageMaskCurrent = warpedImageMask{i} & warpedImageMask{i-1};
    warpedImageMaskCurrent = imsubtract(warpedImageMask{i}, warpedImageMaskCurrent);
else
    warpedImageMaskCurrent = warpedImageMask{i};
end
warpedImageNoBlend = warpedImageNoBlend + warpedImage .* uint8(repmat(warpedImageMaskCurrent, 1, 1, 3));
%}

%{
[sample, map, sample_alpha] = imread('peppers.png');
imshow(sample, 'alphadata', sample_alpha);

d=0;
%}


%{
colorArray = [0 0 0
    1 0 0
    0 1 0
    0 0 1
    1 1 1];

% colorArray  = [0 0 1
%                0 1 0
%                1 1 0
%                1 0 0];

colorArray = rand(15,3);

cmap = customColormap(colorArray, 1, 'linear');

[X,Y,Z] = peaks(25);
surf(X,Y,Z)
colormap(cmap);
colorbar;
xlabel('X'); ylabel('Y'); zlabel('Z');
%}

%{
fileName = 'peppers.png';
image = (imread(fileName));
f = 200;
k1 = 0.0;
k2 = 0.0;
k3 = 0.0;
imageCylindrical = image2cylindrical(image, f, k1, k2, k3);
imageSpherical = image2spherical(image, f, k1, k2, k3);

figure;
subplot(1,3,1); imshow(image); title('Input image')
subplot(1,3,2); imshow(imageCylindrical); title('Cylindrical projection')
subplot(1,3,3); imshow(imageSpherical); title('Spherical projection')
%}


%{
fileName = 'peppers.png';
image = imread(fileName);
f = 200;

[h, w, bypixs] = size(image);

K = [f,0,w/2; 0,f,h/2; 0,0,1];
x = 1:w;
y = 1:h;  

[X,Y] = meshgrid(x,y);

X = X';
Y = Y';

XYZ = [X(:) Y(:) ones(h*w,1)];
XYZnew =  (K \ XYZ')';

A = [sin(XYZnew(:,1)), XYZnew(:,2), cos(XYZnew(:,1))];
B = (K * A')';

% back from homog coords
B = B(:,1:2) ./ B(:,3);

% Make sure warp coords only within image bounds
B((B(:,1) < 0) | (B(:,1) >= w) | (B(:,2) < 0) | (B(:,2) >= h), :) = -1;
               
B = pagetranspose(reshape(B, w, h, 2));

% Image warp
Icy = imwarp(image, B, 'FillValues', [0 0 0]);

% Display warped image
% figure; 
% subplot(1,2,1); imshow(B(:,:,1))
% subplot(1,2,2); imshow(B(:,:,2))

figure;
subplot(1,2,1); imshow(image)
subplot(1,2,2); imshow(Icy)

dst = cv.remap(image, B(:,:,1), B(:,:,2), 'BorderType',  0);

figure; imshow(dst)
%}

%{
fileName = 'peppers.png';
image = imread(fileName);

f = 200;
k1 = 0.0;
k2 = 0.0;

% Get image size
[ydim, xdim, bypixs] = size(image);

% Initialize an array
out = uint8(zeros(ydim, xdim, 3));

% Get the center of image
xc = round(xdim/2);
yc = round(ydim/2);

% Create X and Y coordinates grid
[X,Y] = meshgrid(1:xdim, 1:ydim);

% Perform the cylindrical projection
theta = (X - xc)/f;
h = (Y - yc)/f;
xcap = sin(theta);
ycap = h;
zcap = cos(theta);
xn = xcap ./ zcap;
yn = ycap ./ zcap;
r = xn.^2 + yn.^2;

%   
xd = xn .* (1 + k1 * r + k2 * r.^2);
yd = yn .* (1 + k1 * r + k2 * r.^2);

% Convert to floor
ximg = floor(f * xd + xc);
yimg = floor(f * yd + yc);

% Find boundary of the cylindrical projection
[row,col] = find(ximg > 0 & ximg <= xdim & yimg > 0 & yimg <= ydim);

out(row,col,:) = image(row,col,:);

% for i=1:length(row)    
%    out(row(i),col(i),:) = image(row(i),col(i),:); 
% end
                               
figure;
subplot(1,2,1); imshow(image)
subplot(1,2,2); imshow(out)
%}


%{
fileName = 'peppers.png';
image = imread(fileName);

f = 200;
k1 = 0.;
k2 = 0.;

ydim=size(image, 1);
xdim=size(image, 2);

xc=xdim/2;
yc=ydim/2;

out = uint8(zeros(ydim, xdim, 3));

for y=1:ydim
    for x=1:xdim
        theta = (x - xc)/f;
        h = (y - yc)/f;
        xcap = sin(theta);
        ycap = h;
        zcap = cos(theta);
        xn = xcap / zcap;
        yn = ycap / zcap;
        r = xn^2 + yn^2;
        
        xd = xn * (1 + k1 * r + k2 * r^2);
        yd = yn * (1 + k1 * r + k2 * r^2);
        
        ximg = floor(f * xd + xc);
        yimg = floor(f * yd + yc);
        
        if (ximg > 0 && ximg <= xdim && yimg > 0 && yimg <= ydim)
            out(y, x, :) = [image(yimg, ximg, 1) image(yimg, ximg, 2) image(yimg, ximg, 3)];
        end
                               
    end
end

figure;
subplot(1,2,1); imshow(image)
subplot(1,2,2); imshow(out)
%}

%{
fileName = 'peppers.png';
I = imread(fileName);

imageinfo  = imfinfo(fileName);
[rows, cols, bypixs] = size(I);
xc = round(cols/2);
yc = round(rows/2);

f = 800;

x = 1:cols;
y = 1:rows;

[X,Y] = meshgrid(x,y);

dCy1 = f * tan(X-xc/f) + xc;
dCy2 = yc + (Y-yc ./ cos(X-xc/f));

% Icy = imwarp(I,dCyl);
Icy = cv.remap(image, dCyl, dCy2, 'BorderType',  0);

figure;
subplot(1,2,1); imshow(I)
subplot(1,2,2); imshow(Icy)

%}

%{
folderPath = 'D:\OneDrive\Team Work\Team CrackSTITCH\data\Generic';
folderName = 'crack exp4';
imageDir = fullfile(folderPath,folderName);

% imageDir = fullfile(toolboxdir('vision'),'visiondata','structureFromMotion');
images = imageDatastore(imageDir);
numImages = length(images.Files);

% Initialize variable to hold image sizes.
imageSizeAll = zeros(numImages,2);
imageinfo = imfinfo(char(images.Files(1)));

% Compute features for the first image.
I = rgb2gray(readimage(images,1));
% pointsPrev = detectSURFFeatures(I);
% [featuresPrev, pointsPrev] = extractFeatures(I,pointsPrev);

prevPoints   = detectSURFFeatures(I, 'NumOctaves', 8);

% Extract features. Using 'Upright' features improves matching, as long as
% the camera motion involves little or no in-plane rotation.
prevFeatures = extractFeatures(I, prevPoints, 'Upright', true);

focalLength    = [800, 800]; 
principalPoint = round(size(I,1:2)/2);
imageSize      = size(I,1:2);
intrinsics  = cameraIntrinsics(focalLength,principalPoint,imageSize);

% % Create an image view set object and add feature points of first image.
vSet = imageviewset;
vSet = addView(vSet,1,'Points',prevPoints);
% 
% % Compute features and matches for the rest of the images and add to the image view set.
% for i = 2:numel(images.Files)
%  I = rgb2gray(readimage(images,i));
%  
%  % Save image size.
%  imageSizeAll(i,:) = size(I);
%     
%  points = detectSURFFeatures(I);
%  [features, points] = extractFeatures(I,points);
%  
%  vSet = addView(vSet, i,'Features',features,'Points',points);
%  pairsIdx = matchFeatures(featuresPrev,features);
%  
%  vSet = addConnection(vSet,i-1,i,'Matches',pairsIdx);
%  featuresPrev = features;
% end
% 
% % Find point tracks across the sequence of images.
% pointTracks = findTracks(vSet);

% Add the Rest of the Views
tic
for i = 2:numImages
    I = rgb2gray(readimage(images,i));
 
     % Save image size.
     imageSizeAll(i,:) = size(I);
    
    % Detect, extract and match features.
    currPoints   = detectSURFFeatures(I, 'NumOctaves', 8);
    currFeatures = extractFeatures(I, currPoints, 'Upright', true);    
    indexPairs   = matchFeatures(prevFeatures, currFeatures, ...
        'MaxRatio', .7, 'Unique',  true);
    
    % Select matched points.
    matchedPoints1 = prevPoints(indexPairs(:, 1));
    matchedPoints2 = currPoints(indexPairs(:, 2));
    
    % Estimate the camera pose of current view relative to the previous view.
    % The pose is computed up to scale, meaning that the distance between
    % the cameras in the previous view and the current view is set to 1.
    % This will be corrected by the bundle adjustment.
    [relativeOrient, relativeLoc, inlierIdx] = helperEstimateRelativePose(...
        matchedPoints1, matchedPoints2, intrinsics);
    
    % Get the table containing the previous camera pose.
    prevPose = poses(vSet, i-1).AbsolutePose;
    relPose  = rigid3d(relativeOrient, relativeLoc);
        
    % Compute the current camera pose in the global coordinate system 
    % relative to the first view.
    currPose = rigid3d(relPose.T * prevPose.T);
    
    % Add the current view to the view set.
    vSet = addView(vSet, i, currPose, 'Points', currPoints);

    % Store the point matches between the previous and the current views.
    vSet = addConnection(vSet, i-1, i, relPose, 'Matches', indexPairs(inlierIdx,:));
    
    % Find point tracks across all views.
    tracks = findTracks(vSet);

    % Get the table containing camera poses for all views.
    camPoses = poses(vSet);

    % Triangulate initial locations for the 3-D world points.
    xyzPoints = triangulateMultiview(tracks, camPoses, intrinsics);
    
    % Refine the 3-D world points and camera poses.
    [xyzPoints, camPoses, reprojectionErrors] = bundleAdjustment(xyzPoints, ...
        tracks, camPoses, intrinsics, 'FixedViewId', 1, ...
        'PointsUndistorted', false);

    % Store the refined camera poses.
    vSet = updateView(vSet, camPoses);

    prevFeatures = currFeatures;
    prevPoints   = currPoints;  
end
toc

% Display camera poses.
camPoses = poses(vSet);

figure;
plotCamera(camPoses,'Color','g','Size',0.2);
hold on

% Exclude noisy 3-D points.
goodIdx = (reprojectionErrors < 5);
xyzPoints = xyzPoints(goodIdx, :);

% Display the 3-D points.
pcshow(xyzPoints, [1 0 0], 'VerticalAxis', 'y', 'VerticalAxisDir', 'down', ...
    'MarkerSize', 45);
grid on
hold off

%%
for i = 1:size(refinedPoses,1)
    
    tm = refinedPoses.AbsolutePose(i).T;
    tm(3,:) = [];
    tm(:,3) = [];
    tforms(i) = projective2d(tm);  
end

save tforms.mat tforms

%%
for i = 1:numel(tforms)           
    [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSizeAll(i,2)], [1 imageSizeAll(i,1)]);
end

maxImageSize = max(imageSizeAll);

% Find the minimum and maximum output limits. 
xMin = min([1; xlim(:)]);
xMax = max([maxImageSize(2); xlim(:)]);

yMin = min([1; ylim(:)]);
yMax = max([maxImageSize(1); ylim(:)]);

% Width and height of panorama.
width  = round(xMax - xMin);
height = round(yMax - yMin);

% Initialize the "empty" panorama.
panorama = zeros([height width 3], 'like', I);

%%
blender = vision.AlphaBlender('Operation', 'Binary mask', ...
    'MaskSource', 'Input port');  

% Create a 2-D spatial reference object defining the size of the panorama.
xLimits = [xMin xMax];
yLimits = [yMin yMax];
panoramaView = imref2d([height width], xLimits, yLimits);

% Create the panorama.
for i = 1:numImages
    
    I = readimage(images, i);   
   
    % Transform I into the panorama.
    warpedImage = imwarp(I, tforms(i), 'OutputView', panoramaView);
                  
    % Generate a binary mask.    
    mask = imwarp(true(size(I,1),size(I,2)), tforms(i), 'OutputView', panoramaView);
    
    % Overlay the warpedImage onto the panorama.
    panorama = step(blender, panorama, warpedImage, mask);
end

figure
imshow(warpedImage)

%}

%%
%{
% Refine the camera poses and points.
data = load('sfmGlobe');

[xyzRefinedPoints,refinedPoses] = ...
    bundleAdjustment(data.xyzPoints,data.pointTracks,data.cameraPoses,data.intrinsics);

% Display the refined camera poses and 3-D world points.
pcshow(data.xyzPoints, [1 0 0], 'VerticalAxis','y','VerticalAxisDir',...
    'down','MarkerSize',45);
hold on
pcshow(xyzRefinedPoints, [0 1 0], 'VerticalAxis','y','VerticalAxisDir',...
    'down','MarkerSize',45);
plotCamera(data.cameraPoses, 'Color','r', 'Size',0.1);
plotCamera(refinedPoses,'Color','g','Size',0.1);
hold off
grid on
%}
    
%%
%{
data = load('motionOnlyBA.mat');

impts = data.imagePoints

% Refine the absolute camera poses.  
% refinedPose = bundleAdjustmentStructure(data.xyzPoints,data.pointTracks,data.cameraPoses,data.intrinsics);

refinedPose = bundleAdjustmentMotion(data.xyzPoints,data.imagePoints,data.absPose,data.intrinsics);

% Display the 3-D world points.
pcshow(data.xyzPoints,'VerticalAxis','y','VerticalAxisDir','down','MarkerSize',45);
hold on

% Plot the absolute camera poses before and after refinement.
plotCamera('AbsolutePose',data.absPose,'Color','r','Size',2);
plotCamera('AbsolutePose',refinedPose,'Color','m','Size',2);

%}
%%
% n = 3;
% M = n;
% N = n;
% numMatches = zeros(n);
% parfor i = 1:M*N
%     % IND2SUB converts from a "linear" index into individual
%     % subscripts
%     [ii,jj] = ind2sub([M,N], i);
%     jj=jj+ii-1
%     
%     if (jj>ii)      
%         numMatches(i) = randi(10);
%     end
% end


% pairs = zeros(0,2);
% n = 3;
% for i = 1:n  % or reallyl just to n-1?
%     for j = i+1:n
%         pairs(end+1,:) = [i, j];
%     end   
% end
% npairs = size(pairs,1);
% temp = zeros(npairs,1);
% parfor i = 1:npairs
%     ii = pairs(i,1);
%     jj = pairs(i,2);
%     matches = rand(1);
%     nf = size(matches, 2);
%     temp(i) = rand(1);
% end
% numMatches = zeros(n,n);
% for i=1:npairs
%     ii = pairs(i,1);
%     jj = pairs(i,2);
%     numMatches(ii,jj) = temp(i);
% end

% adjacencyMatrix =[0 1 1 0 1 0 1 0 1 0
%     1 0 1 0 1 0 1 0 1 0
%     1 1 0 0 1 0 1 0 1 0
%     0 0 0 0 0 0 0 1 0 0
%     1 1 1 0 0 0 1 0 1 0
%     0 0 0 0 0 0 0 0 0 1
%     1 1 1 0 1 0 0 0 1 0
%     0 0 0 1 0 0 0 0 0 1
%     1 1 1 0 1 0 1 0 0 0
%     0 0 0 0 0 1 0 1 0 0]
% G = graph(adjacencyMatrix);
% plot(G);  %view the graph
% bins = conncomp(G);
% binnodes = accumarray(bins', 1:numel(bins), [], @(v) {sort(v')});
% fprintf('number of Regions = %d\n\n', numel(binnodes));
% for binidx = 1:numel(binnodes)
%     fprintf('All these nodes are connected:%s\n', sprintf(' %d', binnodes{binidx}));
% end
% fprintf('\n');



% M = 3;
% N = 4;
% P = 5;
% Q = 6;
% % Approach 1: parallelize only the outer loop
% out1 = zeros(M, N, P, Q);
% parfor ii = 1:M
%     for jj = 1:N
%         for kk = 1:P
%             for ll = 1:Q
%                 out1(ii,jj,kk,ll) = ii * 1000 + jj * 100 + kk * 10 + ll;
%             end
%         end
%     end
% end
% 
% % Approach 2: convert to linear indexing
% out2 = zeros(M, N, P, Q);
% parfor idx = 1:(M*N*P*Q)
%     % IND2SUB converts from a "linear" index into individual
%     % subscripts
%     [ii,jj,kk,ll] = ind2sub([M,N,P,Q], idx);
%     out2(idx) = ii * 1000 + jj * 100 + kk * 10 + ll;
% end
% % Check
% assert(isequal(out1, out2))


% n = 3;
% M = n;
% N = n;
% out3 = zeros(n);
% 
% parfor idx = 1:(n * n)
%     % IND2SUB converts from a "linear" index into individual
%     % subscripts
% %     [ii,jj] = ind2sub([M,N], idx);
%     out3(idx) = idx;
% end

% n=3;
% allMatches = cell(n);
% numMatches = zeros(n);
% M = n;
% N = n;

% f=@(a,b,c)([a,b,c]);

% parfor i = 1:M %(M*N)
%     % IND2SUB converts from a "linear" index into individual
%     % subscripts
%     [ii,jj] = ind2sub([M], i);
% 
% %     if (ii~=jj)
%         allMatches{i} = [ii,jj];
%         numMatches(i) = jj;
% %     end
% end


% [1075,28,179;
%    0,792,284;
%  171,274,1204]

% [0,26,178;
%  0,0,284;
%  0,0,0]

% n = 3;
% for i = 1:n
%     for j = i+1:n
%         % Feature matching
%         matches = getMatches(input, allDescriptors{i}, allDescriptors{j});
%         nf = size(matches, 2);
%         numMatches(i,j) = nf;
%     end   
% end

% n = 3;
% M = n;
% N = n;
% parfor i = 1:M*N
%     % IND2SUB converts from a "linear" index into individual
%     % subscripts
%     [ii,jj] = ind2sub([M,N], i);
%     
%     if (ii~=jj)
%         matches = getMatches(input, allDescriptors{ii}, allDescriptors{jj});
%         nf = size(matches, 2);
%         numMatches(i) = nf;
%     end
% end


% pairs = zeros(0,2);
% n = 3;
% for i = 1:n  % or reallyl just to n-1?
%     for j = i+1:n
%         pairs(end+1,:) = [i, j];
%     end   
% end
% 
% % Then your parfor loop can just go through the preset pairs:
% npairs = size(pairs,1);
% parfor i = 1:npairs
%     ii = pairs(i,1);
%     jj = pairs(i,2);
% %         matches = getMatches(input, allDescriptors{ii}, allDescriptors{jj});
% %         nf = size(matches, 2);
%         numMatches(i) = rand(1);
% end

% A=[1 2 3 4; 5 6 7 8; 9 10 11 12;13 14 15 16]
% ii=ones(size(A))
% idx=rot90(tril(rot90(ii)),-1);
% A(~idx)=nan
