function [keypoints, allDescriptors, images, imageinfo, imageFocals, n] = featureMatching(input, imgSetVector, myImg)

    % Number of images in the folder
    n = max(cat(1,imgSetVector(myImg).Count));

    % Initialize the cell arrays
    keypoints = cell(1,n);
    allDescriptors = cell(1,n);
    images = cell(1,n);
    imageinfo = cell(1,n);
    imageFocals = zeros(1,n);
    
    parfor i = 1:n
        imageFile = imgSetVector(myImg).ImageLocation{i};
        image = imread(imageFile);
        
        % Get size of the image
        [imRows, imCols, imChannel] = size(image);
        
        % Camera intrinsics
        K = [input.fx, 0, imCols/2; 0, input.fy, imRows/2; 0, 0, 1];

        if strcmp(input.warpType,'spherical')
            image = image2spherical_v1(image, K, input.DC);
        elseif strcmp(input.warpType,'cylindrical')
            image = image2cylindrical_v1(image, K, input.DC);
        end

        % Replicate the third channel
        if imChannel == 1
            image = repmat(image, 1, 1, 3);
        end
        
        % Image resize
        if input.resizeImage
            if imRows > 480 && imCols > 640
                image = imresize(image,[480, 640]);
            elseif imRows > 480 && imCols < 640
                image = imresize(image,[480, imCols]);
            elseif imRows < 480 && imCols > 640
                image = imresize(image,[imRows, 640]);
            end
        end
        images{i} = image;

        % Get features and valid points
        [descriptors, points] = getFeaturePoints(input, image);
        
        % Concatenate the descriptors and points
        keypoints{i} = points';
        allDescriptors{i} = descriptors; 
    end
end
   
% Get the feature points
function [features, validPts] = getFeaturePoints(input, ImageOriginal)    
if size(ImageOriginal,3) > 1
    grayImage = rgb2gray(ImageOriginal);
else
    grayImage = ImageOriginal;
end
    
switch input.detector
    case 'HARRIS'
        points = detectHarrisFeatures(grayImage);

    case 'SIFT'        
        points = detectSIFTFeatures(grayImage,'NumLayersInOctave',input.NumLayersInOctave, ...
                                    ContrastThreshold=input.ContrastThreshold, ...
                                    EdgeThreshold=input.EdgeThreshold);

    case 'FAST'
        points   = detectFASTFeatures(grayImage);

    case 'SURF'
        points = detectSURFFeatures(grayImage, 'NumOctaves', 8);

    case 'BRISK'
        points  = detectBRISKFeatures(grayImage);

    case 'ORB'
        points   = detectORBFeatures(grayImage);

    case 'KAZE'
        points   = detectKAZEFeatures(grayImage);

    otherwise
        error('Need a valid input!')
end
    [features, validPts] = extractFeatures(grayImage, points);
    validPts = double(validPts.Location);
end

