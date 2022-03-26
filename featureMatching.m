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
        
%         imageinfo{i}  = imfinfo(imageFile);
%         imageFocals(i) = imageinfo{i}.DigitalCamera.FocalLength;
        
%         image = image2cylindrical(image, 500, 0, 0, 0);
%           image = image2spherical(image, 800, 0, 0, 0);
        
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
        [locations, features]  = vl_sift(single(grayImage), 'Octaves', 8);
        features = features';
        if ~isempty(locations)                
            validPts = locations(1:2,:)';
        else
            validPts = [];                
        end
        
%         points = detectSIFTFeatures(grayImage,'NumLayersInOctave',3);

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

if ~(strcmp(input.detector, 'SIFT'))
    [features, validPts] = extractFeatures(grayImage, points);
    validPts = double(validPts.Location);
end

end

