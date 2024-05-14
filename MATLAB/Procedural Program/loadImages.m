function imageFiles = loadImages(imgSetVector, myImg)
    
    %%***********************************************************************%
    %*                   Automatic panorama stitching                       *%
    %*                          Images loader                               *%
    %*                                                                      *%
    %* Code author: Preetham Manjunatha                                     *%
    %* Github link: https://github.com/preethamam                           *%
    %* Date: 05/14/2024                                                     *%
    %************************************************************************%

    % Read images    
    imgFolder = fileparts(imgSetVector(myImg).ImageLocation(1));
    imds = imageDatastore(imgFolder);
    imageFiles = readall(imds);    
end