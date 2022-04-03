function displayPanorama(finalPanoramaTforms, input, images, imageFocals, concomps, myImg, datasetName)

%%***********************************************************************%
%*                   Automatic panorama stitching                       *%
%*                        Display panorama                              *%
%*                                                                      *%
%* Code author: Preetham Manjunatha                                     *%
%* Github link: https://github.com/preethamam                           *%
%* Date: 01/27/2022                                                     *%
%************************************************************************%

for ii = 1:length(finalPanoramaTforms)
        % Create and display panorama
        [panorama, gainpanorama, gainImages, gainRGB, xCorrect, yCorrect] = ...
            renderPanorama(input, images, imageFocals, finalPanoramaTforms{ii}, concomps, ii);
        
        % Final panorama cropper
        croppedImage = panoramaCropper(input, panorama);   
        
        % Plot the bounding boxes
        figure(2);
        imshow(croppedImage)
        ax = gcf;
        exportgraphics(ax,'pano_crop.jpg')

        % Generate a random (bright) color
        lineRGB = mod(rand(size(xCorrect,2),3),100);
        max_channel = max(lineRGB(:,1), max(lineRGB(:,2), lineRGB(:,3)));
        lines_brightRGB = lineRGB ./ max_channel;
        
        % Plot the bounding boxes
        figure(3);
        imshow(panorama)
        if strcmp(input.warpType,'planar')
            hold on
            for i = 1:size(xCorrect,2)
                plot(xCorrect(:,i), yCorrect(:,i), 'Color', lines_brightRGB(i,:), 'LineWidth', 1)
                drawnow;
            end
            hold off
        end
               
        % Export the graphics
%         ax = gca;     

%         if ismac
%             % Code to run on Mac platform
%             exportgraphics(ax, fullfile('../../../Results/MATLAB Stitch/', [input.warpType '_' ...
%             num2str(myImg) '_' num2str(ii) '_' char(datasetName{myImg}) '.png']),'Resolution',300)
%         elseif isunix
%             % Code to run on Linux platform
%             exportgraphics(ax, fullfile(['../../../Results/MATLAB Stitch/', input.warpType '_'  ...
%             num2str(myImg) '_' num2str(ii) '_' char(datasetName{myImg}) '.png']),'Resolution',300)
%         elseif ispc
%             % Code to run on Windows platform
%             exportgraphics(ax, fullfile(['..\..\..\Results\MATLAB Stitch\', input.warpType '_'  ...
%             num2str(myImg) '_' num2str(ii) '_' char(datasetName{myImg}) '.png']),'Resolution',300)
%         else
%             disp('Platform not supported')
%         end

%         ax = gcf;
%         exportgraphics(ax,'pano_full.jpg')

        if ismac
            % Code to run on Mac platform
            imwrite(panorama, fullfile('../../../Results/MATLAB Stitch/', [input.warpType '_' ...
            num2str(myImg) '_' num2str(ii) '_' char(datasetName{myImg}) '.png']))
        elseif isunix
            % Code to run on Linux platform
            imwrite(panorama, fullfile(['../../../Results/MATLAB Stitch/', input.warpType '_'  ...
            num2str(myImg) '_' num2str(ii) '_' char(datasetName{myImg}) '.png']))
        elseif ispc
            % Code to run on Windows platform
            imwrite(panorama, fullfile(['..\..\..\Results\MATLAB Stitch\', input.warpType '_'  ...
            num2str(myImg) '_' num2str(ii) '_' char(datasetName{myImg}) '.png']))
        else
            disp('Platform not supported')
        end
end
end