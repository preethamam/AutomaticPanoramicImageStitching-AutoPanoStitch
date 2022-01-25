function displayPanorama(finalPanoramaTforms, input, images, imageFocals, concomps, myImg, datasetName)

for ii = 1:length(finalPanoramaTforms)
        % Create and display panorama
        [panorama, gainpanorama, gainImages, gainRGB, xCorrect, yCorrect] = ...
            renderPanorama(input, images, imageFocals, finalPanoramaTforms{ii}, concomps, ii);
        
        % Final panorama cropper
%         croppedImage = panoramaCropper(input, panorama);   

        % Generate a random (bright) color
        lineRGB = mod(rand(size(xCorrect,2),3),100);
        max_channel = max(lineRGB(:,1), max(lineRGB(:,2), lineRGB(:,3)));
        lines_brightRGB = lineRGB ./ max_channel;
        
        % Plot the bounding boxes
        figure(2);
        imshow(panorama)
        hold on
        for i = 1:size(xCorrect,2)
            plot(xCorrect(:,i), yCorrect(:,i), 'Color', lines_brightRGB(i,:), 'LineWidth', 1)
            drawnow;
        end
        hold off
               
        % Export the graphics
        ax = gca;      

        if ismac
            % Code to run on Mac platform
        elseif isunix
            % Code to run on Linux platform
            exportgraphics(ax, ['../../../Results/MATLAB Stitch/stitch_kevin_preetham_'  ...
            num2str(myImg) '_' num2str(ii) '_' char(datasetName{myImg}) '.png'],'Resolution',300)
        elseif ispc
            % Code to run on Windows platform
            exportgraphics(ax, ['..\..\..\Results\MATLAB Stitch\stitch_kevin_preetham_'  ...
            num2str(myImg) '_' num2str(ii) '_' char(datasetName{myImg}) '.png'],'Resolution',300)
        else
            disp('Platform not supported')
        end

        % Cropped view
%         figure(3);
%         imshow(croppedImage)
% 
%         figure(4);
%         imshow(gainpanorama)

%         input.canvas_color = 'white';
%         croppedImage2 = panoramaCropper(input, croppedImage);
% 
%         figure(4);
%         imshow(croppedImage2)
end
end