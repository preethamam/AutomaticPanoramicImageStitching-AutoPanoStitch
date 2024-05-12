function [gainpanorama, gainImages, gainRGB] = gainCompensation(input, warpedImages)

%%***********************************************************************%
%*                   Automatic panorama stitching                       *%
%*                        Gain compensation                             *%
%*                                                                      *%
%* Code author: Preetham Manjunatha                                     *%
%* Github link: https://github.com/preethamam                           *%
%* Date: 01/27/2022                                                     *%
%************************************************************************%

    % Initialze
    n = length(warpedImages);
    gainImages = cell(1,n);
    sigmaN = input.sigmaN;
    sigmag = input.sigmag;

    panoramasize = size(warpedImages{1});
    Amat = cell(n);
    Bvec = zeros(n,1);   
    IuppeIdx = nonzeros(triu(reshape(1:numel(Amat), size(Amat))));
    Amat_temp = cell(1,length(IuppeIdx));    
    matSize = size(Amat);     
    
    % Get the Ibarijs and Nijs
    for i = 1:length(IuppeIdx)        
        [ii,jj] = ind2sub(matSize, IuppeIdx(i)); 
        if ii == jj
            diag_val_1 = 0;
            diag_val_2 = 0;
            diag_val = 0;
           
            Z = 1:n;
            Z(Z==ii) = [];
            for d = Z
                [Ibarij, Ibarji, Nij] = getIbarNij(panoramasize, warpedImages, ii, d);
                diag_val_1 = diag_val_1 + ( (Nij + Nij) .*  Ibarij.^2 );
                diag_val_2 = diag_val_2 + Nij;
            end

            diag_val = diag_val_1 + (sigmaN^2/sigmag^2) * diag_val_2;

            Amat_temp{i} = diag_val;
            Bvec(i) = (sigmaN^2/sigmag^2) * diag_val_2;
        end       

        if ii ~= jj
            [Ibarij,Ibarji,Nij] = getIbarNij(panoramasize, warpedImages, ii, jj);
            Amat_temp{i} = -(Nij+Nij) .* (Ibarij .* Ibarji);
        end

    end

    % --------------------------------------------------------------------------------------------------------------
    % Form matrices
    updiagIdx = nonzeros(triu(reshape(1:numel(Amat), size(Amat)),1));
    lowdiagIdx = nonzeros(tril(reshape(1:numel(Amat), size(Amat)),-1));
    
    % Populate A matrix
    Amat(IuppeIdx) = Amat_temp;
    Amat(lowdiagIdx) = Amat(updiagIdx);
    Bvec = nonzeros(Bvec);    
    Amat = cell2mat(Amat);    

    % A matrix gain values
    gainmatR = Amat(:,1:3:size(Amat,2));
    gainmatG = Amat(:,2:3:size(Amat,2));
    gainmatB = Amat(:,3:3:size(Amat,2));

    % --------------------------------------------------------------------------------------------------------------
    % AX = b --> X = A \ b
    gR = gainmatR \ Bvec;
    gG = gainmatG \ Bvec;
    gB = gainmatB \ Bvec;
    
    % Concatenate RGB gains
    gainRGB = [gR, gG, gB]; 

    % Clip gain values to 1
    gainRGB = min(gainRGB,1);

    % --------------------------------------------------------------------------------------------------------------
    % Compensate gains for images        
    gains = num2cell(gainRGB,2);
    gains = cellfun(@(x) reshape(x, [1,1,3]), gains, 'UniformOutput',false);
    gainImages = cellfun(@(x,y) uint8(double(x).*y), warpedImages, gains,'UniformOutput',false);

    % --------------------------------------------------------------------------------------------------------------
    % Construct gain comepensated panorama      
    iterations = size(gainImages{1},1:2);
    panorama = zeros(prod(iterations),3);

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

%--------------------------------------------------------------------------------------------------------
% Auxillary functions
%--------------------------------------------------------------------------------------------------------

%--------------------------------------------------------------------------------------------------------
% [Ibarij,Ibarji,Nij] = getIbarNij(panoramasize, warpedImages, ii, jj)
%--------------------------------------------------------------------------------------------------------
%
% Returns the Ibars and Nij values.
function [Ibarij,Ibarji,Nij] = getIbarNij(panoramasize, warpedImages, ii, jj)
    Ibarij = zeros(panoramasize,'uint8');
    Ibarji = zeros(panoramasize,'uint8');
    
    % Overlay the warpedImage onto the panorama.
    maski = imbinarize(rgb2gray(255 * warpedImages{ii}));  
    maskj = imbinarize(rgb2gray(255 * warpedImages{jj}));
    
    % Find the overlap mask
    Nij_im  = maski & maskj;
    Nij_im  = imfill(Nij_im, 'holes');    
    Nijidx = repmat(Nij_im, 1, 1, 3);
    
    % Warped images
    Imij = warpedImages{ii};
    Imji = warpedImages{jj};
    
    % Get the overlapping region RGB values for two images
    Ibarij(Nijidx) = Imij(Nijidx);
    Ibarji(Nijidx) = Imji(Nijidx);
    
    % Convert to double
    Ibarij_double = double(Ibarij);
    Ibarji_double = double(Ibarji);
    
    % Nij
    Nij = sum(sum(Nij_im));
    
    % Ibar ijs
    Ibarij = reshape(sum(sum(Ibarij_double)) ./ Nij, 1, 3);
    Ibarji = reshape(sum(sum(Ibarji_double)) ./ Nij, 1, 3); 
    
    % Replace NaNs by zeros
    Ibarij(isnan(Ibarij)) = 0;
    Ibarji(isnan(Ibarji)) = 0;
end