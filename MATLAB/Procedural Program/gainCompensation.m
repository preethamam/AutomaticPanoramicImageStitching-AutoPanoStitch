function [gainpanorama, gainImages, gainRGB] = gainCompensation(input, warpedImages, imageNeighbors)

    %%***********************************************************************%
    %*                   Automatic panorama stitching                       *%
    %*                        Gain compensation                             *%
    %*                                                                      *%
    %* Code author: Preetham Manjunatha                                     *%
    %* Github link: https://github.com/preethamam                           *%
    %* Date: 05/12/2024                                                     *%
    %************************************************************************%

    % Initialze
    n = length(warpedImages);
    sigmaN = input.sigmaN;
    sigmag = input.sigmag;
    gainDerivation = input.gainDerivation;
    parforSummation = input.parforSummation;
    Amat = cell(n);
    Bvec = zeros(n,n);
    Bvec_temp = zeros(n,1);
    IuppeIdx = nonzeros(triu(reshape(1:numel(Amat), size(Amat))));
    Amat_temp = cell(1,length(IuppeIdx));    
    matSize = size(Amat);

    % Get the Ibarijs and Nijs
    parfor i = 1:length(IuppeIdx)

        % Index to subscripts
        [ii,jj] = ind2sub(matSize, IuppeIdx(i));                 
        
        % Diagonal entries
        if ii == jj                      
            switch parforSummation
                case true
                    [diag_val_1, diag_val_2] = gainDiagonalSum(warpedImages, ii, imageNeighbors{ii}, gainDerivation);
                case false
                    diag_val_1 = 0;
                    diag_val_2 = 0;
                    for d = imageNeighbors{ii}
                        [Ibarij, Ibarji, Nij] = getIbarNij(warpedImages{ii}, warpedImages{d});
                        switch gainDerivation
                            case 1
                                diag_val_1 = diag_val_1 + ((Nij + Nij) .*  (Ibarij .* Ibarij));
                                diag_val_2 = diag_val_2 + Nij;
                            case 2
                                diag_val_1 = diag_val_1 + ((Nij * Ibarij.^2 + Nij * Ibarij.^2) / sigmaN^2);
                                diag_val_2 = diag_val_2 + (Nij / sigmag^2);
                        end
                    end
            end

            switch gainDerivation
                case 1
                    Amat_temp{i} = diag_val_1 + (sigmaN^2/sigmag^2) * diag_val_2;
                    Bvec_temp(i) = (sigmaN^2/sigmag^2) * diag_val_2;
                case 2
                    Amat_temp{i} = diag_val_1 + diag_val_2;
                    Bvec_temp(i) = diag_val_2;
            end
        end       
        
        % Off-diagonal entries
        if ii ~= jj
            [Ibarij,Ibarji,Nij] = getIbarNij(warpedImages{ii}, warpedImages{jj});
            switch gainDerivation
                case 1
                    Amat_temp{i} = - (Nij+Nij) .* (Ibarij .* Ibarji);
                case 2
                    Amat_temp{i} = - ((Nij+Nij) * (Ibarij .* Ibarji)) /sigmaN^2 ;
            end
        end

    end

    % --------------------------------------------------------------------------------------------------------------    
    % Form A matrix and B vector
    Amat(IuppeIdx) = Amat_temp;
    Amat_t = Amat';
    Amat(tril(true(size(Amat)))) = Amat_t(tril(true(size(Amat))));

    Bvec(IuppeIdx) = Bvec_temp; 
    Bvec = diag(Bvec);
    Amat = cell2mat(Amat);    

    % A matrix gain values
    gainmatR = Amat(:,1:3:size(Amat,2)) + eps;
    gainmatG = Amat(:,2:3:size(Amat,2)) + eps;
    gainmatB = Amat(:,3:3:size(Amat,2)) + eps;

    % --------------------------------------------------------------------------------------------------------------
    % AX = b --> X = A \ b
    gR = gainmatR \ Bvec;
    gG = gainmatG \ Bvec;
    gB = gainmatB \ Bvec;
    
    % Concatenate RGB gains
    gainRGB = [gR, gG, gB];

    % Clip negatives
    gainRGB(gainRGB < 0) = 1;

    % --------------------------------------------------------------------------------------------------------------
    % Compensate gains for images        
    gains = num2cell(gainRGB,2);
    gains = cellfun(@(x) reshape(x, [1,1,3]), gains, 'UniformOutput',false);
    gainImages = cellfun(@(x,y) uint8(double(x).*y), warpedImages, gains,'UniformOutput',false);
    gainImagesMask = cellfun(@(x) repmat(imfill(imbinarize(rgb2gray(255 * x)), 'holes'), 1, 1, size(warpedImages{1},3)), ...
                            gainImages, 'UniformOutput',false);
   
    % --------------------------------------------------------------------------------------------------------------
    % Construct gain comepensated panorama    
    gainpanorama = zeros(size(warpedImages{1}), 'uint8');
    for i = 1:n
        gainpanorama(gainImagesMask{i}) = gainImages{i}(gainImagesMask{i});
    end
end

%--------------------------------------------------------------------------------------------------------
% Auxillary functions
%--------------------------------------------------------------------------------------------------------

%--------------------------------------------------------------------------------------------------------
% [Ibarij,Ibarji,Nij] = getIbarNij(Imij, Imji)
%--------------------------------------------------------------------------------------------------------
%
function [Ibarij,Ibarji,Nij] = getIbarNij(Imij, Imji)

    % Overlay the warpedImage onto the panorama.
    maski = imbinarize(rgb2gray(255 * Imij));  
    maskj = imbinarize(rgb2gray(255 * Imji));

    % Find the overlap mask
    Nij_im  = maski & maskj;
    
    % Get tje Ibarij and Nijs
    if sum(sum(Nij_im)) > 1
        Nij_im  = imfill(Nij_im, 'holes');    
        Nijidx  = repmat(Nij_im, 1, 1, size(Imij,3));
    
        % Get the overlapping region RGB values for two images
        Ibarij = double(reshape(Imij(Nijidx),1,[], size(Imij,3)));
        Ibarji = double(reshape(Imji(Nijidx),1,[], size(Imij,3)));
    
        % Nij
        Nij = sum(sum(Nij_im));
    
        % Ibar ijs
        Ibarij = reshape(sum(Ibarij) ./ Nij, 1, size(Imij,3));
        Ibarji = reshape(sum(Ibarji) ./ Nij, 1, size(Imij,3));
    else
        % Nij
        Nij = 0;

        % Ibar ijs
        Ibarij = zeros(1, size(Imij,3));
        Ibarji = zeros(1, size(Imij,3));
    end
end

%--------------------------------------------------------------------------------------------------------
% [diag_val_1, diag_val_2] = gainDiagonalSum(warpedImages, ii, Z, gainDerivation)
%--------------------------------------------------------------------------------------------------------
%
function [diag_val_1, diag_val_2] = gainDiagonalSum(warpedImages, ii, Z, gainDerivation) %#ok<DEFNU>
    diag_val_1 = 0;
    diag_val_2 = 0;
    parfor d = 1:length(Z)
        dval = Z(d);
        [Ibarij, Ibarji, Nij] = getIbarNij(warpedImages{ii}, warpedImages{dval});
        switch gainDerivation
            case 1
                diag_val_1 = diag_val_1 + ( (Nij + Nij) .*  (Ibarij .* Ibarij) );
                diag_val_2 = diag_val_2 + Nij;
            case 2
                diag_val_1 = diag_val_1 + ((Nij * Ibarij.^2 + Nij * Ibarij.^2) / sigmaN^2);
                diag_val_2 = diag_val_2 + (Nij / sigmag^2);
        end
    end
end