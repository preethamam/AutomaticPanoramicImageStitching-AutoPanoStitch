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
    Iij = cell(n);
    Iji = cell(n);
    Nij = zeros(n);

    matSize = size(Iij);
    IuppeIdx = nonzeros(triu(reshape(1:numel(Iij), size(Iij))));
    Ibarijvalue = cell(length(IuppeIdx),1);
    Ibarjivalue = cell(length(IuppeIdx),1);
    Nijvalue    = zeros(length(IuppeIdx),1);

    parfor i = 1:length(IuppeIdx)
        
        [ii,jj] = ind2sub(matSize, IuppeIdx(i));

        if ii ~= jj
            Ibarij = zeros(panoramasize,'uint8');
            Ibarji = zeros(panoramasize,'uint8');
        
            % Overlay the warpedImage onto the panorama.
            maski = imbinarize(rgb2gray(255 * warpedImages{ii}));  
            maskj = imbinarize(rgb2gray(255 * warpedImages{jj}));
    
            Nij_im  = maski & maskj;
            Nij_im  = imfill(Nij_im, 'holes');
    
            Nijidx = repmat(Nij_im, 1, 1, 3);
    
            Imij = warpedImages{ii};
            Imji = warpedImages{jj};
    
            Ibarij(Nijidx) = Imij(Nijidx);
            Ibarji(Nijidx) = Imji(Nijidx);
            
            Ibarij_double = double(Ibarij);
            Ibarji_double = double(Ibarji);
    
            Nij_sum = sum(sum(Nij_im));
            
            Ibarijvalue{i} = reshape(sum(sum(Ibarij_double)) ./ Nij_sum, 1, 3);
            Ibarjivalue{i} = reshape(sum(sum(Ibarji_double)) ./ Nij_sum, 1, 3); 

            Nijvalue(i)    = Nij_sum;
        end

    end

    % --------------------------------------------------------------------------------------------------------------
    % Form matrices
    Iij(triu(true(n))) = Ibarijvalue;
    Iji(triu(true(n))) = Ibarjivalue;
    Nij(triu(true(n))) = Nijvalue;

    Iijc = Iij;
    Ijic = Iji;

    iii = ones(size(Iij));
    idx = find(tril(iii,-1));

    Iijp = Iij';
    Ijip = Iji';

    Iijc(idx) = Iijp(idx);
    Ijic(idx) = Ijip(idx);

    Iij = Iijc;
    Iji = Ijic;
    Nij = Nij + tril(Nij',-1);

    Iij = cellfun(@(M) subsasgn(M, substruct('()', {isnan(M)}), 0), Iij, 'uniform', 0); %#ok<*SUBSASGN> 
    Iji = cellfun(@(M) subsasgn(M, substruct('()', {isnan(M)}), 0), Iji, 'uniform', 0);

    % --------------------------------------------------------------------------------------------------------------
    % Gain values matrix
    gainmatR = zeros(n);
    gainmatG = zeros(n);
    gainmatB = zeros(n);

    gainR    = zeros(length(IuppeIdx),1);
    gainG    = zeros(length(IuppeIdx),1);
    gainB    = zeros(length(IuppeIdx),1);

    parfor i = 1:length(IuppeIdx)

        [ii,jj] = ind2sub(matSize, IuppeIdx(i));

        if ii ~= jj
            gainR(i) = -(Nij(ii,jj) * Iij{ii,jj}(1) * Iji{jj,ii}(1) + Nij(jj,ii) * Iij{jj,ii}(1) * Iji{jj,ii}(1)) / sigmaN^2;
            gainG(i) = -(Nij(ii,jj) * Iij{ii,jj}(2) * Iji{jj,ii}(2) + Nij(jj,ii) * Iij{jj,ii}(2) * Iji{jj,ii}(2)) / sigmaN^2;
            gainB(i) = -(Nij(ii,jj) * Iij{ii,jj}(3) * Iji{jj,ii}(3) + Nij(jj,ii) * Iij{jj,ii}(3) * Iji{jj,ii}(3)) / sigmaN^2;
        else

            gainRval = 0;
            gainGval = 0;
            gainBval = 0;

            for iii = 1:n
                if iii ~= ii
                    gainRval = gainRval + (((Nij(ii,iii) * Iij{ii,iii}(1)^2 + Nij(iii,ii) * Iji{ii,iii}(1)^2) / sigmaN^2) + (Nij(ii,iii) / sigmag^2));
                    gainGval = gainGval + (((Nij(ii,iii) * Iij{ii,iii}(2)^2 + Nij(iii,ii) * Iji{ii,iii}(2)^2) / sigmaN^2) + (Nij(ii,iii) / sigmag^2));
                    gainBval = gainBval + (((Nij(ii,iii) * Iij{ii,iii}(3)^2 + Nij(iii,ii) * Iji{ii,iii}(3)^2) / sigmaN^2) + (Nij(ii,iii) / sigmag^2));
                end
            end
            gainR(i) = gainRval;
            gainG(i) = gainGval;
            gainB(i) = gainBval;
        end
    end

    % Make matrices
    gainmatR(triu(true(n))) = gainR;
    gainmatG(triu(true(n))) = gainG;
    gainmatB(triu(true(n))) = gainB;

    gainmatR = gainmatR + tril(gainmatR',-1);
    gainmatG = gainmatG + tril(gainmatG',-1);
    gainmatB = gainmatB + tril(gainmatB',-1);

    % --------------------------------------------------------------------------------------------------------------
    % AX = b --> X = A \ b
    b = zeros(n,1);
    parfor j = 1:n
        for i = 1:n
            b(j) = b(j) + (Nij(j,i) / sigmag^2); %#ok<*PFBNS> 
        end
    end

    gR = gainmatR \ b;
    gG = gainmatG \ b;
    gB = gainmatB \ b;

    gainRGB = [gR, gG, gB];

    % --------------------------------------------------------------------------------------------------------------
    % Compensate gains for images
    parfor i = 1:n
        gain = reshape([gR(i), gG(i), gB(i)], [1,1,3]); 
        gainImages{i} = uint8(double(warpedImages{i}) .* gain);
    end

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