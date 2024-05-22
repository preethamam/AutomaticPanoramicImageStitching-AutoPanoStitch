function matches = featureMatching(input, allDescriptors, numImg)
    %%***********************************************************************%
    %*                   Automatic panorama stitching                       *%
    %*                        Feature matching                              *%
    %*                                                                      *%
    %* Code author: Preetham Manjunatha                                     *%
    %* Github link: https://github.com/preethamam                           *%
    %* Date: 05/14/2024                                                     *%
    %************************************************************************%
    
    % Initialize
    matches = cell(numImg);
    matSize = size(matches);
    
    % Use symmetry and run for upper triangular matrix
    IuppeIdx = nonzeros(triu(reshape(1:numel(matches), size(matches))));
    
    % Initialize
    matches_temp = cell(1,length(IuppeIdx));         
    
    % Match features
    parfor i = 1:length(IuppeIdx)
        % IND2SUB converts from a "linear" index into individual
        % subscripts
        [ii,jj] = ind2sub(matSize, IuppeIdx(i));
        
        if (ii~=jj)
            matches_temp{i} = getMatches(input, allDescriptors{ii}, allDescriptors{jj});
        end

    end

    
    % Populate A matrix
    matches(IuppeIdx) = matches_temp;
end

%--------------------------------------------------------------------------------------------------------
% Auxillary functions
%--------------------------------------------------------------------------------------------------------
% [matches_temp] = getMatches(features1,features2)
function matches = getMatches(input,features1,features2)
    matches = matchFeatures(features1,features2, 'Method', input.Matchingmethod, ...
                            'MatchThreshold', input.Matchingthreshold, ...
                            'MaxRatio', input.Ratiothreshold, Unique=true);    
    matches = double(matches)';
end
