% USAGE: FG = utils_get_foreground_map(layers, area_threshold)
%   Tries to 'intelligently' create the foreground/background map, given
%   current layers.
%   'Intelligent' refers to trying to set the threshold in such a way that
%   the foreground occupies less than 'area_threshold' fraction of total
%   pixels.
%   (All this is needed to prevent everything being classified as
%   foreground, if the background layer was erroneously set to 1 due to
%   segmentation.
function FG = utils_get_foreground_map(layers, area_threshold)
    if nargin==1
        area_threshold = 0.8;
    end
    
    % threshold at 0.5 (i.e. everything above is foreground, everythign
    % below is background)
    FG = layers > 0.5;
    kk=1;
    n = numel(FG);
    while ( sum(FG(:))/n > area_threshold && kk<10 ) 
        % in case everything was labeled as foreground...
        FG = layers > (kk+0.5);
        kk=kk+1;
    end
    FG = double(FG);
end