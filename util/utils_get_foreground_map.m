%-----------------------------------------------------------------------------
% utils_get_foreground_map
%
% tries to "intelligently" create the foreground/background map, given
% current layers. "Intelligent" refers to trying to set the threshold in 
% such a way that the foreground occupies less than 'area_threshold' 
% fraction of total pixels. This is handy to prevent everything from 
% being classified as foreground, for instance if the background layer was 
% erroneously set to 1 due to segmentation.
%
% @return: FG: foreground region ("intelligent" layers > 0)
% @param: layers (MxN): input layer values
% @param: area_threshold: threshold for "too much foreground"
%-----------------------------------------------------------------------------
function FG = utils_get_foreground_map(layers, area_threshold)
if ~exist('area_threshold'); area_threshold = 0.8; end

% threshold at 0.5 (i.e. everything above is foreground, below is background)
FG = layers > 0.5;
kk=1;
n = numel(FG);
while ( sum(FG(:)) / n > area_threshold && kk < 10 ) 
  % in case everything was labeled as foreground
  FG = layers > (kk + 0.5);
  kk = kk + 1;
end
FG = double(FG);
end
