%----------------------------------------------------------------------------%
% evaOne
%
% computes some statistics on a single image
%----------------------------------------------------------------------------%
% eva.m just does some evaluation. eventually want to use ap, for now just
% doing some prototyping to make system work. so doing just precision per
% class and recall per class
% function res = evaOne(testBndImg, gtBndImg)
function res = evaOne(gtBndImg, testBndImg)

debug = false;

% if debug;
%   fig(10509); 
%   subplot(1,2,1); imagesc(testImg); title('test image');
%   subplot(1,2,2); imagesc(gtImg); title('ground truth');
% end
% 
% if debug;
%   fig(10509); 
%   subplot(1,5,1); imagesc(gtImg); title('ground truth');
%   subplot(1,5,2); imagesc(testImg); title('test image');
%   pause
% end

dimTest = size(testBndImg);
dimGt = size(gtBndImg);

% metrics (precision, recall, accuracy, f-measure)
testMask = vec(logical(testBndImg));
gtMask = vec(logical(gtBndImg));

tp = sum(gtMask & testMask);
fp = sum(~gtMask & testMask);
fn = sum(gtMask & ~testMask);
tn = sum(~gtMask & ~testMask);

res.p = tp / (tp + fp);
res.r = tp / (tp + fn);
res.a = tp / (tp + fp + fn); % accuracy in terms of iu score as a metric
res.f = 2 * ((res.p * res.r) / (res.p + res.r)); % fmeasure
res.fastf = 2 * tp / (2 * tp + fn + fp);
end
