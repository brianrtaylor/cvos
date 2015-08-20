%-----------------------------------------------------------------------------
% evaOne
%
% computes some statistics between a groundtruth and a test image
%
% @return: res: results structure including...
% * p: precision
% * r: recall
% * a: accuracy in terms of intersection-union score
% * f: f-measure
%
% @param: gtImg: groundtruth label map
% @param: testImg: groundtruth label map
% @note: assumes input images are of the same dimension
%-----------------------------------------------------------------------------
function res = evaOne(gtImg, testImg)

testMask = vec(logical(testImg));
gtMask = vec(logical(gtImg));

tp = sum(gtMask & testMask);
fp = sum(~gtMask & testMask);
fn = sum(gtMask & ~testMask);
tn = sum(~gtMask & ~testMask);

res.p = tp / (tp + fp);
res.r = tp / (tp + fn);
res.a = tp / (tp + fp + fn);
res.f = 2 * ((res.p * res.r) / (res.p + res.r)); % fmeasure
end
