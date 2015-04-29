% out = erode_using_msfm(bwimg, erodesz)
function out = erode_using_msfm(bwimg, erodesz)
[ii, jj] = find(bwimg == 0);
out = msfm(double(bwimg), [ii, jj]', true, true) > erodesz;
end