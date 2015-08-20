% convenience function to find the length of a sequence by counting the
% number of image files in the directory of that particular sequence.
% (relies on each sequence's image files resting in its own directory)
function N = get_sequence_length(sequence)
[~,~,img_path] = dataPaths(sequence);
files = dir([img_path '*.png']);
N = length(files);
end
