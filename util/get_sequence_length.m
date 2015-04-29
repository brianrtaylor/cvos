function N = get_sequence_length(sequence)
[~,~,img_path] = dataPaths(sequence);
files = dir([img_path '*.png']);
N = length(files);
end