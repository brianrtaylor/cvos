%-----------------------------------------------------------------------------
% write_evaluation_tracks_files
%
% prior to running the evaluation suite, convert matlab results into
% ASCII *.dat files (input to evaluation suite)
%-----------------------------------------------------------------------------
for s = 1:length(SEQUENCES);
  sequence = SEQUENCES{s};
  [~, ~, img_path, ~, extension] = dataPaths(sequence);
  %---------------------------------------------------------------------------
  load('/path/where/segmentation/lives/%s/segmentation.mat', sequence);
  
  T = size(segmentation,3);  
  out_fname = sprintf('/path/where/results/live/%s/Tracks%d.dat', sequence, T);
  mkdir( sprintf('/path/where/results/live/%s/', sequence) );
  cmd = sprintf('rm /path/where/results/live/%s/*.dat', sequence);
  unix(cmd);
  write_evaluation_track_file_mex(segmentation, out_fname);
  %---------------------------------------------------------------------------
end
