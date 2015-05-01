function layers_ = postfilter_pff(layers, sigma, K, minsize)
    seg = seg_pff(layers, sigma, K, minsize);
    numsegs = max(seg(:));
    layers_ = layers*0;
    for k=1:numsegs
        layers_( seg==k ) = round( mean( layers( seg==k ) ) );
    end
end

function seg = seg_pff(layers, sigma, K, minsize)
    % turn layers into an RGB image.
    L = uint8( round( repmat(layers,[1 1 3])/max(layers(:))*255 ) );
    out = pff_segment_mex(L, sigma, K, minsize);
    % pff outputs a RGB colored segmentation. here convert it to
    % discrete labels.
    out_ = reshape(out, [numel(out)/3, 3]);
    [~, ~, cj] = unique(out_,'rows');
    seg = reshape(cj, size(layers));
end