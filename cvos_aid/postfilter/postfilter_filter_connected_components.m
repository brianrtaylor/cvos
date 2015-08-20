function layers_ = postfilter_filter_connected_components(layers, SIZE_THRESHOLD, PERCENT_THRESHOLD)
    imsize = size(layers);
    K = max(layers(:));
    
%     ccps = cell(1,K+1);
%     ccps_ds = cell(1,K+1);
    ccps = struct;
    ccps.size = [];
    ccps.layer = [];
    ccps.pixels = {};
    ccps.pixels_ds = {};
    ccps.neighbors = [];
    for k = 0:K
       im = (layers == k);
       conns = bwconncomp(im);
       N = length(conns.PixelIdxList);
       for ll=1:N
           im = zeros(imsize);
           im( conns.PixelIdxList{ll} ) = 1;
           ds = get_ccp_boundary(im);
           ccps.size = [ccps.size, numel(conns.PixelIdxList{ll})];
           ccps.pixels = [ccps.pixels, conns.PixelIdxList(ll)];
           ccps.pixels_ds = [ccps.pixels_ds, {find(ds(:))} ];
           ccps.layer = [ccps.layer, k];
       end
    end
    % ---------------------------------------------------------------------
    % construct adjacency matrix
    M = length(ccps.layer);
    ADJ = zeros(M);
    for ii=1:M
        pix = ccps.pixels{ii};
        for jj=1:M
            adjpix = ccps.pixels_ds{jj};
            if ~isempty(intersect(pix, adjpix))
                ADJ(ii,jj) = 1;
            end
        end
    end
    % ---------------------------------------------------------------------
    % for each component, find the smalles adjacent component. if the
    % smallest adjacent component is much larger than the given one, then 
    % delete the smallest..
    ccps.invalid = zeros(1,M);
    ccps.invalid = (ccps.size < SIZE_THRESHOLD);
    
    for ii=1:M
        thissz = ccps.size(ii);
        bad_idx = (ccps.invalid) | (ADJ(ii,:)==0);
        othersz = min( ccps.size + 1e10*( bad_idx ) );
        if thissz/othersz*100 < PERCENT_THRESHOLD && othersz < 1e10
            ccps.invalid(ii) = 1;
        end
    end
    % construct layers:
    layers_ = zeros(imsize) + nan;
    for ii=1:M
        if ccps.invalid(ii), continue; end;
        layers_( ccps.pixels{ii} ) = ccps.layer(ii);
    end
end

function ds = get_ccp_boundary(imbw)
    se = strel('disk',3);
    im_dilated = imdilate(imbw, se);
    ds = im_dilated > 0;
    ds(imbw==1)=0;
end