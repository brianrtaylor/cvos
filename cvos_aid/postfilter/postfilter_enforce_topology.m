% every region assigned to a nonzero layer (i.e. c=1,2,3,...) must be 
% connected to a region with layer value smaller by 1.
function layers_ = postfilter_enforce_topology(layers)
    imsize = size(layers);
    layers_ = layers;
    nlayers = max(layers(:));
    se = strel('disk',3);
    for ii=1:nlayers
        L = (layers==ii);
        conns = bwconncomp( L );
        N = length(conns.PixelIdxList);
        for ll=1:N
            im = zeros(imsize);
            im( conns.PixelIdxList{ll} ) = 1;
            im_dilated = imdilate(im, se);
            ds = im_dilated > 0;
            ds(im==1)=0;
            % which layer values compose the boundary of the region?
            layer_vals = unique(layers(ds>0));
            layer_vals = layer_vals(~isnan(layer_vals));
            layer_vals = layer_vals(layer_vals < ii);
            maxlayerval = max(layer_vals); % largest value strictly smaller than current layer.
            if isempty(layer_vals), continue; end;
            if (maxlayerval < ii-1) && maxlayerval>=0
                layers_(im>0) = maxlayerval+1;
                fprintf('enforcing depth-layer topology of a ccp...\n');
            end
        end
    end
    
    % smallest value is always zero.
    layers_ = layers_ - min(layers_(:));
end