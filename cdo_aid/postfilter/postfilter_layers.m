function layers = postfilter_layers(layers, I, lambda)
    PERCENT_THRESHOLD = 0; % percent of another.
    SIZE_THRESHOLD = 20; % pixels
    if ~isa(I,'double'), I = im2double(I); end;
    if max(layers(:)) < 0.5, 
        layers = round(layers);
        return; 
    end;
    layers_og = layers;
    
    layers_rounded = round(layers);
    pff_sigma = 0.1;
    pff_K = 10000;
    pff_minsize = 0;
    layers = postfilter_pff(layers, pff_sigma, pff_K, pff_minsize);

    %figure(20); imagesc([layers_og, layers_rounded, layers]); axt;
%     figure(2); 
%     vl_tightsubplot(1,3,1); imagesc(layers_og); axt;
%     vl_tightsubplot(1,3,2); imagesc(layers_rounded); axt;
%     vl_tightsubplot(1,3,3); imagesc(layers); axt;
    
    % removes component 'x', if the SMALLEST adjacent component 'y' size is
    % s.t. x/y*100 < PERCENT_THRESHOLD
    % or if x < SIZE_THREHSOLD
    layers = postfilter_filter_connected_components(layers, SIZE_THRESHOLD, PERCENT_THRESHOLD);
    layers = postfilter_enforce_topology(layers);
    
    fprintf('Post filtering (TV segmentation)\n');
    beta = 5;
    max_iterations = 100;
    tic;
    unary = postfilter_layers2unary(layers);

    if size(unary,3)==1
        % return the original layers (without any filtering or nans)
        layers = layers_og;       
        return;
    end;
    
    layers = bcv_tvsegment_mex( I, -unary, lambda, beta, max_iterations);
    toc;
% 
%     figure(1);
%     vl_tightsubplot(1,2,1);
%     imagesc(layers_og); axt;
%     vl_tightsubplot(1,2,2);
%     imagesc(layers); axt;
end
