% ensure there are no 'skipped' layers
function layers_ = postfilter_ensure_no_skipped_layers(layers)
    layers_ = layers;
    unique_layers = unique(layers(:));
    unique_layers = unique_layers(~isnan(unique_layers));
    for ii=1:length(unique_layers)
        qq = unique_layers(ii);
        layers_(layers==qq) = (ii-1);
    end
end