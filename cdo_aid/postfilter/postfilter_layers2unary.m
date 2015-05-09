function unary = layers2unary(layers)
    K = max(layers(:))+1;
    unary = zeros(size(layers,1), size(layers,2), K);
    for k=1:K
        unary(:,:,k) = (layers == (k-1));
    end
    unary = double(unary);
end
