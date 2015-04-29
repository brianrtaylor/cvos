function max_iterations = pd_max_iteration_wizard(imsize)
    n = prod(imsize);
    knots = [200*200, 320*240, 360*640, 640*480];
    niters = [5000, 8000, 16000, 25000];
    idx2 = find(knots > n,1,'first');
    idx1 = find(knots < n,1,'first');
    if isempty(idx2)
        max_iterations = niters(end); 
    elseif isempty(idx1)
        max_iterations = niters(1); 
    else
        d1 = abs(n-knots(idx1)) / abs(knots(idx2)-knots(idx1));
        d2 = abs(n-knots(idx2)) / abs(knots(idx2)-knots(idx1));
        max_iterations = round( d2*niters(idx1) + d1*niters(idx2) );
    end
    fprintf('Iteration wizard hath spoken: "you want %d iterations"\n', max_iterations);
end