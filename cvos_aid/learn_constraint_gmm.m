% USAGE:
% [local_gmm_bg, local_gmm_fg] = ...
%                 learn_constraint_gmm(constraints, I0, weights, imsize, params, layers)
%   constraints - Nx2 vector of constraints
%   I0 - image
%   weights -
%   imsize - [1x2] vector
%   params - 
%   layers - optional parameter determining whether layer mask is used or
%   not.
function [local_gmm_bg, local_gmm_fg] = ...
                learn_constraint_gmm( ...
                    constraints, groups, I0d, weights, imsize, params, layers)
    if nargin==6
        layers = [];
    end
    if ~isa(I0d,'double'), I0d = im2double(I0d); end
    
    EPS = 15; % how big of a ball to look at, for estimating GMMs.

    % create a boundary image and binarize it. we will not consider
    % occluders/occluded points in the vicinity of the edge.
    W = 0.5 * sum(weights, 3);
    W_binary = W <= 0.75;

    occluded = zeros(imsize);
    occluder = zeros(imsize);
    occluded( constraints(:,2) ) = 1;
    occluder( constraints(:,1) ) = 1;

    se = strel('disk', 1);
    % this MUST introduce MORE constrants:
    occluded_smoothed = imclose(occluded, se);
    occluded_smoothed = occluded_smoothed > 0;
    occluded_smoothed(logical(occluded)) = 1;

    occluder_smoothed = imclose(occluder, se);
    occluder_smoothed = occluder_smoothed > 0;
    occluder_smoothed(logical(occluder)) = 1;

    [YY1,XX1] = ind2sub( imsize, constraints(:,1) );
    [YY2,XX2] = ind2sub( imsize, constraints(:,2) );    

    PTS_OCCD = [YY2, XX2];
    PTS_OCCR = [YY1, XX1];
    
    [mu_fg, cov_fg, pi_fg, mu_bg, cov_bg, pi_bg] = ...
      learn_constraint_gmm_mex( I0d, occluded_smoothed, occluder_smoothed, ...
        W_binary, double(layers), PTS_OCCD, PTS_OCCR, groups, EPS, params );

    local_gmm_bg = struct;
    local_gmm_bg.mu = mu_bg;
    local_gmm_bg.cov = cov_bg;
    local_gmm_bg.pi = pi_bg;
    
    local_gmm_fg = struct;    
    local_gmm_fg.mu = mu_fg;
    local_gmm_fg.cov = cov_fg;
    local_gmm_fg.pi = pi_fg;
end