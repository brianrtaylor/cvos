function groups = utils_group_constraints(constraints, imsize, group_size)
    groups = utils_group_constraints_mex(constraints, imsize, group_size);
    % NOTE: the code below was used in the publication.
    % It is somewhat slow when the number of constraints is large... the
    % code above performs faster, but may produce slightly different
    % results.
end

% % occr = constraints(:,1);
% % occd = constraints(:,2);
% % [yy_occr, xx_occr] = ind2sub(imsize, occr);
% % [yy_occd, xx_occd] = ind2sub(imsize, occd);
% % 
% % X = double([yy_occr, xx_occr, yy_occd, xx_occd])';
% % numcenters = ceil(size(X,2)/group_size);
% % if numcenters < size(X,2)
% %     [~, groups] = vl_kmeans(X, numcenters);
% %     groups = verify_groups(X, groups, group_size);
% % else
% %     groups = 1:size(X,2);
% %     fprintf('bypassing clustering (with %d groups)\n', size(X,2) );
% % end
% % groups = double(groups);
% % end
% % 
% % function groups = verify_groups(X, groups, group_size)
% %     ngroups = max(groups);
% %     lastgroupsofar=ngroups;
% %     
% %     dist_threshold = sqrt(2)* (group_size*2);
% %     
% %     for grp=1:ngroups
% %         this = (groups==grp);
% %         D = vl_alldist( X(:,this) );
% %         maxdist = sqrt( max(D(:)) );
% %         if maxdist > dist_threshold
% %             % declare group invalid and split it into individual points..
% %             idx = find(this);
% %             for k = idx(2:end)
% %                 groups(k) = lastgroupsofar+1;
% %                 lastgroupsofar=lastgroupsofar+1;
% %             end        
% %         end
% %     end
% %     
% % end