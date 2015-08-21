%-----------------------------------------------------------------------------
% cvos_show_timing_info
%
% show timing information
%
% @return: total_time: total time everything took to run
% @return: rest_time: total time everything BUT SOLVER took to run
% @param: tm: struct with timing info
%-----------------------------------------------------------------------------
function [total_time, rest_time] = cvos_show_timing_info(tm)
rest_time = tm.preprocess + tm.layersnow + tm.box_probs ...
  + tm.constraints_now + tm.weights_now + tm.constraint_weights_now ...
  + tm.weights_propagate + tm.unity + tm.constraints_propagate ...
  + tm.fgprior + tm.problem_setup ...
  + tm.postfilter + tm.past + tm.box_models + tm.objects_update ...
  + tm.object_uvf + tm.vis + tm.flowocc_cbf + tm.flowocc_bf;
total_time = rest_time + tm.problem_solve;

fprintf('--------- timing info ---------\n');
fprintf('total time:              %02.3f\n', total_time);
fprintf('solver:                  %02.3f\n', tm.problem_solve);
fprintf('rest:                    %02.3f\n', rest_time);
fprintf('--------- ----------- ---------\n');
fprintf('preprocessing images:    %02.3f\n', tm.preprocess);
fprintf('flow & occlusions cbf:   %02.3f\n', tm.flowocc_cbf);
fprintf('flow & occlusions bf:    %02.3f\n', tm.flowocc_bf);
fprintf('layers warped to now:    %02.3f\n', tm.layersnow);
fprintf('calc prob box images:    %02.3f\n', tm.box_probs);
fprintf('constraints now:         %02.3f\n', tm.constraints_now);
fprintf('weights now:             %02.3f\n', tm.weights_now);
fprintf('constraint weight now:   %02.3f\n', tm.constraint_weights_now);
fprintf('weights propagation:     %02.3f\n', tm.weights_propagate);
fprintf('unity prior:             %02.3f\n', tm.unity);
fprintf('constraints weight prop: %02.3f\n', tm.constraints_propagate);
fprintf('foreground prior:        %02.3f\n', tm.fgprior);
fprintf('problem setup:           %02.3f\n', tm.problem_setup);
fprintf('problem solving:         %02.3f\n', tm.problem_solve);
fprintf('postfiltering:           %02.3f\n', tm.postfilter);
fprintf('update past:             %02.3f\n', tm.past);
fprintf('update box models:       %02.3f\n', tm.box_models);
fprintf('update objects:          %02.3f\n', tm.objects_update);
fprintf('update object uvf maps:  %02.3f\n', tm.object_uvf);
fprintf('visualization + saving:  %02.3f\n', tm.vis);
fprintf('----------- ----------- -------\n');
end
