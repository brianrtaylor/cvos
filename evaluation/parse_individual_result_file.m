% Parses result files of evaluation software in
%   http://lmb.informatik.uni-freiburg.de/resources/datasets/moseg.en.html
function [P, R, F, num_obj_extracted, num_obj_in_scene] = parse_individual_result_file(fname)
    fid = fopen(fname,'r');
    
    % read file until encounter the string below:
    while (1)
        str = fgetl(fid);
        if ~isempty(strfind(str, 'Average Precision, Recall, F-measure:'))
            fprintf('found string.\n');
            break;
        end
    end
	dat = fgetl(fid);
    
    fgetl(fid);
    str = fgetl(fid);
    if isempty(strfind(str, 'Visible objects in the evaluated part of the shot:'))
        fprintf('trouble, unexpected string.\n');
        fprintf('saw:\n%s\n', str);        
    end
    num_obj_in_scene = fgetl(fid);
    fgetl(fid);
    str = fgetl(fid);
    if isempty(strfind(str, 'Extracted objects (#{F-measure > 0.75}) (excluding background):'))
        fprintf('trouble, unexpected string.\n');
        fprintf('saw:\n%s\n', str);
    end
    num_obj_extracted = fgetl(fid);

    dat = str2double(strsplit(dat,' '));
    num_obj_in_scene = str2double(num_obj_in_scene);
    num_obj_extracted = str2double(num_obj_extracted);
    P = dat(1);
    R = dat(2);
    F = dat(3);
    fclose(fid);
    fprintf('%s\n', fname);
    fprintf('P:%f\nR:%f\nF:%f\n', P, R, F);
    fprintf('obj %d/%d\n', num_obj_extracted, num_obj_in_scene);
end
