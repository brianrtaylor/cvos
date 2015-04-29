% Evaluation suite saves individual result files into the same folders
% where the Tracks.dat files are.
p = '/path/where/the/results/live/';

% -------------------------------------------------------------------------
TESTSET = {'cars1', 'cars4', 'cars5', 'cars10', 'marple2', 'marple4',...
'marple6', 'marple7', 'marple9', 'marple12', 'people1', 'people2',...
'tennis', 'camel1', 'cats1', 'cats3', 'cats6', 'dogs1', 'dogs2',...
'farm1', 'goats1', 'horses2', 'horses4', 'horses5', 'lion1',...
'people3', 'giraffes1', 'rabbits2', 'rabbits3', 'rabbits4'};         
              
TRAINSET = {'cars2', 'cars3', 'cars6', 'cars7', 'cars8', 'cars9', 'marple1', 'marple3', ...
'marple5', 'marple8', 'marple10', 'marple11', 'marple13', 'bear1', 'bear2', ...
'cats2', 'cats4', 'cats5', 'cats7', 'ducks1', 'horses1', 'horses3', ...
'horses6', 'lion2', 'meerkats1', 'people4', 'people5', 'rabbits1', 'rabbits5'};
% -------------------------------------------------------------------------
% Read each file.
SET = TRAINSET;
N = length(SET);
P = zeros(1,N);
R = zeros(1,N);
F = zeros(1,N);
num_obj_extracted = zeros(1,N);
num_obj = zeros(1,N);
for k=1:N
    sequence = SET{k};
    files = dir([p, sequence '/*Numbers.txt']);
    fname = files(1).name;
    fprintf('%s: %s\n', sequence, fname);
    [P(k), R(k), F(k), num_obj_extracted(k), num_obj(k)] = ...
            parse_individual_result_file([p, sequence, '/', fname]);
end
% final scores:
Pavg = mean(P);
Ravg = mean(R);
numextracted = sum(num_obj_extracted);
numobjects = sum(num_obj);
Favg = 2*Pavg*Ravg/(Pavg+Ravg);
Find = (2*P.*R)./(P+R);
fprintf('Final scores\n');
fprintf('P: %f\n', Pavg);
fprintf('R: %f\n', Ravg);
fprintf('F: %f\n', Favg);
fprintf('n: %d/%d\n', numextracted, numobjects);
