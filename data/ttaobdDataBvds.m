function [problem, params] = ttaobdDataBvds(DATA, paramsIn)

outpath = '/pad_local/btay/projects/cdo/';

%----------------------------------------------------------------------------%
% Sequence Specific 
%----------------------------------------------------------------------------%

% only ones that work right now are as follows:
% Test:  8001, 8002, 8003, 8004, 8005, 8006, 8009, 8016, 8029, 8033, 8039, 8041, 8058, 8059
% Train: 8078, 8080, 8089, 8178 

% new: 8041 (planet_earth1), 8058 (vw), 8080 (kia), 8089 (sailing)

%----------------------------------------------------------------------------%
% Test 
%----------------------------------------------------------------------------%
switch DATA;
case 8001; id = 123; WB = 59; WF = 59; seq = 'airplane';
case 8002; id = 126; WB = 59; WF = 59; seq = 'angkor_wat';
case 8003; id = 120; WB = 41; WF = 41; seq = 'animal_chase';
case 8004; id = 116; WB = 14; WF = 14; seq = 'arctic_kayak';
case 8005; id = 071; WB = 59; WF = 59; seq = 'ballet';
case 8006; id = 170; WB = 59; WF = 59; seq = 'baseball';
case 8007; id = 161; WB = 59; WF = 59; seq = 'beach_volleyball';
case 8008; id = 114; WB = 59; WF = 59; seq = 'belly_dancing';
case 8009; id = 120; WB = 59; WF = 59; seq = 'beyonce';
case 8010; id = 128; WB = 59; WF = 59; seq = 'bicycle_race';
case 8011; id = 120; WB = 59; WF = 59; seq = 'birds_of_paradise';
case 8012; id = 120; WB = 59; WF = 59; seq = 'buck';
case 8013; id = 161; WB = 59; WF = 59; seq = 'buffalos';
case 8014; id = 168; WB = 59; WF = 59; seq = 'capoeira';
case 8015; id = 136; WB = 59; WF = 59; seq = 'chameleons';
case 8016; id = 096; WB = 59; WF = 59; seq = 'fisheye';
case 8017; id = 120; WB = 29; WF = 29; seq = 'fish_underwater';
case 8018; id = 068; WB = 51; WF = 51; seq = 'freight_train';
case 8019; id = 090; WB = 59; WF = 59; seq = 'frozen_lake';
case 8020; id = 105; WB = 59; WF = 59; seq = 'gokart';
case 8021; id = 040; WB = 28; WF = 28; seq = 'harley_davidson';
case 8022; id = 120; WB = 59; WF = 59; seq = 'hockey';
case 8023; id = 120; WB = 59; WF = 59; seq = 'hockey_goals';
case 8024; id = 175; WB = 59; WF = 59; seq = 'horse_gate';
case 8025; id = 118; WB = 59; WF = 59; seq = 'hummingbird';
case 8026; id = 103; WB = 59; WF = 59; seq = 'humpback';
case 8027; id = 097; WB = 59; WF = 59; seq = 'jungle_cat';
case 8028; id = 137; WB = 59; WF = 59; seq = 'kangaroo_fighting';
case 8029; id = 120; WB = 59; WF = 59; seq = 'kim_yu_na';
case 8030; id = 123; WB = 54; WF = 54; seq = 'koala';
case 8031; id = 054; WB = 53; WF = 53; seq = 'monkeys_behind_fence';
case 8032; id = 103; WB = 59; WF = 59; seq = 'nba_commercial';
case 8033; id = 104; WB = 59; WF = 59; seq = 'new_york';
case 8034; id = 147; WB = 59; WF = 59; seq = 'nordic_skiing';
case 8035; id = 161; WB = 59; WF = 59; seq = 'octopus';
case 8036; id = 140; WB = 59; WF = 59; seq = 'palm_tree';
case 8037; id = 121; WB = 59; WF = 59; seq = 'panda';
case 8038; id = 111; WB = 59; WF = 59; seq = 'panda_cub';
case 8039; id = 034; WB = 33; WF = 33; seq = 'penguins';
case 8040; id = 120; WB = 59; WF = 59; seq = 'pepsis_wasps';
case 8041; id = 120; WB = 06; WF = 06; seq = 'planet_earth_1';
case 8042; id = 120; WB = 28; WF = 28; seq = 'planet_earth_2';
case 8043; id = 084; WB = 42; WF = 42; seq = 'riverboat';
case 8044; id = 114; WB = 14; WF = 14; seq = 'rock_climbing';
case 8045; id = 153; WB = 38; WF = 38; seq = 'rock_climbing2';
case 8046; id = 100; WB = 59; WF = 59; seq = 'salsa';
case 8047; id = 138; WB = 59; WF = 59; seq = 'samba_kids';
case 8048; id = 075; WB = 59; WF = 59; seq = 'shark_attack';
case 8049; id = 146; WB = 59; WF = 59; seq = 'sled_dog_race';
case 8050; id = 124; WB = 59; WF = 59; seq = 'slow_polo';
case 8051; id = 120; WB = 59; WF = 59; seq = 'snowboarding';
case 8052; id = 144; WB = 47; WF = 47; seq = 'snowboarding_crashes';
case 8053; id = 136; WB = 59; WF = 59; seq = 'snow_leopards';
case 8054; id = 063; WB = 44; WF = 44; seq = 'street_food';
case 8055; id = 120; WB = 59; WF = 59; seq = 'swimming';
case 8056; id = 157; WB = 35; WF = 35; seq = 'up_dug';
case 8057; id = 121; WB = 32; WF = 32; seq = 'up_trailer';
case 8058; id = 120; WB = 21; WF = 04; seq = 'vw_commercial';
case 8059; id = 140; WB = 59; WF = 59; seq = 'white_tiger';
case 8060; id = 152; WB = 59; WF = 59; seq = 'yosemite';

case 8121; id = 040; WB = 28; WF = 28; seq = 'harley_kr';

%----------------------------------------------------------------------------%
% Train
%----------------------------------------------------------------------------%
case 8061; id = 160; WB = 39; WF = 39; seq = 'alec_baldwin';
case 8062; id = 090; WB = 59; WF = 59; seq = 'anteater';
case 8063; id = 180; WB = 59; WF = 59; seq = 'avalanche';
case 8064; id = 090; WB = 41; WF = 41; seq = 'big_wheel';
case 8065; id = 090; WB = 59; WF = 59; seq = 'bowling';
case 8066; id = 090; WB = 59; WF = 59; seq = 'campanile';
case 8067; id = 072; WB = 59; WF = 59; seq = 'car_jump';
case 8068; id = 072; WB = 59; WF = 59; seq = 'chrome';
case 8069; id = 090; WB = 59; WF = 59; seq = 'deoksugung';
case 8070; id = 090; WB = 59; WF = 59; seq = 'dominoes'; % was last one running
case 8071; id = 090; WB = 59; WF = 59; seq = 'drone';
case 8072; id = 090; WB = 59; WF = 59; seq = 'excavator';
case 8073; id = 075; WB = 59; WF = 59; seq = 'floorhockey';
case 8074; id = 112; WB = 30; WF = 30; seq = 'galapagos';
case 8075; id = 065; WB = 59; WF = 59; seq = 'gray_squirrel';
case 8076; id = 180; WB = 59; WF = 59; seq = 'guitar';
case 8077; id = 075; WB = 59; WF = 59; seq = 'hippo_fight';
case 8078; id = 075; WB = 59; WF = 59; seq = 'horse_riding';
case 8079; id = 170; WB = 59; WF = 59; seq = 'juggling';
case 8080; id = 076; WB = 24; WF = 24; seq = 'kia_commercial';
case 8081; id = 180; WB = 59; WF = 59; seq = 'knot';
case 8082; id = 088; WB = 59; WF = 59; seq = 'lion';
case 8083; id = 090; WB = 59; WF = 59; seq = 'lion2';
case 8084; id = 090; WB = 59; WF = 59; seq = 'lukla_airport';
case 8085; id = 090; WB = 59; WF = 59; seq = 'pouring_tea';
case 8086; id = 078; WB = 29; WF = 29; seq = 'rock_climbing'; % CAREFUL, THIS ALSO EXISTS WITH THIS NAME IN THE TESTING DATA. EVEN A ROCK_CLIMBING2 EXISTS
case 8087; id = 080; WB = 49; WF = 49; seq = 'roller_coaster';
case 8088; id = 090; WB = 59; WF = 59; seq = 'rolling_pin';
case 8089; id = 068; WB = 56; WF = 56; seq = 'sailing'; % stolen from pad/avinash/cvpr13 --> ayvaci-kr
case 8090; id = 075; WB = 59; WF = 59; seq = 'sea_snake';
case 8091; id = 090; WB = 59; WF = 59; seq = 'sea_turtle';
case 8092; id = 090; WB = 31; WF = 31; seq = 'sitting_dog';
case 8093; id = 065; WB = 59; WF = 59; seq = 'snow_shoes';
case 8094; id = 100; WB = 49; WF = 49; seq = 'soccer';
% case 8094; id = 50; WB = 49; WF = 49; seq = 'soccer';
case 8095; id = 090; WB = 59; WF = 59; seq = 'space_shuttle';
case 8096; id = 090; WB = 59; WF = 59; seq = 'swing';
case 8097; id = 120; WB = 59; WF = 59; seq = 'tarantula';
case 8098; id = 075; WB = 59; WF = 59; seq = 'tennis';
case 8099; id = 090; WB = 59; WF = 59; seq = 'trampoline';
case 8100; id = 110; WB = 43; WF = 43; seq = 'zoo';

case 8178; id = 075; WB = 59; WF = 59; seq = 'horse_riding_huge';
end

switch DATA;
case 8001; iid = 123; WWB = 1; WWF = 1; seq = 'airplane';
case 8002; iid = 126; WWB = 2; WWF = 2; seq = 'angkor_wat';
case 8003; iid = 120; WWB = 3; WWF = 3; seq = 'animal_chase';
case 8004; iid = 116; WWB = 1; WWF = 1; seq = 'arctic_kayak';
case 8005; iid = 071; WWB = 1; WWF = 1; seq = 'ballet';
case 8006; iid = 170; WWB = 1; WWF = 1; seq = 'baseball';
case 8007; iid = 161; WWB = 1; WWF = 1; seq = 'beach_volleyball';
case 8008; iid = 114; WWB = 1; WWF = 1; seq = 'belly_dancing';
case 8009; iid = 120; WWB = 1; WWF = 1; seq = 'beyonce';
case 8010; iid = 128; WWB = 1; WWF = 1; seq = 'bicycle_race';
case 8011; iid = 120; WWB = 1; WWF = 1; seq = 'birds_of_paradise';
case 8012; iid = 120; WWB = 1; WWF = 1; seq = 'buck';
case 8013; iid = 161; WWB = 1; WWF = 1; seq = 'buffalos';
case 8014; iid = 168; WWB = 1; WWF = 1; seq = 'capoeira';
case 8015; iid = 136; WWB = 1; WWF = 1; seq = 'chameleons';
case 8016; iid = 096; WWB = 1; WWF = 1; seq = 'fisheye';
case 8017; iid = 120; WWB = 1; WWF = 1; seq = 'fish_underwater';
case 8018; iid = 068; WWB = 1; WWF = 1; seq = 'freight_train';
case 8019; iid = 090; WWB = 1; WWF = 1; seq = 'frozen_lake';
case 8020; iid = 105; WWB = 1; WWF = 1; seq = 'gokart';
case 8021; iid = 040; WWB = 1; WWF = 1; seq = 'harley_davidson';
case 8022; iid = 120; WWB = 1; WWF = 1; seq = 'hockey';
case 8023; iid = 120; WWB = 1; WWF = 1; seq = 'hockey_goals';
case 8024; iid = 175; WWB = 1; WWF = 1; seq = 'horse_gate';
case 8025; iid = 118; WWB = 1; WWF = 1; seq = 'hummingbird';
case 8026; iid = 103; WWB = 1; WWF = 1; seq = 'humpback';
case 8027; iid = 097; WWB = 1; WWF = 1; seq = 'jungle_cat';
case 8028; iid = 137; WWB = 1; WWF = 1; seq = 'kangaroo_fighting';
case 8029; iid = 120; WWB = 1; WWF = 1; seq = 'kim_yu_na';
case 8030; iid = 123; WWB = 1; WWF = 1; seq = 'koala';
case 8031; iid = 054; WWB = 1; WWF = 1; seq = 'monkeys_behind_fence';
case 8032; iid = 103; WWB = 1; WWF = 1; seq = 'nba_commercial';
case 8033; iid = 104; WWB = 1; WWF = 1; seq = 'new_york';
case 8034; iid = 147; WWB = 1; WWF = 1; seq = 'nordic_skiing';
case 8035; iid = 161; WWB = 1; WWF = 1; seq = 'octopus';
case 8036; iid = 140; WWB = 1; WWF = 1; seq = 'palm_tree';
case 8037; iid = 121; WWB = 1; WWF = 1; seq = 'panda';
case 8038; iid = 111; WWB = 1; WWF = 1; seq = 'panda_cub';
case 8039; iid = 034; WWB = 5; WWF = 5; seq = 'penguins';
case 8040; iid = 120; WWB = 1; WWF = 1; seq = 'pepsis_wasps';
case 8041; iid = 120; WWB = 1; WWF = 1; seq = 'planet_earth_1';
case 8042; iid = 120; WWB = 1; WWF = 1; seq = 'planet_earth_2';
case 8043; iid = 084; WWB = 1; WWF = 1; seq = 'riverboat';
case 8044; iid = 114; WWB = 1; WWF = 1; seq = 'rock_climbing';
case 8045; iid = 153; WWB = 1; WWF = 1; seq = 'rock_climbing2';
case 8046; iid = 100; WWB = 1; WWF = 1; seq = 'salsa';
case 8047; iid = 138; WWB = 1; WWF = 1; seq = 'samba_kids';
case 8048; iid = 075; WWB = 1; WWF = 1; seq = 'shark_attack';
case 8049; iid = 146; WWB = 1; WWF = 1; seq = 'sled_dog_race';
case 8050; iid = 124; WWB = 1; WWF = 1; seq = 'slow_polo';
case 8051; iid = 120; WWB = 1; WWF = 1; seq = 'snowboarding';
case 8052; iid = 144; WWB = 1; WWF = 1; seq = 'snowboarding_crashes';
case 8053; iid = 136; WWB = 1; WWF = 1; seq = 'snow_leopards';
case 8054; iid = 063; WWB = 1; WWF = 1; seq = 'street_food';
case 8055; iid = 120; WWB = 1; WWF = 1; seq = 'swimming';
case 8056; iid = 157; WWB = 1; WWF = 1; seq = 'up_dug';
case 8057; iid = 121; WWB = 1; WWF = 1; seq = 'up_trailer';
case 8058; iid = 120; WWB = 1; WWF = 1; seq = 'vw_commercial';
case 8059; iid = 140; WWB = 1; WWF = 1; seq = 'white_tiger';
case 8060; iid = 152; WWB = 1; WWF = 1; seq = 'yosemite';

case 8121; iid = 040; WWB = 1; WWF = 1; seq = 'harley_kr';

%----------------------------------------------------------------------------%
% Train
%----------------------------------------------------------------------------%
case 8061; iid = 160; WWB = 1; WWF = 1; seq = 'alec_baldwin';
case 8062; iid = 090; WWB = 1; WWF = 1; seq = 'anteater';
case 8063; iid = 180; WWB = 1; WWF = 1; seq = 'avalanche';
case 8064; iid = 090; WWB = 1; WWF = 1; seq = 'big_wheel';
case 8065; iid = 090; WWB = 1; WWF = 1; seq = 'bowling';
case 8066; iid = 090; WWB = 1; WWF = 1; seq = 'campanile';
case 8067; iid = 072; WWB = 1; WWF = 1; seq = 'car_jump';
case 8068; iid = 072; WWB = 1; WWF = 1; seq = 'chrome';
case 8069; iid = 090; WWB = 1; WWF = 1; seq = 'deoksugung';
case 8070; iid = 090; WWB = 1; WWF = 1; seq = 'dominoes'; % was last one running
case 8071; iid = 090; WWB = 1; WWF = 1; seq = 'drone';
case 8072; iid = 090; WWB = 1; WWF = 1; seq = 'excavator';
case 8073; iid = 075; WWB = 1; WWF = 1; seq = 'floorhockey';
case 8074; iid = 112; WWB = 1; WWF = 1; seq = 'galapagos';
case 8075; iid = 065; WWB = 1; WWF = 1; seq = 'gray_squirrel';
case 8076; iid = 180; WWB = 1; WWF = 1; seq = 'guitar';
case 8077; iid = 075; WWB = 1; WWF = 1; seq = 'hippo_fight';
case 8078; iid = 075; WWB = 1; WWF = 1; seq = 'horse_riding';
case 8079; iid = 170; WWB = 1; WWF = 1; seq = 'juggling';
case 8080; iid = 076; WWB = 1; WWF = 1; seq = 'kia_commercial';
case 8081; iid = 180; WWB = 1; WWF = 1; seq = 'knot';
case 8082; iid = 088; WWB = 1; WWF = 1; seq = 'lion';
case 8083; iid = 090; WWB = 1; WWF = 1; seq = 'lion2';
case 8084; iid = 090; WWB = 1; WWF = 1; seq = 'lukla_airport';
case 8085; iid = 090; WWB = 1; WWF = 1; seq = 'pouring_tea';
case 8086; iid = 078; WWB = 1; WWF = 1; seq = 'rock_climbing'; % CAREFUL, THIS ALSO EXISTS WITH THIS NAME IN THE TESTING DATA. EVEN A ROCK_CLIMBING2 EXISTS
case 8087; iid = 080; WWB = 1; WWF = 1; seq = 'roller_coaster';
case 8088; iid = 090; WWB = 1; WWF = 1; seq = 'rolling_pin';
case 8089; iid = 068; WWB = 1; WWF = 1; seq = 'sailing'; % Flow stolen from pad/avinash/cvpr13 --> ayvaci-kr
case 8090; iid = 075; WWB = 1; WWF = 1; seq = 'sea_snake';
case 8091; iid = 090; WWB = 1; WWF = 1; seq = 'sea_turtle';
case 8092; iid = 090; WWB = 1; WWF = 1; seq = 'sitting_dog';
case 8093; iid = 065; WWB = 1; WWF = 1; seq = 'snow_shoes';
case 8094; iid = 100; WWB = 1; WWF = 1; seq = 'soccer';
% case 8094; iid = 050; WWB = 1; WWF = 1; seq = 'soccer';
case 8095; iid = 090; WWB = 1; WWF = 1; seq = 'space_shuttle';
case 8096; iid = 090; WWB = 1; WWF = 1; seq = 'swing';
case 8097; iid = 120; WWB = 1; WWF = 1; seq = 'tarantula';
case 8098; iid = 075; WWB = 1; WWF = 1; seq = 'tennis';
case 8099; iid = 090; WWB = 1; WWF = 1; seq = 'trampoline';
case 8100; iid = 110; WWB = 1; WWF = 1; seq = 'zoo';

case 8178; iid = 075; WWB = 1; WWF = 1; seq = 'horse_riding_huge';
end

%----------------------------------------------------------------------------%
% data path 
%----------------------------------------------------------------------------%
dataset = 'bvds'; ext = 'png';
dpath = sprintf('/plot/vasiliy/CVPR15/data/%s/', seq);
% dpath = '/plot/vasiliy/CVPR15/data/';
uvpath  = sprintf('/plot/vasiliy/CVPR15/data/%s/flow/', seq); % ayvaci flow

nameStr = '%s_%03d';
% imgNameStr = ['%s/image%03d.' ext];
imgNameStr = [nameStr '.' ext];
uvNameStr = [nameStr '.mat'];

expName = 'bvds';
if exist('paramsIn', 'var'); expName = paramsIn.expName; end;

%----------------------------------------------------------------------------%
% dataset parameters
%----------------------------------------------------------------------------%
params.LAMBDA = 2.0;
params.TAU = 2;
params.OCCPROB = 0.50;

%----------------------------------------------------------------------------%
% Outputs
%----------------------------------------------------------------------------%
if exist('paramsIn', 'var');
  params = structmerge(params, paramsIn);
end

% build problem for dataset
problem = v2struct(seq, ext, id, WB, WF, iid, WWB, WWF, ...
  nameStr, imgNameStr, uvNameStr, dpath, uvpath, dataset, expName);
end
