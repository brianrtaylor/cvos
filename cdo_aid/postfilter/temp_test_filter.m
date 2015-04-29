layers = zeros([100, 200]);
layers(20:80,50:150)=2;
layers(79:82,149:151)=1;
layers(81,51)=1;
layers(85,55)=1;
for k=1:10
    rr = randi(100);
    cc = randi(200);

    layers(rr:min(100,rr+5), cc:min(200,cc+5) )=1;
end

SIZE_THRESHOLD = 0;
PERCENT_THRESHOLD =10;

layers_ = postfilter_filter_connected_components(layers, SIZE_THRESHOLD, PERCENT_THRESHOLD);
layers_ = postfilter_enforce_topology(layers_);

figure(1);
vl_tightsubplot(1,2,1);
imagesc(layers); axt;
vl_tightsubplot(1,2,2);
imagesc(layers_); axt;
