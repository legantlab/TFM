nan_count = 0;
for i = 1:69
    image_temp = imread("Z:\Yu\TFMimage_trial\C2-u2os_LifeactGFP_560Beads_40xwT_wells_zstack001.tif",i);
    nan_count = nan_count + sum(isnan(image_temp),'all');
end

print(nan_count);
