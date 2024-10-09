function disp_vec = calc_corr_vec(img1, img2)
img_thresh_1 = int8(imbinarize(img1,'adaptive','ForegroundPolarity','dark','Sensitivity',0.7));
img_thresh_2 = int8(imbinarize(img2,'adaptive','ForegroundPolarity','dark','Sensitivity',0.7));
corr_mat = xcorr2(img_thresh_1,img_thresh_2);
% corr_mat = xcorr2(img1,img2);
[~,index] = max(corr_mat(:));
[i,j] = ind2sub(size(corr_mat),index);
disp_vec = [i-size(img1,1),j-size(img1,2)];
end


