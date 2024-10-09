function trans_vec = calc_trans_vec(img,stdimg)
img_thresh = int8(~imbinarize(img));
corr_mat = xcorr2(img_thresh,stdimg);
[~,index] = max(corr_mat(size(img,1),:));
[i,j] = ind2sub(size(corr_mat),index);
trans_vec = [i-size(img,1),j-size(img,2)];
end