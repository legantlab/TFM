window = 7;
blur_kernel = ones(window)/window^2;
blur_img = imfilter(imgsection_1,blur_kernel);
figure
imagesc(blur_img);