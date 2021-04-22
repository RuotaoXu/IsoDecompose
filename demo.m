addpath(genpath('./func/'))
im          = imread('img/House256.png');
[cartoon,texture]=decom_isotropy(im);
imshow([cartoon,3*texture+0.5]);
imwrite([cartoon,3*texture+0.5],'results.jpg')
rmpath(genpath('./func/'))