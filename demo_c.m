addpath(genpath('./func/'))
im          = imread('img/Flower.jpg');
[cartoon,texture]=decom_isotropy_color(im);
imshow([cartoon,3*texture+0.5]);
imwrite([cartoon,3*texture+0.5],'results_c.jpg')
rmpath(genpath('./func/'))