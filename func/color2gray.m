function I=color2gray(im)
if(size(im,3)>1)
    I=rgb2gray(im);
else
    I=im;
end