function [wei,idx]=gpuMatch6(im,opts) 
im=single(im);
%% Parameters
blksize=opts.blksize;
blkradius=floor(blksize/2);
if(isfield(opts,'mask'))
    mask=opts.mask;
    midY=ceil(size(mask,1)/2);
    midX=ceil(size(mask,2)/2);
    neighbor=sum(mask(:)>0);
else
    radius=opts.radius;
    win=radius*2+1;
    mask=ones(win,win);    
    midY=ceil(size(mask,1)/2);
    midX=ceil(size(mask,2)/2);
    neighbor=sum(mask(:)>0);
end

ext_mode=opts.ext_mode;

nblk=opts.nblk;

height=size(im,1);
width=size(im,2);
%% Padding
switch ext_mode
    case 'symmetric'
        pad=padarray(im,[blkradius blkradius],'symmetric','both');
    otherwise
        error('Unknown extension mode!!!\n');
end
%% Pixel index
gpuX=reshape(gpuArray.colon(1,width),[1 1 width]);
gpuY=reshape(gpuArray.colon(1,height),[1 height]);
%% Neighborhood index
[DX,DY]=ind2sub(size(mask),find(mask));
DY=DY-midY;
DX=DX-midX;
gpuDX=gpuArray(DX);
gpuDY=gpuArray(DY);
%% Block summation
gpuRes=zeros(neighbor,height,width,'single','gpuArray');
for i_sub=0:blksize-1
    for j_sub=0:blksize-1
        gpuRes=gpuRes+arrayfun(@pixelDif,gpuY,gpuX,gpuDX,gpuDY);
    end
end
gpuRes=reshape(gpuRes,[neighbor,height*width]);
res=gather(gpuRes);
%% Sorting
% [wei,idx]=sort(res,win*win);
if(exist('mink','builtin'))
    [wei,idx]=mink(res,nblk);
else
    [wei,idx]=sort(res,1,'ascend');
    wei=wei(1:nblk,:);
    idx=idx(1:nblk,:);
end
%% Indexing
T=idxTranslate_(size(im),DX,DY);
idx=T(idx)+repmat(1:height*width,[size(idx,1),1]);
%% Sub functions
    %% pixel differences
    function res=pixelDif(r,c,x,y)
        if(r+y>0&&r+y<=height&&c+x>0&&c+x<=width)
            res=(pad(r+i_sub,c+j_sub)-pad(r+y+i_sub,c+x+j_sub)).^2;
        else
            res=single(inf);
        end
    end
    %%
end
function T=idxTranslate_(imsize,DX,DY)
T=DX*imsize(1)+DY;
T=T(:);
end