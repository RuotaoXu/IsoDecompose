%stationary frame transform for 2D image
%im---image
%order---spline order: 1 or 3
%level---decomposition level
function [wfilt, wfilt_cell, wfilt_norm] = swt2_filt(order,level)

%get the mask
allmask = splmask(order,1,'U');
[maskno, ~] = size(allmask);

%space allocation
%|  |   |   |...| ......|   |    |...|
%app    det(le)           det(1)
%ie, first array is approx part, followed by detail part on level le,
%then detail part on level (le-1) all the way up to detail on level 1
wfilt_cell = cell(1+level*(maskno^2-1),1);
for le=1:level
    if le == 1
        mask0 = 1;
        s0 = 1;
    else
        mask0 = kron(allmask(1,:)', allmask(1,:));
        s0 = s0*norm(allmask(1,:))*norm(allmask(1,:));
    end
    allmask = splmask(order,le,'U');
    [maskno, ~] = size(allmask);
    windex = (level-le)*(maskno^2-1);
    for fi=1:maskno
        for fj=1:maskno            
            wfilt_cell{windex+maskno*(fi-1)+fj} = conv2(mask0,kron(allmask(fi,:)', allmask(fj,:)));
            wfilt_norm(windex+maskno*(fi-1)+fj) = norm(allmask(fi,:))*norm(allmask(fj,:))*s0;            
        end
    end
end

max_len = 0;
for i = 1: numel(wfilt_cell)
    max_len = max(max_len,sqrt(numel(wfilt_cell{i})));
end
wfilt = zeros(max_len^2,numel(wfilt_cell));
for i = 1: numel(wfilt_cell)
    wf = wfilt_cell{i}';
    wf = padarray(wf, [1 1]*(max_len-size(wf,1)),'post');
    wfilt(:,i) = wf(:);
end

%=========%