function TI=transinfo(imdata,Atform,Dthresh,Pmax)
%imdata - original image
%Atform - Matrix defining 2D affine transformation
%Dthresh - anything above min(Dthresh,0) of original image will be counted in domain
%Pmax - normalize by TI Pmax. if Pmax=0, normalize by max value of imdata.

[M,N] = size(imdata);
tform = affine2d(Atform);

centerOutput = affineOutputView([M,N],tform,'BoundsStyle','CenterOutput');
Timdata = imwarp(imdata,tform,'OutputView',centerOutput);

Dim  = double(imdata);
DTim = double(Timdata);

if Pmax==0
    Pmax = max(max(Dim));
end

Acount    = sum(sum(Dim>0 & DTim>0));
integrand = Dim.*log(Dim./DTim);
TI        = sum(sum(integrand(~isinf(integrand))));


TI=TI/Acount/Pmax;


