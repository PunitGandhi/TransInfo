%% Tranformation Information - asymmetry analysis of flower image FlowerPinwheel.jpg
%The sections of this matlab script go through the steps for applying TI to
%identify the approximate symmetries of the provided flower in the image. 
%It also compares to other symmetry measures based on identification of the
%center of the flower and the edges of the petals.


%% Read image

tmp=imread('FlowerPinwheel.jpg');
imdata0=tmp;


%emphasize greenness of pixel, and add 1 to keep away from 0;
imdata=2*tmp(:,:,1)-tmp(:,:,2)-tmp(:,:,3)+1;

Pmax=256;
[M,N]=size(imdata);

% center of image 
pyc=(M+1)/2; 
pxc=(N+1)/2; 

%index vectors
%(1,1) is at the top left of the image
indx=1:N;
indy=1:M;

% coordinates
xcoord=indx-pxc; %positive right, same as indices
ycoord=pyc-indy; %positivue up, backwards from indices


figure;imshow(imdata)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up standard conventions for coordinates and transformations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test rotate about center 
 theta=pi/6;
 %rotate theta ccw about center of image
 Arot=[cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
 %ArotI=[cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];

tform=affine2d(Arot);
centerOutput = affineOutputView([M,N],tform,'BoundsStyle','CenterOutput');
Timdata = imwarp(imdata,tform,'OutputView',centerOutput);

figure;imshow(Timdata)

%% Test Translate
b=[500, 1000]; %left and up
%rotate theta ccw about center of image
Atb=[1,0,0;0,1,0;b(1),-b(2),1];
AtbI=[1,0,0;0,1,0;-b(1),b(2),1];

tform=affine2d(Atb);
centerOutput = affineOutputView([M,N],tform,'BoundsStyle','CenterOutput');
Timdata = imwarp(imdata,tform,'OutputView',centerOutput);

figure;imshow(Timdata)


%% Test rotate about point [b(1), b(2)] up and left froms from center
%transforms are right multiplied, so leftmost is carried out first
tform=affine2d(AtbI*Arot*Atb); %rotate about point b(1) left of and b(2) above center
centerOutput = affineOutputView([M,N],tform,'BoundsStyle','CenterOutput');
Timdata = imwarp(imdata,tform,'OutputView',centerOutput);

figure;imshow(Timdata)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify center and edge of flower by thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Threshold to create binary

%use r-g for thresholding on to indentify flower from background
imdataThr=tmp(:,:,1)-tmp(:,:,2)+1;
figure;histogram(imdataThr)
%threshold level 20 chosen by experimentation
imbin=zeros(M,N);% y-indices (rows) are first
imbin(imdataThr>=20)=1; 

figure;imshow(imbin)


%% Find Center

%input window for center of flower
cwlim=[-200, 200];  %range about center of image to consider
indCW=zeros(M,N); % y-indices (rows) are first
indCW((ycoord>cwlim(1) ) & (ycoord<cwlim(2)) ,(xcoord>cwlim(1) ) & (xcoord<cwlim(2)))=1;

%enhance yellow pixels in center of flower
imdataC=3*tmp(:,:,2)-tmp(:,:,1)-tmp(:,:,3);
%figure;imshow(imdataC)
%figure;histogram(imdataC)

%threshold 40 tuned by experimentation
imbinC=zeros(M,N);
imbinC(imdataC>=40 & (indCW))=1; %look in window for center
imbinE=zeros(M,N);
imbinE(imdataC>=40 & imbin & (~indCW) )=1; %look out of window but in flower
%figure;imshow(imbinC)


%Index meshgrid
[indX,indY]=meshgrid(indx,indy);

%coordinate meshgrid
[Xcoord, Ycoord]=meshgrid(xcoord,ycoord);


xFc=sum(Xcoord(logical(imbinC)))/sum(imbinC(:));
yFc=sum(Ycoord(logical(imbinC)))/sum(imbinC(:));


%need to subtract yFc because axis is reversed
figure;imshow(imbinE);hold on; plot(pxc+xFc,pyc-yFc,'wx');


%get shape of flower vs. angle
%first convert to polar coordinates with origin at (xFc,yFc)
[thetaE,rhoE]=cart2pol(Xcoord(logical(imbinE)) - xFc, Ycoord(logical(imbinE))-yFc);


rhoMax=1450;
indtmp=rhoE<rhoMax;
rhoE=rhoE(indtmp);
thetaE=thetaE(indtmp);



%get approximate petal sections
[~,nEmin]=min(rhoE); %find start of first petal
petalintervals=[thetaE(nEmin)]+linspace(0,2*pi,6); %approximate edges of other petal



%shift theta to align with petal intervals
thetaEtmp=thetaE;
thetaEtmp(thetaE<0)=thetaE(thetaE<0)+2*pi;
indtmp= thetaEtmp > ( 2*pi + petalintervals(1));
thetaEtmp(indtmp)=thetaEtmp(indtmp) - 2*pi;

%find petal tips
for ii=1:5
 
    thetaPtmp=petalintervals([ii, ii+1]);
    indtmp=thetaEtmp>thetaPtmp(1) & thetaEtmp<thetaPtmp(2); %theta values on petal inteval
    [rPtmp, indPtmp]=max(rhoE(indtmp ));
    rPtip(ii)=rPtmp;
    thetaPtip(ii)=thetaE(rhoE==rPtmp & indtmp);
end


figure;polarplot(thetaE,rhoE,'.');
hold on;
polarplot([petalintervals; petalintervals],[zeros(1,6); 1500*ones(1,6)],'--k'); 

polarplot([thetaPtip; thetaPtip],[zeros(1,5);rPtip],'-k','LineWidth',2); 







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure of Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ary=[1,0,0;0,-1,0;0,0,1];
Arx=[-1,0,0;0,1,0;0,0,1];
lw=2;


fig1=figure('Position',[100 100 800 200 ]);

% original image
subplot(1,3,1); imshow(imdata0); hold on;

% Image for TI analysis
ax=subplot(1,3,2); image(imdata),colormap(ax,flipud(gray(256))),axis image, axis off;


% Identified petal vectors for comparison
[xPtip,yPtip]=pol2cart(-thetaPtip,rPtip);
[xEdge,yEdge]=pol2cart(-thetaE,rhoE);
%Get Flower without background
imdataP=imdata0;
imdataP(logical(repmat(~imbin,1,1,3)))=255;
subplot(1,3,3); imshow(imdataP),axis image, axis off;hold on;
plot(pxc+xFc +[zeros(5,1), xPtip']', pyc-yFc + [zeros(5,1), yPtip']','LineStyle','-','Color','k','LineWidth',lw)
plot(pyc+xFc +xEdge, pyc-yFc + yEdge,'LineStyle','None','Color','k','LineWidth',lw,'Marker','.','MarkerSize',3)

print('FlowerImages.png','-dpng')



%%%%%%%%%%%%%%%
%% Find Center
%%%%%%%%%%%%%%%

%% Find center by rotation about x,y

%y and y shifts to consider
%Use x,y frame that starts at center of image with up and right postiive
NXmax=21; %51;
NYmax=21; %51;
yshifts=linspace(-100,400,NYmax); 
xshifts=linspace(-250,250,NXmax);

%angles of rotation to consider
dTh=2*pi/12; 
thetavals=dTh:dTh:(2*pi-dTh); %Measure angle conterclockwise from X

%downsample image
npskip=16;

b=zeros(1,2);
TIcenter=zeros(NXmax,NYmax,length(thetavals));


for llx=1:NXmax
    b(1)=xshifts(llx)/npskip;
    for kky=1:NYmax
        b(2)=yshifts(kky)/npskip;
        
        %shift right and up
        Atb=[1,0,0;0,1,0;b(1),-b(2),1];
        AtbI=[1,0,0;0,1,0;-b(1),b(2),1];
        
        for mth=1:length(thetavals)
            theta=thetavals(mth);
            %rotate theta ccw about origin
            Arot=[cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
            %ArotI=[cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];

            
            Atform=AtbI*Arot*Atb; %transforms right multiplied, so leftmost is carried out first
            TIcenter(llx,kky,mth)=transinfo(imdata(1:npskip:end,1:npskip:end),Atform,0,Pmax);
        end
        
    end
    
end

figure;imagesc(xshifts,yshifts,TIcenter(:,:,1)')
figure;imagesc(xshifts,yshifts,(max(TIcenter,[],3)-min(TIcenter,[],3))')
figure;imagesc(xshifts,yshifts,sum(abs(diff(TIcenter,[],3)).^2,3)')

%identify from plot
rxc=0;
ryc=100;


%% translate and reflect x,y to find center by TI


TIx=zeros(1,NXmax);
TIy=zeros(NYmax,1);

Ary=[1,0,0;0,-1,0;0,0,1];
Arx=[-1,0,0;0,1,0;0,0,1];

%shift in x and reflect about y
for llx=1:NXmax
        b=[xshifts(llx),0];
        Atb=[1,0,0;0,1,0;b(1),-b(2),1];
        AtbI=[1,0,0;0,1,0;-b(1),b(2),1];
        Atx=AtbI*Arx*Atb;
        
        
        TIx(llx)=transinfo(imdata,Atx,0,Pmax);
 
end

for kky=1:NYmax
        b=[0,yshifts(kky)];
        
        Atb=[1,0,0;0,1,0;b(1),-b(2),1];
        AtbI=[1,0,0;0,1,0;-b(1),b(2),1];
        Aty=AtbI*Ary*Atb;
        
        TIy(kky)=transinfo(imdata,Aty,0,Pmax);
end
    

%center will be at min of TI
[xpks,xlocs,w,proms]=findpeaks(-TIx);
figure;findpeaks(-TIx,xshifts,'Annotate','extents')
title('x center')
[~,nxc]=max(proms);
bxc=xshifts(xlocs(nxc)) %xshift in pixels

[ypks,ylocs,w,proms]=findpeaks(-TIy);
figure;findpeaks(-TIy,yshifts,'Annotate','extents')
title('y center')
[~,nyc]=max(proms);
byc=yshifts(ylocs(nyc)) %yshift in pixels

%get min and max y locations of peaks on plateau
%bycMin=min(yshifts(ylocs));
%bycMax=max(yshifts(ylocs));
byc2=120; %identify from graph


%% Find center by area method

%get area above and below 
for llx=1:NXmax
        
        txc=pxc+xshifts(llx);
        
        Aleft=sum(sum(imbin(:,indx<txc)));
        Aright=sum(sum(imbin(:,indx>txc)));
        
        Atot=Aleft+Aright;
        Adiff=abs(Aleft-Aright);
        indsi=Atot~=0;

        SIx(llx)= sum(Adiff(indsi)./Atot(indsi))/sum(indsi);
        
end

for kky=1:NYmax
    
            tyc=pyc-yshifts(kky);
        
        Atop=sum(sum(imbin(indy<tyc,:)));
        Abot=sum(sum(imbin(indy>tyc,:)));
        
        Atot=Atop+Abot;
        Adiff=abs(Atop-Abot);
        indsi=Atot~=0;

        SIy(kky)= sum(Adiff(indsi)./Atot(indsi))/sum(indsi);

   end


%center will be at min of SI
[xpks,xlocs,w,proms]=findpeaks(-SIx);
figure;findpeaks(-SIx,xshifts,'Annotate','extents')
title('x center')
[~,nxc]=max(proms);
axc=xshifts(xlocs(nxc)) %xshift in pixels

[ypks,ylocs,w,proms]=findpeaks(-SIy);
figure;findpeaks(-SIy,yshifts,'Annotate','extents')
title('y center')
[~,nyc]=max(proms);
ayc=yshifts(ylocs(nyc)) %yshift in pixels


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot comparison of centers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot of center Finding
ltcolor=[0.0 0.75 0.25];
lrcolor=[0.75 0.25 0];
lw=2;
ms=8;
fntsz=12;
plylims=[0.12 0.17];
xlabshift=-160;
ylabshift=50;

xClim=pxc+xshifts([1 end]);
yClim=pyc-yshifts([end 1]);
yClimR=yshifts([1 end]);

fig1=figure('Position',[100,100,800,235],'DefaultAxesFontSize',fntsz,'DefaultTextInterpreter','Latex');


ax=subplot(1,3,1);
hold on; imagesc(xshifts,yshifts,sum(abs(diff(TIcenter,[],3)).^2,3)');
%axis image, axis off; 
set(ax,'YDir','normal','DataAspectRatio',[1 1 1]);
line((rxc)*[1,1], yClimR,'Color',lrcolor,'Linewidth',lw/2,'LineStyle','--')
line(xClim-pxc,(ryc)*[1 1],'Color',lrcolor,'Linewidth',lw/2,'LineStyle','--')
plot((rxc),(ryc),'Color',lrcolor,'Marker','*','MarkerSize',ms)
plot((rxc),(ryc),'Color',lrcolor,'Marker','o','MarkerSize',ms)
xlim(xClim-pxc); ylim(yClimR);
xlabel('$x$ shift (pix)')
ylabel('$y$ shift (pix)')

ax=subplot(1,3,2);
hold on;
image(imdata),colormap(ax,flipud(gray(256))), 
axis image, axis off; 
set(ax,'YDir','reverse','DataAspectRatio',[1 1 1]);
xlim(xClim); ylim(yClim);

 

line(pxc*[1,1], yClim,'Color','k','Linewidth',lw,'LineStyle',':')
line(xClim,pyc*[1 1],'Color','k','Linewidth',lw,'LineStyle',':')
%rotation center
line(xClim,(pyc-ryc)*[1 1],'Color',lrcolor,'Linewidth',lw/2,'LineStyle','--')
plot((pxc+rxc),(pyc-ryc),'Color',lrcolor,'Marker','*','MarkerSize',ms)
plot((pxc+rxc),(pyc-ryc),'Color',lrcolor,'Marker','o','MarkerSize',ms)
%translation center
line(xClim,(pyc-byc)*[1 1],'Color',ltcolor,'Linewidth',lw/2,'LineStyle','--')
plot((pxc+bxc),(pyc-byc),'Color',ltcolor,'Marker','x','MarkerSize',ms)
plot((pxc+bxc),(pyc-byc),'Color',ltcolor,'Marker','o','MarkerSize',ms)
%area center
plot((pxc+axc),(pyc-ayc),'Color','k','Marker','^','MarkerSize',ms/2)
plot((pxc+axc),(pyc-ayc),'Color','k','Marker','o','MarkerSize',ms)
%manual center
plot((pxc+xFc),(pyc-yFc),'Color','k','Marker','.','MarkerSize',ms)
plot((pxc+xFc),(pyc-yFc),'Color','k','Marker','o','MarkerSize',ms)

ax=subplot(1,3,3);
yyaxis right;
hold on,plot(TIy,yshifts,'LineWidth',lw,'Color',ltcolor,'LineStyle','-')
plot(plylims,byc*[1 1],'LineWidth',lw,'Color',ltcolor,'LineStyle','--')
plot(plylims,byc2*[1 1],'LineWidth',lw,'Color',ltcolor,'LineStyle',':')
%plot(bycMax*[1 1],plylims,'LineWidth',lw,'Color',lycolor)
xlabel('$TI$'),
ylabel('$y$ shift (pix)')
xlim(plylims)
set(gca,'FontSize',fntsz);
ax.YAxis(1).Visible='off';
ax.YAxis(2).Color='k';
print('FlowerRotCenter.png','-dpng')





%% Compute Frey et al 2007 measure of asymmetry

%shift angle to be relative to longest
[rPmax ntmp]=max(rPtip);
thetaPmax=thetaPtip(ntmp);
thetaPtip2pi=mod(thetaPtip-thetaPmax,2*pi);
%thetaPtip(thetaPtip>pi)=thetaPtip(thetaPtip>pi)-2*pi;

%sort conterclockwise
[thetaPtip2pi,indtmp]=sort(thetaPtip2pi);
rPtip2pi=rPtip(indtmp);

%get avg angle and length
rPavg=sum(rPtip)/numel(rPtip);
thetaPavg=sum(thetaPtip2pi -(0:4)*2*pi/5)/numel(thetaPtip2pi);

%formula from paper
Qdist= sum(rPavg^2+rPtip.^2-2*rPavg*rPtip.*cos(thetaPtip-thetaPavg))/rPmax^2/numel(rPtip)


%% rotate tips and compare sum of squared distances


Nthmax=361;
thetavals=linspace(0,360,Nthmax)*pi/180;


ZIrot=zeros(Nthmax,1);

pTip0=[xPtip;yPtip];
pTip=pTip0;
pTip1=circshift(pTip,1,2);


for llth=1:Nthmax
    
    
    %Rotate ccw by theta by left multiply
    theta=thetavals(llth);
    ArotL=[cos(theta), -sin(theta); sin(theta), cos(theta)];
    %ArotLI=[cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
    
    TpTip=ArotL*pTip0;
    
    ZItmp=sqrt(sum(sum( (pTip-TpTip).^2 )));
    ZItmp1=sqrt(sum(sum( (pTip1-TpTip).^2 )));
    
    if ZItmp>ZItmp1
        pTip=pTip1;
        pTip1=circshift(pTip,1,2);
        ZIrot(llth)=ZItmp1;
    else
        ZIrot(llth)=ZItmp;
    end
    
    
    
end

figure;polarplot(thetavals,ZIrot/rPmax)




%% Test SI about center y axis
indtop=indy< (pyc+ayc);
indbot=indy> (pyc+ayc);


for ii=1:N
    
    Atop(ii)=sum(imbin(indtop,ii));
    Abot(ii)=sum(imbin(indbot,ii));
end

Atot=Atop+Abot;
Adiff=abs(Atop-Abot);
indsi=Atot~=0;

SI= sum(Adiff(indsi)./Atot(indsi))/sum(indsi);

%% rotate test TI

b=[bxc,byc];
Atb=[1,0,0;0,1,0;b(1),-b(2),1];
AtbI=[1,0,0;0,1,0;-b(1),b(2),1];

theta=20 *pi/180;
Arot=[cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
ArotI=[cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];


%construct affine trasformation
Atform=AtbI*Arot*Atb;
tform=affine2d(Atform);
centerOutput = affineOutputView(size(imdata),tform,'BoundsStyle','CenterOutput');

Timdata=imwarp(imdata,tform,'OutputView',centerOutput);
figure;imshow(Timdata)
%imshow(imbin)

%% rotate and reflect about center

%b=[bxc,byc];
%b=[rxc,ryc];
b=[xFc,yFc];
Atb=[1,0,0;0,1,0;b(1),-b(2),1];
AtbI=[1,0,0;0,1,0;-b(1),b(2),1];

Arx=[-1,0,0;0,1,0;0,0,1];
Ary=[1,0,0;0,-1,0;0,0,1];

Nthmax=361;
thetavals=linspace(0,360,Nthmax)*pi/180;


TIrot=zeros(Nthmax,1);
TIref=zeros(Nthmax,1);

for llth=1:Nthmax

    theta=thetavals(llth);
    Arot=[cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
    ArotI=[cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];

    
    %rotate by theta
    tform=AtbI*Arot*Atb;
    TIrot(llth)=transinfo(imdata,tform,0,Pmax);
    
    %reflect about axis at theta
    tform=AtbI*ArotI*Ary*Arot*Atb; %bring axis to x-axis, reflect y->-y, rotate back
    TIref(llth)=transinfo(imdata,tform,0,Pmax);
    
   
    
     
    
 end
    
figure;
subplot(1,2,1),polarplot(thetavals,TIrot'),title('rotation');
subplot(1,2,2),polarplot(thetavals,TIref'),title('reflection');


figure;
subplot(2,1,1),plot(thetavals*180/pi,TIrot'),title('rotation');
subplot(2,1,2),plot(thetavals*180/pi,TIref'),title('reflection');


%%
%plot Nprom-1 most prominent minima for reflections from 0 to pi
Nprm=5;

thetaexclude=[pi/2,3*pi/2];

[pks,locs,w,proms]=findpeaks(-TIref);
figure;findpeaks(-TIref,thetavals*180/pi,'Annotate','extents')

    
    [~,indprom]=sort(proms,'descend');

    figure1=figure('Position',[100 100 800 200]);
    
    kk=0;
    kp=1;
    subplot(1,Nprm,kp) ,polar(thetavals,TIref'),title('TI');
    
while kp<Nprm
    kk=kk+1;
    indk=locs(indprom(kk));
    theta=thetavals(indk);
        
    while (theta < thetaexclude(1)) || (theta > thetaexclude(2)) % exclude 0 translation
        kk=kk+1;
        indk=locs(indprom(kk));
        theta=thetavals(indk);
    end
    
    Arot=[cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
    ArotI=[cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];

    
    %reflect about axis at theta
    Atform=AtbI*ArotI*Ary*Arot*Atb;
    
    TI=TIref(indk);
    
    tform=affine2d(Atform);
    centerOutput = affineOutputView(size(imdata),tform,'BoundsStyle','CenterOutput');
    Timdata=imwarp(imdata,tform,'OutputView',centerOutput);
    
    kp=kp+1;
    subplot(1,Nprm,kp),
    colormap(gray(256)),imagesc(imdata, 'AlphaData', 0.5),
    hold on, imagesc(Timdata, 'AlphaData', 0.5),axis image;
    %title(['dir=' num2str(round(angle*180/pi)) ', pix=' num2str(round(bnorm)) ', TI=' num2str(TI)]);
    title([ 'angle=' num2str(round(theta*180/pi)) ', TI=' num2str(TI)]);
end


%%
%plot Nprom most prominent minima for rotations from 0 to 2pi
Nprm=5;

thetaexclude=[0,2*pi];

[pks,locs,w,proms]=findpeaks(-TIrot);
figure;findpeaks(-TIrot,thetavals*180/pi,'Annotate','extents')

    
    [~,indprom]=sort(proms,'descend');

    figure1=figure('Position',[100 100 800 200]);
    
    kk=0;
    kp=1;
    subplot(1,Nprm,kp) ,polar(thetavals,TIrot'),title('TI');

while kp<Nprm
    kk=kk+1;
    indk=locs(indprom(kk));
    theta=thetavals(indk);
        
    while (theta < thetaexclude(1)) || (theta > thetaexclude(2)) % exclude 0 translation
        kk=kk+1;
        indk=locs(indprom(kk));
        theta=thetavals(indk);
    end
    
    Arot=[cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
    ArotI=[cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];

    
    %rotate by angle theta
    Atform=AtbI*Arot*Atb;
    
    TI=TIrot(indk);
    
    tform=affine2d(Atform);
    centerOutput = affineOutputView(size(imdata),tform,'BoundsStyle','CenterOutput');
    Timdata=imwarp(imdata,tform,'OutputView',centerOutput);
    
    kp=kp+1;
    subplot(1,Nprm,kp),
    colormap(gray(256)),imagesc(imdata, 'AlphaData', 0.5),
    hold on, imagesc(Timdata, 'AlphaData', 0.5),axis image;
    %title(['dir=' num2str(round(angle*180/pi)) ', pix=' num2str(round(bnorm)) ', TI=' num2str(TI)]);
    title([ 'angle=' num2str(round(theta*180/pi)) ', TI=' num2str(TI)]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures of asymmetry analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Top symmetries
TopRotInd=[73 288];
TopRotRnk=[1 2];
TopRefInd=[96];
TopRefRnk=[3];


%% Plot TI and ZI, rank top symmetries

figure;
subplot(1,2,1),polarplot(thetavals,TIrot'),title('rotation');
subplot(1,2,2),polarplot(thetavals,TIref'),title('reflection');

lrfcolor=[0.0 0.75 0.25];
lrtcolor=[0.75 0.25 0];
lw=2;
ms=8;
fntsz=12;

fig1=figure('Position',[100,100,800,235],'DefaultAxesFontSize',fntsz,'DefaultTextInterpreter','Latex');

subplot(1,3,1); 
polarplot(thetavals(thetavals<pi),TIref(thetavals<pi)','LineWidth',lw/2,'Color','k','LineStyle','-');hold on;
polarplot(thetavals(thetavals>pi),TIref(thetavals>pi)','LineWidth',lw/2,'Color','k','LineStyle','--');
polarplot( [thetavals(TopRefInd); thetavals(TopRefInd)],...
           [zeros(size(TopRefInd)); 0.15],...
          'LineWidth',lw,'Color',lrfcolor,'LineStyle','-');
polarplot( [thetavals(TopRefInd); thetavals(TopRefInd)]+pi,...
           [zeros(size(TopRefInd)); 0.15],...
          'LineWidth',lw,'Color',lrfcolor,'LineStyle','--');  
text(thetavals(TopRefInd(1))+4*pi/180,0.16,'(c)','Color',lrfcolor);      
title('$TI_{ref}$');


subplot(1,3,2);
polarplot(thetavals,TIrot','LineWidth',lw/2,'Color','k'),hold on;
polarplot( [thetavals(TopRotInd); thetavals(TopRotInd)],...
           [zeros(size(TopRotInd)); 0.15*ones(size(TopRotInd))],...
          'LineWidth',lw,'Color',lrtcolor,'LineStyle','-');
text(thetavals(TopRotInd(1))+4*pi/180,0.18,'(a)','Color',lrtcolor);
text(thetavals(TopRotInd(2))-4*pi/180,0.18,'(b)','Color',lrtcolor);
title('$TI_{rot}$');

subplot(1,3,3);
polarplot(thetavals,ZIrot/rPmax,'LineWidth',lw/2,'Color','k'),title('$ZI_{rot}$')

print('FlowerTI.eps','-depsc')



%% Plot Approximate Symmetries

fig1=figure('Position',[100,100,800,235],'DefaultAxesFontSize',fntsz,'DefaultTextInterpreter','Latex');

% (a)
theta=thetavals(TopRotInd(1));
Arot=[cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
ArotI=[cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];    
 
Atform=AtbI*Arot*Atb;
tform=affine2d(Atform);
centerOutput = affineOutputView(size(imdata),tform,'BoundsStyle','CenterOutput');
Timdata=imwarp(imdata,tform,'OutputView',centerOutput);
imdiffa=double(Timdata)-double(imdata);  

Atform2=[cos(theta), sin(theta); -sin(theta), cos(theta)];
ax=subplot(1,3,1);hold on;
imagesc(imdiffa),colormap(turbo),axis image, axis off;
plot([pxc+b(1), pxc+xFc+xPtip(1) ], [pyc-b(2), pyc-yFc+yPtip(1) ],'k--')
lrot=Atform2*[xPtip(1); yPtip(1)];
plot([pxc+b(1), pxc+xFc+lrot(1) ], [pyc-b(2), pyc-yFc+lrot(2) ],'k--')
set(ax,'YDir','reverse','DataAspectRatio',[1 1 1]);


% (b)
theta=thetavals(TopRotInd(2));
Arot=[cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
ArotI=[cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];    
 
Atform=AtbI*Arot*Atb;
tform=affine2d(Atform);
centerOutput = affineOutputView(size(imdata),tform,'BoundsStyle','CenterOutput');
Timdata=imwarp(imdata,tform,'OutputView',centerOutput);
imdiffb=double(Timdata)-double(imdata);  

Atform2=[cos(theta), sin(theta); -sin(theta), cos(theta)];
ax=subplot(1,3,2);hold on;
imagesc(imdiffb),colormap(turbo),axis image, axis off;
%plot([pxc+b(1), N ], [pyc-b(2), pyc-b(2) ],'k--')
plot([pxc+b(1), pxc+xFc+xPtip(1) ], [pyc-b(2), pyc-yFc+yPtip(1) ],'k--')
lrot=Atform2*[xPtip(1); yPtip(1)];
plot([pxc+b(1), pxc+xFc+lrot(1) ], [pyc-b(2), pyc-yFc+lrot(2) ],'k--')
set(ax,'YDir','reverse','DataAspectRatio',[1 1 1]);

% (c)
theta=thetavals(TopRefInd(1));
Arot=[cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
ArotI=[cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];    
 
Atform=AtbI*ArotI*Ary*Arot*Atb;
tform=affine2d(Atform);
centerOutput = affineOutputView(size(imdata),tform,'BoundsStyle','CenterOutput');
Timdata=imwarp(imdata,tform,'OutputView',centerOutput);
imdiffc=double(Timdata)-double(imdata);  

ax=subplot(1,3,3);hold on;
imagesc(imdiffc),colormap(turbo),axis image, axis off;
plot( -([1 M ] -pyc+yFc)/tan(theta)+pxc+xFc, [1 M],'k--')
set(ax,'YDir','reverse','DataAspectRatio',[1 1 1]);

print('FlowerSymmetries.png','-dpng')

