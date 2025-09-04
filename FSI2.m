clc
clear all
close all

% addpath('C:\Users\aj0036\OneDrive\Downloads\Telegram Desktop\139912140');
% cd 'C:\Users\aj0036\OneDrive\Downloads\Telegram Desktop\139912140'


addpath('C:\Users\mofid\Downloads\1-Code2\Code\139912140');
cd 'C:\Users\mofid\Downloads\1-Code2\Code\139912140'


scans = avantesSpectrumRead;
bscans = scans;
cc = 'gauss4';

N = length(scans) ;  % Total number of files

for i = 1 : N
    bscans(i).signal(1050:1090)=1550;
    bscans(i).signal(1525:1580)=1330;
end

figure; for kk=1:10 plot(bscans(kk).wvl, bscans(kk).signal); hold on; end

xx = bscans(1).wvl;
yy = bscans(1).signal;
f = fit(xx.',yy.',cc);
plot(f,xx,yy)


cscans = scans;

N = length(scans) ;  % Total number of files

for i = 1 : N
    cscans(i).signal(585:614)=14000;
    cscans(i).signal(1050:1090)=1550;
    cscans(i).signal(1525:1580)=1330;
end

figure; for kk=1:10 plot(cscans(kk).wvl, cscans(kk).signal); hold on; end

xx = cscans(1).wvl;
yy = cscans(1).signal;
f = fit(xx.',yy.',cc);
plot(f,xx,yy)


x=59;
y=29;
varA=[]; 
varB=[];
varC=[];

first = 614;
last = 700;

for i = 1 : N
    scans(i).signal(1:first)=0;
    scans(i).signal(last:length(scans(1).signal))=0;
    scans(i).signal(1050:1090)=0;
    scans(i).signal(1525:1580)=0;

    [val,loc] = max(scans(i).signal);
    varA=[varA;val];
    varC=[varC;scans(i).wvl(loc)];
    
%plot(scans(5).wvl,scans(5).signal) 
end

figure; for kk=1:10 plot(scans(kk).wvl, scans(kk).signal); hold on; end

xx = scans(1).wvl;
yy = scans(1).signal;
f = fit(xx.',yy.',cc);
plot(f,xx,yy)


j=1; k=y;    
rowS=[];
rowI=[];
for T= 1 : x
    if mod(T,2)==1
    rowS=[rowS varC(j:k)];
    rowI=[rowI varA(j:k)];
    else
    rowS=[rowS flipud(varC(j:k))];
    rowI=[rowI flipud(varA(j:k))];
    end
    
    k=k+y;
    j=j+y;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imageX = [1:y];
imageY = [1:x];
% colormap(jet)
%%%%%%%
%  imagesc(imageX,imageY,rowS')
%  colorbar
 figure
 IblurS = imgaussfilt(rowS,2);
norma = IblurS;% - max(IblurS(:));
IblurS = norma;% ./ min(norma(:));
colormap(jet)
clims = [545 565]; 
imagesc(imageX,imageY,IblurS', clims)
saveas(gcf,'IXshift','tiff');

% imagesc(imageX,imageY,IblurS');
colorbar
% caxis([0 1])
title('Shift')

savefig('Xsfift.fig');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
%  imagesc(imageX,imageY,rowI')
%  colorbar
 figure
 IblurI = imgaussfilt(rowI,2);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 norma = IblurI - min(IblurI(:));
IblurI = norma ./ max(norma(:));
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap(jet)
 imagesc(imageX,imageY,IblurI')
 saveas(gcf,'IXintensity','tiff');
colorbar
% caxis([0 1])
title('Intensity')

savefig('Xintensity.fig');

Is = imread('IXshift.tif');
Ii = imread('IXintensity.tif');

IsR=(Is(:,:,1));
IsG=(Is(:,:,2));
IsB=(Is(:,:,3));

[xR xP]=imhist(IsR);
[xG xP]=imhist(IsG);
[xB xP]=imhist(IsB);


figure; 
subplot(1,3,1); plot(xP,xR);
subplot(1,3,2); plot(xP,xG);
subplot(1,3,3); plot(xP,xB);

figure; 
subplot(1,3,1); bar(xP,xR,'r');
subplot(1,3,2); bar(xP,xG,'g');
subplot(1,3,3); bar(xP,xB,'b');

entropy(IsR)
entropy(IsG)
entropy(IsB)