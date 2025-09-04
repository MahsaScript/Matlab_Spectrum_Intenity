clc
clear all
close all

addpath('C:\Users\aj0036\OneDrive\Downloads\Telegram Desktop\139912140');
cd 'C:\Users\aj0036\OneDrive\Downloads\Telegram Desktop\139912140'

scans = avantesSpectrumRead;
frequency = scans(1).wvl;

cc = 'gauss4';
clims = [540 570]; 
tune = 7;

bscans = scans;

first = 450;
last = 882;

picshow = 603;
ppf = picshow - first + 1;
freq = bscans(1).wvl(first:last);
pfreq = bscans(1).wvl(picshow:last);


x=59;
y=29;
varA=[]; 
varC=[];

varAP=[]; 
varCP=[];


N = length(scans) ;  % Total number of files

for i = 1 : N
    LaserRemove1(i,:) = medfilt1(bscans(i).signal(first:last),5);
    LaserRemove2(i,:) = filloutliers(LaserRemove1(i,:),'nearest','mean');
    LaserRemove3(i,:) = smoothdata(LaserRemove2(i,:));
    laserf(i,:) = LaserRemove3(i,ppf:length(LaserRemove3(i,:)));
    data(i,:) = bscans(i).signal(first:last);

    [val,loc] = max(laserf(i,:));
    varA=[varA;val];
    varC=[varC;pfreq(loc)];
    
    P = findpeaksG(bscans(i).wvl(first:last),LaserRemove3(i,:),1,1,tune,tune,3);
    [numRows,numCols] = size(P);
    
    valP = P(numRows,3);
    pp = P(numRows,2);
    
    varAP=[varAP;valP];
    varCP=[varCP;pp];
    
end


% figure(1); mesh(data);
% figure(2); mesh(LaserRemove1);
% figure(3); mesh(LaserRemove2);
% figure(4); mesh(LaserRemove3);
% figure(5); mesh(laserf);


figure(5); for kk=1:10 plot(bscans(kk).wvl(first:last), LaserRemove3(kk,:)); hold on; end

xx = bscans(1).wvl(first:last);
yy = LaserRemove3(1,:);
f = fit(xx.',yy.',cc);
figure(6); plot(f,xx,yy);


M = [bscans(kk).wvl(first:last) LaserRemove3(kk,:)];
% % P = ipeak(M,1000)
% % peakfit([M])

figure(21);
[Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit([xx' yy'],0,0,2,33,0,10);
figure(20); plot(xi,yi);
title('Plot of model peaks evaluated at 600 x-values')



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

figure(10); 
IblurS = imgaussfilt(rowS,2);
norma = IblurS;% - max(IblurS(:));
IblurS = norma;% ./ min(norma(:));
colormap(jet)
imagesc(imageX,imageY,IblurS', clims)
saveas(gcf,'IXshift','tiff');

colorbar
title('Shift')

savefig('Xsfift.fig');

figure(11); 
IblurI = imgaussfilt(rowI,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norma = IblurI - min(IblurI(:));
IblurI = norma ./ max(norma(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap(jet)
imagesc(imageX,imageY,IblurI')
saveas(gcf,'IXintensity','tiff');
colorbar
title('Intensity')




j=1; k=y;    
rowSP=[];
rowIP=[];
for T= 1 : x
    if mod(T,2)==1
    rowSP=[rowSP varCP(j:k)];
    rowIP=[rowIP varAP(j:k)];
    else
    rowSP=[rowSP flipud(varCP(j:k))];
    rowIP=[rowIP flipud(varAP(j:k))];
    end
    
    k=k+y;
    j=j+y;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imageX = [1:y];
imageY = [1:x];

figure(100); 
IblurSP = imgaussfilt(rowSP,2);
normaP = IblurSP;% - max(IblurS(:));
IblurSP = normaP;% ./ min(norma(:));
colormap(jet)
imagesc(imageX,imageY,IblurSP', clims)
saveas(gcf,'IXshift','tiff');

colorbar
title('Shift')

savefig('Xsfift.fig');

figure(101); 
IblurIP = imgaussfilt(rowIP,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normaP = IblurIP - min(IblurIP(:));
IblurIP = normaP ./ max(normaP(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap(jet)
imagesc(imageX,imageY,IblurIP')
saveas(gcf,'IXintensity','tiff');
colorbar
title('Intensity')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% cscans = scans;
% 
% x=59;
% y=29;
% varA1=[]; 
% varB1=[];
% varC1=[];
% 
% first1 = 617;
% last1 = 700;
% 
% freq1 = cscans(1).wvl(first1:last1);
% 
% 
% for i = 1 : N
%     cscans(i).signal(1:first1)=0;
%     cscans(i).signal(1050:1090)=0;
%     cscans(i).signal(1525:1580)=0;
%     
%     LaserRemove22(i,:) = filloutliers(cscans(i).signal(first1:last1),'nearest','mean');
%     LaserRemove33(i,:) = smoothdata(LaserRemove22(i,:));
% 
%     [val1,loc1] = max(LaserRemove33(i,:));
%     varA1=[varA1;val1];
%     varC1=[varC1;freq1(loc1)];
%     
% end
% 
% figure(12); for kk=1:10 plot(cscans(kk).wvl(first1:last1), LaserRemove33(kk,:)); hold on; end
% 
% xx1 = cscans(1).wvl(first1:last1);
% yy1 = LaserRemove33(1,:);
% f1 = fit(xx1.',yy1.',cc);
% figure(13); plot(f1,xx1,yy1)
% 
% figure(14); plot(freq1,cscans(i).signal(first1:last1)); hold on; plot(freq1,LaserRemove33(i,:)); hold on; plot(freq1,LaserRemove22(i,:));
% 
% j=1; k=y;    
% rowS1=[];
% rowI1=[];
% for T= 1 : x
%     if mod(T,2)==1
%     rowS1=[rowS1 varC1(j:k)];
%     rowI1=[rowI1 varA1(j:k)];
%     else
%     rowS1=[rowS1 flipud(varC1(j:k))];
%     rowI1=[rowI1 flipud(varA1(j:k))];
%     end
%     
%     k=k+y;
%     j=j+y;
%     
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% imageX1 = [1:y];
% imageY1 = [1:x];
% 
% figure(15); 
% IblurS1 = imgaussfilt(rowS1,2);
% norma1 = IblurS1;% - max(IblurS(:));
% IblurS1 = norma1;% ./ min(norma(:));
% colormap(jet)
% imagesc(imageX1,imageY1,IblurS1', clims)
% saveas(gcf,'IXshift','tiff');
% 
% colorbar
% title('Shift')
% 
% savefig('Xsfift.fig');
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% dscans = scans;
% 
% x=59;
% y=29;
% varA2=[]; 
% varB2=[];
% varC2=[];
% 
% 
% first2 = 450;
% last2 = 882;
% 
% freq2 = dscans(1).wvl(first2:last2);
% 
% laser1 = 585;
% laser2 = 595;
% 
% 
% 
% for i = 1 : N
%     
%     dscans(i).signal(laser1:laser2)=0;
%     LaserRemove111(i,:) = medfilt1(dscans(i).signal(first2:last2),5);
%     LaserRemove222(i,:) = filloutliers(LaserRemove111(i,:),'nearest','mean');
%     LaserRemove333(i,:) = smoothdata(LaserRemove222(i,:));
%     [val2,loc2] = max(LaserRemove333(i,:));
%     varA2=[varA2;val2];
%     varC2=[varC2;freq2(loc2)];
% 
%        
% end
% 
% figure(16); for kk=1:10 plot(dscans(kk).wvl(first2:last2), LaserRemove333(kk,:)); hold on; end
% 
% xx2 = dscans(1).wvl(first2:last2);
% yy2 = LaserRemove333(1,:);
% f2 = fit(xx2.',yy2.',cc);
% figure(17); plot(f2,xx2,yy2)
% 
% 
% figure(18); plot(freq2,dscans(i).signal(first2:last2)); hold on; plot(freq2,LaserRemove333(i,:)); hold on; plot(freq2,LaserRemove222(i,:));
% 
% j=1; k=y;    
% rowS2=[];
% rowI2=[];
% for T= 1 : x
%     if mod(T,2)==1
%     rowS2=[rowS2 varC2(j:k)];
%     rowI2=[rowI2 varA2(j:k)];
%     else
%     rowS2=[rowS2 flipud(varC2(j:k))];
%     rowI2=[rowI2 flipud(varA2(j:k))];
%     end
%     
%     k=k+y;
%     j=j+y;
%     
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% imageX2 = [1:y];
% imageY2 = [1:x];
% 
% figure(19);
% IblurS2 = imgaussfilt(rowS2,2);
% norma2 = IblurS2;% - max(IblurS(:));
% IblurS2 = norma2;% ./ min(norma(:));
% colormap(jet)
% imagesc(imageX2,imageY2,IblurS2', clims)
% saveas(gcf,'IXshift','tiff');
% 
% colorbar
% title('Shift')
% 
% savefig('Xsfift.fig');