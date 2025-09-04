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

Lambda_Tresh = 538;

xstd=59;
y=29;

varA=[]; 
varC=[];

varAP=[]; 
varCP=[];

varAPNew=[]; 
varCPNew=[];

varCCNew=[];


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
    
    xx = bscans(i).wvl(first:last);
    yy = LaserRemove3(i,:);

    [Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit([xx' yy'],0.5,0,5,42,0,2);
%     [Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit([xx' yy'],0.5,0,5,42,0,2);
    [numRowsNew,numColsNew] = size(Results);
    
    [t1,t2]=max(Results(:,5));
        
        valPNew = Results(t2,3);
        ppNew = Results(t2,2);     
    
    varAPNew=[varAPNew;valPNew];
    varCPNew=[varCPNew;ppNew];

    XX = ppNew;
    if XX>560
        XX = 560;
    end
    coefficients = [9.23e-245 1.0028 0 0.000625 0.34]; 
    ccNew = coefficients(1) * exp(coefficients(2)*XX) + coefficients(3) + coefficients(4)*XX + coefficients(5);
    varCCNew=[varCCNew;ccNew];

    
end


xx = bscans(1).wvl(first:last);
yy = LaserRemove3(1,:);
% f = fit(xx.',yy.',cc);
kk = length(LaserRemove3);
M = [bscans(kk).wvl(first:last) LaserRemove3(kk,:)];

figure(14);
[Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit([xx' yy'],0.5,0,3,6,0,2);

yiEstim = yi(1,:)+yi(2,:)+yi(3,:); 
freq = xx';
yiMeas = yy';
save spectrum14.mat freq yiMeas Results LowestError baseline BestStart coeff xi yi yiEstim


figure(15);
[Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit([xx' yy'],0.5,0,3,33,0,2);

yiEstim = yi(1,:)+yi(2,:)+yi(3,:); 
freq = xx';
yiMeas = yy';
save spectrum15.mat freq yiMeas Results LowestError baseline BestStart coeff xi yi yiEstim


figure(16);
[Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit([xx' yy'],0.5,0,3,32,0,2);

yiEstim = yi(1,:)+yi(2,:)+yi(3,:); 
freq = xx';
yiMeas = yy';
save spectrum16.mat freq yiMeas Results LowestError baseline BestStart coeff xi yi yiEstim


figure(17); 
% plot(xi,yi); title('Plot of model peaks evaluated at 600 x-values')
[Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit([xx' yy'],0.5,0,3,13,0,2);

yiEstim = yi(1,:)+yi(2,:)+yi(3,:); 
freq = xx';
yiMeas = yy';
save spectrum17.mat freq yiMeas Results LowestError baseline BestStart coeff xi yi yiEstim


figure(18); plot(xi,yi);
[Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit([xx' yy'],0.5,0,3,3,0,2);

yiEstim = yi(1,:)+yi(2,:)+yi(3,:); 
freq = xx';
yiMeas = yy';
save spectrum18.mat freq yiMeas Results LowestError baseline BestStart coeff xi yi yiEstim


figure(19); plot(xi,yi);
[Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit([xx' yy'],0.5,0,3,2,0,2);

yiEstim = yi(1,:)+yi(2,:)+yi(3,:); 
freq = xx';
yiMeas = yy';
save spectrum19.mat freq yiMeas Results LowestError baseline BestStart coeff xi yi yiEstim


figure(97); 
% plot(xi,yi); title('Plot of model peaks evaluated at 600 x-values')
[Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit([xx' yy'],0.5,0,3,1,0,2);

yiEstim = yi(1,:)+yi(2,:)+yi(3,:); 
freq = xx';
yiMeas = yy';
save spectrum97.mat freq yiMeas Results LowestError baseline BestStart coeff xi yi yiEstim


figure(98); plot(xi,yi);
[Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit([xx' yy'],0.5,0,3,7,0,2);

yiEstim = yi(1,:)+yi(2,:)+yi(3,:); 
freq = xx';
yiMeas = yy';
save spectrum98.mat freq yiMeas Results LowestError baseline BestStart coeff xi yi yiEstim


figure(99); plot(xi,yi);
[Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit([xx' yy'],0.5,0,3,42,0,2);

yiEstim = yi(1,:)+yi(2,:)+yi(3,:); 
freq = xx';
yiMeas = yy';
save spectrum99.mat freq yiMeas Results LowestError baseline BestStart coeff xi yi yiEstim
% B = load('spectrum99.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%          First Algorithm


j=1; k=y;    
rowS=[];
rowINew=[];
for T= 1 : xstd
    if mod(T,2)==1
    rowS=[rowS varC(j:k)];
    rowINew=[rowINew varA(j:k)];
    else
    rowS=[rowS flipud(varC(j:k))];
    rowINew=[rowINew flipud(varA(j:k))];
    end
    
    k=k+y;
    j=j+y;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pixelsize = 16.67;



imageX = [1:y]*Pixelsize;
imageY = [1:xstd]*Pixelsize;

figure(10); 
IblurSNew = imgaussfilt(rowS,2);
normaNew = IblurSNew;% - max(IblurS(:));
IblurSNew = normaNew;% ./ min(norma(:));
colormap(jet)
imagesc(imageX,imageY,IblurSNew', clims)
saveas(gcf,'IXshift','tiff');

colorbar
title('Fluorescence Wavelength Shift')
xlabel('X [\mu m]');
ylabel('Y [\mu m]');


savefig('Xsfift.fig');

figure(11); 
IblurINew = imgaussfilt(rowINew,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normaNew = IblurINew - min(IblurINew(:));
IblurINew = normaNew ./ max(normaNew(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap(jet)
imagesc(imageX,imageY,IblurINew')
saveas(gcf,'IXintensity','tiff');
colorbar
title('Fluorescence Intensity')
xlabel('X [\mu m]');
ylabel('Y [\mu m]');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%          Second Algorithm



j=1; k=y;    
rowSP=[];
rowIP=[];
for T= 1 : xstd
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

figure(20); 
IblurSP = imgaussfilt(rowSP,2);
normaP = IblurSP;% - max(IblurS(:));
IblurSP = normaP;% ./ min(norma(:));
colormap(jet)
imagesc(imageX,imageY,IblurSP', clims)
saveas(gcf,'IXshift','tiff');

colorbar
title('Fluorescence Wavelength Shift')
xlabel('X [\mu m]');
ylabel('Y [\mu m]');

savefig('Xsfift.fig');

figure(21); 
IblurIP = imgaussfilt(rowIP,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normaP = IblurIP - min(IblurIP(:));
IblurIP = normaP ./ max(normaP(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap(jet)
imagesc(imageX,imageY,IblurIP')
saveas(gcf,'IXintensity','tiff');
colorbar
title('Fluorescence Intensity')
xlabel('X [\mu m]');
ylabel('Y [\mu m]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    Third Algorithm

j=1; k=y;    
rowSNew=[];
rowINew=[];
for T= 1 : xstd
    if mod(T,2)==1
    rowSNew=[rowSNew varCPNew(j:k)];
    rowINew=[rowINew varAPNew(j:k)];
    else
    rowSNew=[rowSNew flipud(varCPNew(j:k))];
    rowINew=[rowINew flipud(varAPNew(j:k))];
    end
    
    k=k+y;
    j=j+y;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(30); 
IblurSNew = imgaussfilt(rowSNew,2);
normaNew = IblurSNew;% - max(IblurS(:));
IblurSNew = normaNew;% ./ min(norma(:));
colormap(jet)
imagesc(imageX,imageY,IblurSNew', clims)
saveas(gcf,'IXshift','tiff');

colorbar
title('Fluorescence Wavelength Shift')
xlabel('X [\mu m]');
ylabel('Y [\mu m]');

savefig('Xsfift.fig');

SpectrumWImage = IblurSNew';

save WImage.mat imageX imageY SpectrumWImage
% C = load('WImage.mat');


figure(31); 
IblurINew = imgaussfilt(rowINew,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normaNew = IblurINew - min(IblurINew(:));
IblurINew = normaNew ./ max(normaNew(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap(jet)
imagesc(imageX,imageY,IblurINew')
saveas(gcf,'IXintensity','tiff');
colorbar
title('Fluorescence Intensity')
xlabel('X [\mu m]');
ylabel('Y [\mu m]');

SpectrumIImage = IblurINew';

save IImage.mat imageX imageY SpectrumIImage
% D = load('IImage.mat');


STDP = std(IblurSNew'); 
STDP2 = std2(IblurSNew')
MeanP = mean(IblurSNew'); 
MeanP2 = mean2(IblurSNew')
Error = MeanP-MeanP2; 
E = std(IblurSNew').*ones(size(MeanP)); 

xxstd = (1:length(MeanP))*Pixelsize;
figure(40); errorbar(xxstd,MeanP,E);
% ylim([540 570]);
xlabel('\mu m');
ylabel('Wavelength [nm]');

save statist.mat STDP STDP2 MeanP MeanP2 Error E xxstd
% A = load('statist.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    Concentration Algorithm

j=1; k=y;    
rowCNew=[];
for T= 1 : xstd
    if mod(T,2)==1
    rowCNew=[rowCNew varCCNew(j:k)];
    else
    rowCNew=[rowCNew flipud(varCCNew(j:k))];
    end
    
    k=k+y;
    j=j+y;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

climsC = [0 2000]; 

figure(100); 
IblurCNew = imgaussfilt(rowCNew,2);
normaNew = IblurCNew;% - max(IblurS(:));
IblurCNew = normaNew;% ./ min(norma(:));
colormap(jet)
imagesc(imageX,imageY,IblurCNew*2000', climsC)

colorbar
title('Fluorescence Concentration')
xlabel('X [\mu m]');
ylabel('Y [\mu m]');


SpectrumCImage = IblurCNew*2000;

save CImage.mat imageX imageY SpectrumCImage
% E = load('CImage.mat');
