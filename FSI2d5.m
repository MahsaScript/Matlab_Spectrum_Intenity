clc
clear all
close all

% addpath('C:\Users\aj0036\OneDrive\Downloads\Telegram Desktop\139912140');
% cd 'C:\Users\aj0036\OneDrive\Downloads\Telegram Desktop\139912140'
addpath('C:\Users\mofid\Downloads\1-Code2\Code\139912140');
cd 'C:\Users\mofid\Downloads\1-Code2\Code\139912140'
scans = avantesSpectrumRead;
frequency = scans(1).wvl;

cc = 'gauss4';
clims = [540 570]; 
tune = 7;

bscans = scans;

first = 450;
last = 882;

Lambda_Tresh = 540;     % Threshold
g1 = 4;                   % Number of Shapes
h1 = 42;

g2 = 2;                   % Number of Shapes
h2 = 7;

Cof = 0.6; % 0.6 for SCC and Health, 0.85 for BCC and Melanoma

picshow = 603;
ppf = picshow - first + 1;
freq = bscans(1).wvl(first:last);
pfreq = bscans(1).wvl(picshow:last);

xstd=59;
y=29;

varAPNew=[]; 
varCPNew=[];

varAPNew_No=[]; 
varCPNew_No=[];

varCCNew=[];
varCCNew_No=[];


N = length(scans) ;  % Total number of files

% xstd=5;
% y=4;
% N = xstd*y;

for i = 1 : N

    LaserRemove1(i,:) = medfilt1(bscans(i).signal(first:last),5);
    LaserRemove2(i,:) = filloutliers(LaserRemove1(i,:),'nearest','mean');
    LaserRemove3(i,:) = smoothdata(LaserRemove2(i,:));
  
    xx = bscans(i).wvl(first:last);
    yy = LaserRemove2(i,:);

    yy_No = LaserRemove1(i,:);
    yy_No(131:170) = 0; 
    yy_No(131:170) = max(yy_No(171:210))*Cof;
    yyc = yy_No;

    [Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit([xx' yy'],0.5,0,g1,h1,0,2);
    [numRowsNew,numColsNew] = size(Results);
    
    ResultsN = Results;
    
    for h = 1:g1
        if ResultsN(h,2) < Lambda_Tresh
            ResultsN(h,5) = 0;
        end
        if ResultsN(h,4) > 100
            ResultsN(h,5) = 0;
        end
        if ResultsN(h,2) > 600
            ResultsN(h,5) = 0;
        end
    end

    [t1,t2]=max(ResultsN(:,5));
        
    valPNew = ResultsN(t2,3);
    ppNew = ResultsN(t2,2);     
    
    varAPNew=[varAPNew;valPNew];
    varCPNew=[varCPNew;ppNew];


    [Results_No,LowestError_No,baseline_No,BestStart_No,coeff_No,xi_No,yi_No]=peakfit([xx' yyc'],0.5,0,g2,h2,0,2);
    [numRowsNew_No,numColsNew_No] = size(Results_No);


    for h = 1:g2
        if Results_No(h,2) < Lambda_Tresh
            Results_No(h,5) = 0;
        end
        if Results_No(h,4) > 100
            Results_No(h,5) = 0;
        end
        if Results_No(h,2) > 600
            Results_No(h,5) = 0;
        end
    end

    [t1,t2]=max(Results_No(:,5));
        
    valPNew_No = Results_No(t2,3);
    ppNew_No = Results_No(t2,2);     
    
    varAPNew_No=[varAPNew_No;valPNew_No];
    varCPNew_No=[varCPNew_No;ppNew_No];

    XX = ppNew;
    if XX>560
        XX = 560;
    end
    coefficients = [9.23e-245 1.0028 0 0.000625 0.34]; 
    ccNew = coefficients(1) * exp(coefficients(2)*XX) + coefficients(3) + coefficients(4)*XX + coefficients(5);
    varCCNew=[varCCNew;ccNew];


    XX = ppNew_No;
    if XX>560
        XX = 560;
    end
    coefficients = [9.23e-245 1.0028 0 0.000625 0.34]; 
    ccNew_No = coefficients(1) * exp(coefficients(2)*XX) + coefficients(3) + coefficients(4)*XX + coefficients(5);
    varCCNew_No=[varCCNew_No;ccNew_No];

    
end

close all

yiEstim = yi(1,:)+yi(2,:)+yi(3,:)+yi(4,:); 
freq = xx';
yiMeas = yy';
save spectrum.mat freq yiMeas Results LowestError baseline BestStart coeff xi yi yiEstim

yiEstim_No = yi_No(1,:)+yi_No(2,:); 
freq_No = xx';
yiMeas_No = yyc';
save spectrum_No.mat freq_No yiMeas_No Results_No LowestError_No ...
    baseline_No BestStart_No coeff_No xi_No yi_No yiEstim_No

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pixelsize = 16.67;

imageX = [1:y]*Pixelsize;
imageY = [1:xstd]*Pixelsize;

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
%%%%%%%%    Forth Algorithm

j=1; k=y;    
rowSNew_No=[];
rowINew_No=[];
for T= 1 : xstd
    if mod(T,2)==1
    rowSNew_No=[rowSNew_No varCPNew_No(j:k)];
    rowINew_No=[rowINew_No varAPNew_No(j:k)];
    else
    rowSNew_No=[rowSNew_No flipud(varCPNew_No(j:k))];
    rowINew_No=[rowINew_No flipud(varAPNew_No(j:k))];
    end
    
    k=k+y;
    j=j+y;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(50); 
IblurSNew_No = imgaussfilt(rowSNew_No,2);
normaNew_No= IblurSNew_No;% - max(IblurS(:));
IblurSNew_No = normaNew_No;% ./ min(norma(:));
colormap(jet)
imagesc(imageX,imageY,IblurSNew_No', clims)
saveas(gcf,'IXshift','tiff');

colorbar
title('Fluorescence Wavelength Shift')
xlabel('X [\mu m]');
ylabel('Y [\mu m]');

savefig('Xsfift.fig');

SpectrumWImage_No = IblurSNew_No';

save WImage_No.mat imageX imageY SpectrumWImage
% C = load('WImage.mat');


figure(51); 
IblurINew_No = imgaussfilt(rowINew_No,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normaNew_No = IblurINew_No - min(IblurINew_No(:));
IblurINew_No = normaNew_No ./ max(normaNew_No(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap(jet)
imagesc(imageX,imageY,IblurINew_No')
saveas(gcf,'IXintensity','tiff');
colorbar
title('Fluorescence Intensity')
xlabel('X [\mu m]');
ylabel('Y [\mu m]');

SpectrumIImage_No = IblurINew_No';

save IImage_No.mat imageX imageY SpectrumIImage
% D = load('IImage.mat');


STDP_No = std(IblurSNew_No'); 
STDP2_No = std2(IblurSNew_No')
MeanP_No = mean(IblurSNew_No'); 
MeanP2_No = mean2(IblurSNew_No')
Error_No = MeanP_No-MeanP2_No; 
E_No = std(IblurSNew_No').*ones(size(MeanP_No)); 

xxstd_No = (1:length(MeanP_No))*Pixelsize;
figure(52); errorbar(xxstd_No,MeanP_No,E_No);
% ylim([540 570]);
xlabel('\mu m');
ylabel('Wavelength [nm]');

save statist_No.mat STDP_No STDP2_No MeanP_No MeanP2_No Error_No E_No xxstd_No
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

climsC = [1500 3000]; 

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j=1; k=y;    
rowCNew_No=[];
for T= 1 : xstd
    if mod(T,2)==1
    rowCNew_No=[rowCNew_No varCCNew_No(j:k)];
    else
    rowCNew_No=[rowCNew_No flipud(varCCNew_No(j:k))];
    end
    
    k=k+y;
    j=j+y;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

climsC = [1500 3000]; 

figure(101); 
IblurCNew_No = imgaussfilt(rowCNew_No,2);
normaNew_No = IblurCNew_No;% - max(IblurS(:));
IblurCNew_No = normaNew_No;% ./ min(norma(:));
colormap(jet)
imagesc(imageX,imageY,IblurCNew_No*2000', climsC)

colorbar
title('Fluorescence Concentration')
xlabel('X [\mu m]');
ylabel('Y [\mu m]');


SpectrumCImage_No = IblurCNew_No*2000;

save CImage_No.mat imageX imageY SpectrumCImage_No
% E = load('CImage.mat');


[Results_No,LowestError_No,baseline_No,BestStart_No,coeff_No,xi_No,yi_No]=peakfit([xx' yyc'],0.5,0,g2,h2,0,2);
[numRowsNew_No,numColsNew_No] = size(Results_No);
% hold on; plot(xx,LaserRemove1(i,:));
% xlabel('Wavelength [nm]');
% ylabel('Fluorescence Intensity');
% ylim([0 2.5e4]);