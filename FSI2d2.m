clc
clear all
close all

addpath(['C:\Users\aj0036\OneDrive\Downloads\Telegram Desktop\139912140']);
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

Lambda_Tresh_L = 539;
Lambda_Tresh_H = 600;

xstd=59;
y=29;

varA=[]; 
varC=[];

varAP=[]; 
varCP=[];

varAPNew=[]; 
varCPNew=[];


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

%     [Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit([xx' yy'],0.5,0,3,6,0,2);
    [Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit([xx' yy'],0.5,0,3,33,0,2);
    [numRowsNew,numColsNew] = size(Results);
    
    [t1,t2]=max(Results(:,5));
        
        valPNew = Results(t2,3);
        ppNew = Results(t2,2);   


    if  Results(t2,2) < Lambda_Tresh_L
        
        valPNew = Results(t2+1,3);
        ppNew = Results(t2+1,2);
    
    end
    if Results(t2,2) > Lambda_Tresh_H
            
        valPNew = Results(t2-1,3);
        ppNew = Results(t2-1,2);
        
    end
    
    varAPNew=[varAPNew;valPNew];
    varCPNew=[varCPNew;ppNew];

    
end


% xx = bscans(1).wvl(first:last);
% yy = LaserRemove3(1,:);
% f = fit(xx.',yy.',cc);
% 
% M = [bscans(kk).wvl(first:last) LaserRemove3(kk,:)];
% 
% figure(22);
% [Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit([xx' yy'],.4,0,2,1,0,10);
% figure(21);
% [Results,LowestError,baseline,BestStart,coeff,xi,yi]=peakfit([xx' yy'],0,0,2,33,0,10);
% figure(20); plot(xi,yi);
% title('Plot of model peaks evaluated at 600 x-values')


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

imageX = [1:y];
imageY = [1:xstd];

figure(10); 
IblurSNew = imgaussfilt(rowS,2);
normaNew = IblurSNew;% - max(IblurS(:));
IblurSNew = normaNew;% ./ min(norma(:));
colormap(jet)
imagesc(imageX,imageY,IblurSNew', clims)
saveas(gcf,'IXshift','tiff');

colorbar
title('Wavelength')
xlabel('Sample Rows');
ylabel('Sample Columns');


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
title('Intensity')
xlabel('Sample Rows');
ylabel('Sample Columns');


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

imageX = [1:y];
imageY = [1:xstd];

figure(20); 
IblurSP = imgaussfilt(rowSP,2);
normaP = IblurSP;% - max(IblurS(:));
IblurSP = normaP;% ./ min(norma(:));
colormap(jet)
imagesc(imageX,imageY,IblurSP', clims)
saveas(gcf,'IXshift','tiff');

colorbar
title('Wavelength')
xlabel('Sample Rows');
ylabel('Sample Columns');

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
title('Intensity')
xlabel('Sample Rows');
ylabel('Sample Columns');


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

imageX = [1:y];
imageY = [1:xstd];

figure(30); 
IblurSNew = imgaussfilt(rowSNew,2);
normaNew = IblurSNew;% - max(IblurS(:));
IblurSNew = normaNew;% ./ min(norma(:));
colormap(jet)
imagesc(imageX,imageY,IblurSNew', clims)
saveas(gcf,'IXshift','tiff');

colorbar
title('Wavelength')
xlabel('Sample Rows');
ylabel('Sample Columns');

savefig('Xsfift.fig');

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
title('Intensity')
xlabel('Sample Rows');
ylabel('Sample Columns');

STDP = std(IblurSNew');
STDP2 = std2(IblurSNew')
MeanP = mean(IblurSNew');
MeanP2 = mean2(IblurSNew')
Error = MeanP-MeanP2;
E = std(IblurSNew').*ones(size(MeanP));

xstd = (1:length(MeanP));
figure(40); errorbar(xstd,MeanP,E);
% ylim([540 570]);
xlabel('Number of Samples');
ylabel('Wavelength');
