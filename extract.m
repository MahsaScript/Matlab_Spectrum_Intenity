% Get coordinates of peaks
% The first line below is another option
% [~, locs] = findpeaks(signal, 'MinPeakHeight', range(signal)*.4 + min(signal), 'MinPeakDistance', 200);
load sig.mat
signal=varCPNew;
signalSmooth = smooth(signal,20); 
[~, locs] = findpeaks(signalSmooth, 'MinPeakHeight', range(signal)*.4 + min(signal), 'MinPeakDistance', 200); 

% Look for 4 consecutive increases
mm = movmean([inf;diff(signal)] > 0, [0,3]) == 1;
mmIdx = find(mm); 
firstOnes  = find([mm(1);diff(mmIdx)>1]);
startIdx = mmIdx(firstOnes); 

clf()
plot(signal, 'b-')
hold on
plot(locs, signal(locs), 'r*')
plot(startIdx, signal(startIdx), 'mx', 'LineWidth',2)

% *******
% Now exact the segments and plot them separately
% Extract
segements = arrayfun(@(start,stop){signal(start:stop)},startIdx,locs); 

figure()
hold on
cellfun(@(c)plot(c),segements)

% *****
% Extract
segements = arrayfun(@(start,stop){signal(start:stop)},startIdx,locs); 
segement_xvalues = arrayfun(@(start,stop){start:stop},startIdx,locs); 

figure()
hold on
cellfun(@(x,y)plot(x,y),segement_xvalues, segements)

