function scans = avantesSpectrumRead(dataFiles)
% Imports scan data from binary or text files saved by Avantes software
%
% DATAFILES is either
% - nonexistent, then roh, trt, tat and ttt (in this order) files in the
%   current directory are used or
% - a list of data files to read as if returned by the 'dir' command or
% - a type of data files to use, 'trt', 'tat', or 'ttt'
%
% Copyright: (cc-by) Kotya Karapetyan, 2011.
% kotya.karapetyan@gmail.com
% http://agmeschede.iap.uni-bonn.de


scans = [];

timer = tic;

%% Create a random name for temporary file
tempfilename = ['scans_' dec2hex(randi(1e9,1)) '.dat'];
while exist(tempfilename, 'file')
    tempfilename = ['scans_' dec2hex(randi(1e9,1)) '.dat'];
end

%% Process input argument
if nargin == 0
    dataFiles = dir;
    if strcmpi(dataFiles(2).name, '..')
        dataFiles(2) = [];
    end
    if strcmpi(dataFiles(1).name, '.')
        dataFiles(1) = [];
    end
elseif iscell(dataFiles) % list of strings provided
    filenames = dataFiles;
    dataFiles = [];
    for i = 1:numel(filenames)
        dataFiles = [dataFiles dir(filenames{i})]; %#ok<AGROW>
    end
elseif ischar(dataFiles) % string provided
    filelist = dir(dataFiles);
    if isempty(filelist)
        addtolog(['Input argument is a string (' dataFiles ...
            ') but no matching files found. May be a case problem']);
        fprintf('File missing. Return. May be a case problem\n')
        savelogfile;
        return
    end
    dataFiles = filelist;
    clear filelist
elseif strcmpi(dataFiles, 'trt') || strcmpi(dataFiles, 'tat') ...
        || strcmpi(dataFiles, 'ttt') % file type is given
    dataFileExtension = dataFiles;
    dataFiles = dir(['*.' dataFileExtension]);
else
    error('Something wrong with input arguments')
end

%% Read data files
scans = struct('sample', '');
scans(numel(dataFiles)) = scans;

i = 0;
for k=1:numel(dataFiles)
    filename = dataFiles(k).name;
    if dataFiles(k).isdir % skip subdirectories
        continue
    end
    fprintf('%s. %g files to go\n', filename, numel(dataFiles)-k);
    
    [~, ~, dataFileExtension] = fileparts(filename);
    dataFileExtension(1) = []; % remove leading dot
    switch lower(dataFileExtension)
        case {'trt', 'roh'}
            commentFileExtension = 'RCM';
            %         column = 2;
        case {'tat', 'abs'}
            commentFileExtension = 'ACM';
            %         column = 5;
        case {'ttt', 'trm'}
            commentFileExtension = 'TCM';
            %         column = 5;
        otherwise
            continue;
    end
    i = i + 1;
    
    filetime = dataFiles(k).datenum;
    if strcmpi(dataFileExtension, 'ROH') ...
            || strcmpi(dataFileExtension, 'TRM') ...
            || strcmpi(dataFileExtension, 'ABS')
        s = avantesBinaryRead(filename);
        scans(i).serialNumber = s.serialnr;
        scans(i).integrationTime = s.inttime;
        scans(i).avg = s.average;
        scans(i).smoothingPixels = s.smoothpix;
        scans(i).spectrometerName = s.userfriendlyname;
        scans(i).avantesTimestamp = s.timestamp;
        scans(i).smpl = numel(s.scope);

        scans(i).signal = s.scope'; 
        scans(i).reference = s.reference';
        scans(i).dark = s.dark';
        
        scans(i).wvl = zeros(1, scans(i).smpl);
        j = 0:numel(scans(i).signal)-1;
        scans(i).wvl = s.WLIntercept + s.WLX1 * j + s.WLX2 * j.^2 + ...
            s.WLX3 * j.^3 + s.WLX4 * j.^4;
        assert(numel(scans(i).signal) == numel(scans(i).wvl));
        
        scans(i).timeStart = filetime;
        scans(i).timeEnd = scans(i).timeStart;
        
        commentFileName = [filename(1:end-4) '.' commentFileExtension];
        if exist(commentFileName, 'file')
            f = fopen(commentFileName);
            text = textscan(f, '%s');
            text = text{1};
            comment = [];
            for j = 1:numel(text)
                comment = [comment text{j} ' ']; %#ok<AGROW>
            end
            comment(end) = []; % remove trailing space
            assert(strcmpi(comment(1:9), scans(i).serialNumber));
            if numel(comment) < 12 % no comment present
                scans(i).comment = [];
            else
                assert(strcmpi(comment(10:11), '- '));
                scans(i).comment = comment(12:end);
                fclose(f);
                scans(i).timeStart = filetime;
            end
        else
            scans(i).timeStart = [];
            addtolog(['Comment file ' commentFileName ' not found']);
        end
        
    else % not ROH, ABS or TRM files, assume a T*T (ascii) file
        fid = fopen(filename);
        s = textscan(fid, '%s', 'delimiter', '\n'); % skip headerLinesNumber lines, read the rest into string cells
        assert(feof(fid) == true);
        fclose(fid);
        s = s{1};
        
        % Read header
        % 1st line: serial number and comment
        textline = s{1};
        scans(i).serialNumber = textline(1:9);
        scans(i).comment = textline(12:end);
        
        % 2nd line: integration time
        textline = s{2};
        textline = textline(19:end);
        spacePos = strfind(textline, ' ');
        textline(spacePos(1):end) = []; % remove 'ms' at the end
        scans(i).integrationTime = str2double(textline);
        
        % 3rd line: averaging
        textline = s{3};
        textline = textline(10:end);
        spacePos = strfind(textline, ' ');
        textline(spacePos(1):end) = []; % remove 'scans' at the end
        scans(i).avg = str2double(textline);
        
        % 4th line: smoothing
        textline = s{4};
        textline = textline(34:end);
        scans(i).smoothingPixels = str2double(textline);
        
        % 5th line: spectrometer name (always serial number?)
        textline = s{5};
        scans(i).spectrometerName = textline(39:end);
        
        % 6th line: timestamp
        textline = s{6};
        scans(i).avantesTimestamp = str2double(textline(30:end));
        commentFileName = [filename(1:end-4) '.' commentFileExtension];
        if exist(commentFileName, 'file')
            filemeta = dir(commentFileName);
            timestamp = filemeta.datenum;
            scans(i).timeStart = timestamp;
        else
            scans(i).timeStart = [];
            addtolog(['Comment file ' commentFileName ' not found']);
        end
        scans(i).timeEnd = scans(i).timeStart;
        
        headerLinesNumber = 8;
        
        % Other data
        
        % Read data
        wvl = NaN * ones(1, numel(s) - headerLinesNumber);
        signal = wvl;
        for j = 1:(numel(s) - headerLinesNumber)
            textline = s{j + headerLinesNumber};
            textline(textline == ',') = '.'; % replace all decimal commas, if any, with points
            data = regexp(textline, ';', 'split');
            try
                wvl(j) = str2double(data{1});
                signal(j) = str2double(data{2});
            catch ME %#ok<NASGU>
                fprintf('Cannot read wavelength and signal, check if the file header contains %n lines\n', headerLinesNumber)
            end
        end
        
        scans(i).wvl = wvl;
        scans(i).signal = signal;
        scans(i).smpl = numel(wvl);
    end
    scans(i).filename = filename;
    scans(i).sample = avantesBinaryRead(filename);
    scans(i).resln = 'Depends on installed slit or fibre core diameter';
    
    saveeach = 300;
    if mod(i, saveeach) == 0
        try
            save(tempfilename, 'scans', '-mat')
        catch ME  %#ok<NASGU>
            fprintf('Saving temporary file failed. This may be due to the lack of disk space\n')
            fprintf('The scans variable now contains %g spectra. BE CAREFUL: scans have not been saved\n', numel(scans))
            if exist(tempfilename, 'file')
                delete(tempfilename) % delete the incomplete temp file
                fprintf('Temporary file %s deleted\n', tempfilename);
            end
            return
        end
        fprintf('Saved\n'); % confirm saving to temp file
    end
end

scans(i+1:end) = [];

%% If many scans read, save to a file with the directory name
if numel(scans) > 50
    [~, scansFileName] = fileparts(cd);
    scansFileName = [scansFileName '.dat'];
    [~, IX] = sort([scans(:).timeStart]); % sort scans by start time
    scans = scans(IX);
    save(scansFileName, 'scans', '-mat'); % save the final file
    fprintf('%g spectra processed\n', numel(scans))
    fprintf('Spectra saved to %s\n', scansFileName);
    if exist(tempfilename, 'file')
        delete(tempfilename) % delete the temporary file
        fprintf('Temporary file %s deleted\n', tempfilename);
    end
end

% savelogfile;
fprintf('Time: %g s\n', ceil(toc(timer)))


%% List of comment files
% commentFiles = [];
% for i = 1:numel(dataFiles)
%
%         commentFiles(j) = dir([int2str(includeExperiments) '*' zero int2str(includeSpectra(1,j)) commentFileExtension]); % get list of experiment comment files in directory to analyse
%     end
% end;
