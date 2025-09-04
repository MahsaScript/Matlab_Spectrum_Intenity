function s = avantesBinaryRead(filename)
% Reads a single binary file saved by Avantes AvaSoft7USB2 software. 

% Copyright: (cc-by) Kotya Karapetyan, 2011.
% kotya.karapetyan@gmail.com
% http://agmeschede.iap.uni-bonn.de

s = [];

if exist(filename, 'file')
    f = fopen(filename);
else
    warning('FileError:NoFile', 'File %s does not exist\n', filename);
    return
end

%% Header
s.versionID = fread(f, 1, 'single');
s.serialnr = char(fread(f, 9, 'single')); %#ok<FREAD>
s.serialnr = s.serialnr';
s.userfriendlyname = char(fread(f, 64, 'single'));  %#ok<FREAD>
s.userfriendlyname = s.userfriendlyname';

s.WLIntercept = fread(f, 1, 'single');
s.WLX1 = fread(f, 1, 'single');
s.WLX2 = fread(f, 1, 'single');
s.WLX3 = fread(f, 1, 'single');
s.WLX4 = fread(f, 1, 'single');
s.ipixfirst = fread(f, 1, 'single');
s.ipixlast = fread(f, 1, 'single');
s.measuremode = fread(f, 1, 'single');
s.dummy1 = fread(f, 1, 'single');
s.laserwavelength = fread(f, 1, 'single');
s.laserdelay = fread(f, 1, 'single');
s.laserwidth = fread(f, 1, 'single');
s.strobercontrol = fread(f, 1, 'single');
s.dummy2 = fread(f, 1, 'single');
s.dummy3 = fread(f, 1, 'single');
s.timestamp = fread(f, 1, 'single');
s.dyndarkcorrection = fread(f, 1, 'single');

s.smoothpix = fread(f, 1, 'single'); %			//see measurement structure
s.smoothmodel = fread(f, 1, 'single'); %		//see measurement structure
s.triggermode = fread(f, 1, 'single'); %		//see measurement structure
s.triggersource = fread(f, 1, 'single'); %		//see measurement structure
s.triggersourcetype = fread(f, 1, 'single'); %	//see measurement structure
s.NTC1 = fread(f, 1, 'single'); %			//onboard temp in degrees Celsius
s.NTC2 = fread(f, 1, 'single'); %			//NTC2 in Volt (not connected)
s.Thermistor = fread(f, 1, 'single'); %		//detector temp in degr Celsius (only TEC, NIR)
s.dummy4 = fread(f, 1, 'single'); %			//

%% Data
[~, ~, ext] = fileparts(filename);
if strcmpi(ext, '.ABS') || strcmpi(ext, '.TRM')
    data = fread(f, (s.ipixlast - s.ipixfirst + 1) * 3, 'single');
    s.scope = data(1:3:end);
    s.reference = data(2:3:end);
    s.dark = data(3:3:end);
elseif strcmpi(ext, '.ROH') % scope mode
    s.scope = fread(f, s.ipixlast - s.ipixfirst + 1, 'single');
    s.reference = NaN;
    s.dark = NaN;
else
    error('Unspported file'); 
end

%% Footer
s.inttime = fread(f, 1, 'single'); %			//integration time [ms] during file-save
s.average = fread(f, 1, 'single'); %			//nr of average during file-save
s.integrationdelay = fread(f, 1, 'single'); %	//see measurement structure

fclose(f);

