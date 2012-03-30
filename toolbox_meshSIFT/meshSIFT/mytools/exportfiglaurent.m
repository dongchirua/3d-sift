function exportfig(fname, f,keepeps, PATH)

%%%%%%%%%%%%%%%%%
% SAVE LOCATION %
%%%%%%%%%%%%%%%%%


% Define full path.
if nargin < 4
    PATH = ['/mic2/cmaes4/Matlab/Thesis/exportedfigs/'];
    if nargin < 3
        keepeps = 0;
        if nargin < 2
            f = gcf;
        end
    end
end
if nargin < 1
	fname = get(f, 'Name');
end
shortpath = [PATH, fname];

%%%%%%%%%%%%%%
% PREPROCESS %
%%%%%%%%%%%%%%
saveas(f,shortpath,'fig');

% Define objects.
allLines  = findall(f, 'type', 'line');
allText   = findall(f, 'type', 'text');
allAxes   = findall(f, 'type', 'axes');
allImages = findall(f, 'type', 'image');
allLights = findall(f, 'type', 'light');
allPatch  = findall(f, 'type', 'patch');
allSurf   = findall(f, 'type', 'surface');
allRect   = findall(f, 'type', 'rectangle');
allFont   = [allText; allAxes];
allColor  = [allLines; allText; allAxes; allLights];
allMarker = [allLines; allPatch; allSurf];
allEdge   = [allPatch; allSurf];
allCData  = [allImages; allPatch; allSurf];
allLWidth = [allLines; allAxes; allSurf; allPatch];

% Scale lines.
scale = 3;
oldlines = LocalGetAsCell(allLWidth,'LineWidth');
newlines = LocalScale(oldlines, scale);
set(allLWidth,{'LineWidth'},newlines);

% Scale fontsize.
oldaxes = LocalGetAsCell(allAxes,'FontSize');
set(allAxes,{'FontSize'},{11});
oldtext = LocalGetAsCell(allText,'FontSize');
set(allText,{'FontSize'},{13});

%%%%%%%%%
% PRINT %
%%%%%%%%%

shortpath_eps = [shortpath, '.eps'];
print(f,'-depsc','-r600',shortpath_eps);
if ~ispc
    system(['epstopdf ',shortpath_eps]);
    if exist([shortpath_eps(1:end-4) '.pdf'],'file')
        file = dir([shortpath_eps(1:end-4) '.pdf']);
        if file.bytes > 1200
            if ~keepeps
                delete(shortpath_eps);
            end
        else
            disp('empty pdf produced... keeping .eps file')
        end
    end
else
    % TODO vul opties in zoals voorbeeld in comments hieronder.
    % -dEmbedAllFonts=false -dUseFlateCompression=true
    % -dAutoRotatePages=/None -dHaveTrueTypes -r300-dSubsetFonts=true
    % -dNOPLATFONTS -dUseCIEColor=true
    % -dColorConversionStrategy=/UseDeviceIndependentColor
    % -dProcessColorModel=/DeviceRGB -dAntiAliasColorImages=false
    % -sDEVICE=pdfwrite
    % -sOutputFile="output.pdf" -dUseFlateCompression=true -dLZWEncodePages=true
    % -dCompatibilityLevel=1.6 -dAutoFilterColorImages=false
    % -dAutoFilterGrayImages=false  -dColorImageFilter=/FlateEncode
    % -dGrayImageFilter=/FlateEncode -f "output.eps"
    system(['gswin32c.exe -opties']);
end

%%%%%%%%%
% RESET %
%%%%%%%%%

set(allText,{'FontSize'},oldtext);
set(allAxes,{'FontSize'},oldaxes);
set(allLWidth,{'LineWidth'},oldlines);
    
end

function cellArray = LocalGetAsCell(fig,prop,allowemptycell)
    cellArray = get(fig,prop);
    if nargin < 3
      allowemptycell = 0;
    end
    if ~iscell(cellArray) && (allowemptycell || ~isempty(cellArray))
      cellArray = {cellArray};
    end
end

function newArray = LocalScale(inArray, scale)
    n = length(inArray);
    newArray = cell(n,1);
    for k=1:n
      newArray{k} = scale*inArray{k}(1);
    end
end