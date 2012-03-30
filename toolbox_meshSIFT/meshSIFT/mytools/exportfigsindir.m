function exportfigsindir(PATH,keepeps)
if nargin < 2
    keepeps = 0;
    if nargin < 1
        PATH = '.';
    end
end
w = cd;
cd(PATH)
eps = dir('*.eps');
for i = 1:size(eps,1)
    fname = eps(i).name(1:end-4);
    if ~ispc
        v = [cd '/'];
        epsfile = [v fname '.eps'];
        system(['epstopdf ',epsfile]);
        if exist([epsfile(1:end-4) '.pdf'],'file')
            file = dir([epsfile(1:end-4) '.pdf']);
            if file.bytes > 1000
                if ~keepeps
                    delete(epsfile);
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
        v = [cd '\'];
        epsfile = [v fname '.eps'];
        shortpath_eps = ['"' epsfile '"'];
        shortpath_pdf = ['"' v fname '.pdf"'];
        system(['"C:\Program Files\gs\gs8.70\bin\gswin32c.exe" -sDEVICE=pdfwrite -dEPSCrop -o ' shortpath_pdf ' ' shortpath_eps]);
         if exist([epsfile(1:end-4) '.pdf'],'file')
            file = dir([epsfile(1:end-4) '.pdf']);
            if file.bytes > 1000
                if ~keepeps
                    delete(epsfile);
                end
            else
                disp('empty pdf produced... keeping .eps file')
            end
        end
    end
end
cd(w);
end