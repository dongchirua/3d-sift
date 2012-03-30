function lmstruct = BO_readlm(filename)
%   lmstruct = BO_readlm(filename)
%
% INPUT
%   filename: .lm2 or .lm3 file
%
% OUTPUT
%   lmstruct: structure variable
%       .Location: coordinates of the landmarks
%       .Type: type of the landmark
%
% author: Jeroen Hermans, Dirk Smeets
  lmstruct.Location = [];
  lmstruct.Type = [];
  
  fid = fopen(filename,'r');
  
  % Look for first enumeration = line ending with ':'
  
  stop = 0;  
  while(~stop)    
    line = fgetl(fid);
    if(~isempty(line) & strcmp(line(end),':'))
      stop = 1;
    end
  end
  
  % Read data, separate coordinates from types and enumeration
  % announcements
  
  stop = 0;
  while(~stop)
    line = fgetl(fid);
    
    if(line==-1)
      stop = 1;
    elseif(~isempty(line))
      
      if(~strcmp(line(end),':'))
        % no enumeration announcement
        
        if(isempty(str2num(line)))
          % line = landmark type          
          lmstruct.Type = strvcat(lmstruct.Type,deblank(line));          
        else
          % line = landmark coordinate          
          lmstruct.Location = [lmstruct.Location str2num(line)'];
        end
      
      end
    end
    
  end
  
  fclose(fid);