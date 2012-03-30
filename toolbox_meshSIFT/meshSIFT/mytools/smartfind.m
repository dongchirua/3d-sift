function i = smartfind(filename,poses)
% this function returns the index of the pose from the poses-cell 
% that is present in the string 'filename'
% returns zero when there is no pose present
% author: Chris Maes
% 2009/10

len = length(poses);
i = len;
while i > 0
    if findstr(filename,poses{i})
        return
    else
        i = i-1;
    end
end
end