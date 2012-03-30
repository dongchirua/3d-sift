function exportallfigs(keepeps)
if nargin < 1
    keepeps = 0;
end
if ~exist('/mic2/cmaes4/Matlab/Thesis/exportedfigs/','dir')
    error 'adapt export folder'
end
stop = 0;
while ~stop
    if gcf == 1
        stop = 1;
    end
    exportfiglaurent(num2str(gcf),gcf,keepeps)
    close(gcf)
end