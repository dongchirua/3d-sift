function lamdi
a=dir;
for i=3:length(a)
    if length(a(i).name)-2>0
    w=[a(i).name(length(a(i).name)-2) a(i).name(length(a(i).name)-1) a(i).name(length(a(i).name))];
    if w=='cpp'
        mex(a(i).name);
    end
    end
end