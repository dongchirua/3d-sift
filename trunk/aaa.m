li=dir('colordata');
li1=dir('bindata');
for a=3:442
    flat=0;
    for b=3:441
        if(strcmp(li(a).name,li1(b).name)==1)
        flat=1;
        end
    end
    if(flat==0) g=li(a).name;
    end
end
