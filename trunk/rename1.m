function [ output_args ] = rename1( x,i,l )
j=1;
name=num2str(i);
list=dir('colordata');
for k=l:(l+x-1)
    movefile(strcat('colordata/',list(l).name),strcat('colordata/',name,'_',num2str(j),'.jpg'));
    movefile(strcat('bindata/',list(l).name),strcat('bindata/',name,'_',num2str(j),'.jpg'));
    j=j+1;
    l=l+1;
end
i=i+1;
save i;
save l;
end

