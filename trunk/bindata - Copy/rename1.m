function [ output_args ] = rename1( x,i,l )
j=1;
name=num2str(i);
list=dir;
for k=l:(l+x-1)
    movefile(list(l).name,strcat(name,'_',num2str(j),'.jpg'));
    j=j+1;
    l=l+1;
end
i=i+1;
save i;
save l;
end

