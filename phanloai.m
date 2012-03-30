list=dir('data/bindata/');
u=p;
for i=u:811107
    p=i;
    k=0;
    c=num2str(i);
    aa=[1 2 3 4]
    l=1
    t1=ceil(rand(1,1)*4);
    aa=[1:t1-1 t1+1:4]
    t=[t1]
    for j=1:l
        c
        t(j)
        strcat('test/',c,'_',num2str(t(j)),'.jpg.mat')
        strcat('data/bindata/',c,'_',num2str(t(j)),'.jpg')
        movefile(strcat('data/bindata/',c,'_',num2str(t(j)),'.jpg'),strcat('test/',c,'_',num2str(t(j)),'.jpg'))
        movefile(strcat('matdata10/',c,'_',num2str(t(j)),'.jpg.mat'),strcat('test/',c,'_',num2str(t(j)),'.jpg.mat'))
    end
end