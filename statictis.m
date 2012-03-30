list1=dir('test/bin');
mysize1=size(list1);
count=0;
ii=0;
[names xx, yy] = textread('Nose.txt','%s %d %d');
count=0;
tt=0;
ttt=0;
for j=3:mysize1(1)
    tic
    I1=imreadbw(strcat('test/bin/',list1(j).name)) ; 
    disp(list1(j).name);disp(count);
    I1=I1-min(I1(:)) ;
    I1=I1/max(I1(:)) ;
    [frames1,descr1,gss1,dogss1] = sift( I1, 'Verbosity', 1 ) ;


    descr1=uint8(512*descr1) ;

    list=dir('matdata10');
mysize=size(list);
maxx=1;
maxxx=0;
maxname='';
maxxname='';
I=zeros(480,1280,1);
[names xx, yy] = textread('Nose.txt','%s %d %d');
for i=3:mysize(1)
    load(strcat('matdata10/',list(i).name),'descr2','I2','frames2');
     matches=siftmatch( descr1, descr2 ) ;

     ynose1=yy(find(strcmp(names,list1(j).name)));
     xnose1=xx(find(strcmp(names,list1(j).name)));
     ynose2=yy(find(strcmp(names,list(i).name(1:length(list(i).name)-4))));
     xnose2=xx(find(strcmp(names,list(i).name(1:length(list(i).name)-4))));
     out=plotmatches(ynose1,xnose1,ynose2,xnose2,I1,I2,frames1(1:2,:),frames2(1:2,:),matches) ;
     if(out>maxxx & out<maxx) maxxx=out; maxxxname=list(i).name;end;
     if(out>=maxx) maxxx=out; maxxx=maxx;maxx=out;maxxxname=maxxname;maxxname=list(i).name;end;
     
     progressbar(j/mysize1(1),i/mysize(1),count/213,ttt/100);
%drawnow ;
end
     if(strcmp(maxxname(1:length(maxxname)-10),list1(j).name(1:length(list1(j).name)-6)))
         count=count+1;end;
    clc
    tt=toc
    ttt=((j-3)*ttt+tt)/(j-2);
end
