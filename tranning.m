clear all
list=dir('data/bindata');
mysize=size(list);
for i=3:mysize(1)
       disp(list(i).name);
     
     I2=imreadbw(strcat('data/bindata/',list(i).name)) ; 
      


%I1=imreadbw('bindata/90405d34.jpg') ; 


%I1=I1-min(I1(:)) ;
%I1=I1/max(I1(:)) ;
I2=I2-min(I2(:)) ;
I2=I2/max(I2(:)) ;

fprintf('Computing frames and descriptors.\n') ;
%[frames1,descr1,gss1,dogss1] = sift( I1, 'Verbosity', 1 ) ;
[frames2,descr2,gss2,dogss2] = sift( I2, 'Verbosity', 1 ) ;

% figure(11) ; clf ; 
%plotss(dogss1) ; colormap gray ;
% figure(12) ; clf ; 
plotss(dogss2) ; colormap gray ;
drawnow ;

% figure(2) ; clf ;
%subplot(1,2,1) ; imagesc(I1) ; colormap gray ;
%hold on ;
%h=plotsiftframe( frames1 ) ; set(h,'LineWidth',2,'Color','g') ;
%h=plotsiftframe( frames1 ) ; set(h,'LineWidth',1,'Color','k') ;

%subplot(1,2,2) ; imagesc(I2) ; colormap gray ;
%hold on ;
%h=plotsiftframe( frames2 ) ; set(h,'LineWidth',2,'Color','g') ;
%h=plotsiftframe( frames2 ) ; set(h,'LineWidth',1,'Color','k') ;

fprintf('Computing matches.\n') ;
% By passing to integers we greatly enhance the matching speed (we use
% the scale factor 512 as Lowe's, but it could be greater without
% overflow)
%descr1=uint8(512*descr1) ;
descr2=uint8(512*descr2) ;
tic ; 
%matches=siftmatch( descr1, descr2 ) ;
%fprintf('Matched in %.3f s\n', toc) ;

clf ;
[ynose xnose]=getnose(I2);
save(strcat('matdata8/',list(i).name,'.mat'),'descr2','I2','frames2','ynose','xnose');
%load('matdata/abc.mat','descr2','I2','frames2');
%descr2=descr;
%I2=I;
%frames2=frames;
%plotmatches(I1,I2,frames1(1:2,:),frames2(1:2,:),matches) ;
%drawnow ;
    clc;
    disp(i);
    fprintf('xong %d %% \n', round(i/mysize(1)*100));
end
disp('xong');