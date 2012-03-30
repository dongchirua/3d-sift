function varargout = chuongtrinh(varargin)
% CHUONGTRINH MATLAB code for chuongtrinh.fig
%      CHUONGTRINH, by itself, creates a new CHUONGTRINH or raises the existing
%      singleton*.
%
%      H = CHUONGTRINH returns the handle to a new CHUONGTRINH or the handle to
%      the existing singleton*.
%
%      CHUONGTRINH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHUONGTRINH.M with the given input arguments.
%
%      CHUONGTRINH('Property','Value',...) creates a new CHUONGTRINH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before chuongtrinh_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to chuongtrinh_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help chuongtrinh

% Last Modified by GUIDE v2.5 09-Mar-2012 13:22:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @chuongtrinh_OpeningFcn, ...
                   'gui_OutputFcn',  @chuongtrinh_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before chuongtrinh is made visible.
function chuongtrinh_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to chuongtrinh (see VARARGIN)

% Choose default command line output for chuongtrinh
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes chuongtrinh wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = chuongtrinh_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)axes(handles.axes5)

[file_name file_path] = uigetfile ('*.abs;*.bmp;*.BMP;*.tif;*.TIF;*.jpg','Chon anh kiem tra ','data/bindata');
        if file_path ~= 0
            if(strcmp(file_name(length(file_name)-2:length(file_name)),'abs')==1)
                file_name=strcat(file_name(1:length(file_name)-3),'jpg');
            end
            I1 = imreadbw (strcat('data/bindata/',file_name));
            I2 = imread(strcat('data/colordata/',file_name));
        end
        
axes(handles.axes1)
imshow(I1);
axes(handles.axes2)
imshow(I2);
save I1;
save file_name;



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
%I1=imreadbw('bindata/90054d115.jpg') ; 
load I1;load file_name;
progressbar(0);

I1=I1-min(I1(:)) ;
I1=I1/max(I1(:)) ;
%I2=I2-min(I2(:)) ;
%I2=I2/max(I2(:)) ;

fprintf('Computing frames and descriptors.\n') ;
[frames1,descr1,gss1,dogss1] = sift( I1, 'Verbosity', 1 ) ;
%[frames2,descr2,gss2,dogss2] = sift( I2, 'Verbosity', 1 ) ;

% figure(11) ; clf ; 
%plotss(dogss1) ; colormap gray ;
% figure(12) ; clf ; 
%plotss(dogss2) ; colormap gray ;
%drawnow ;

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
descr1=uint8(512*descr1) ;
%descr2=uint8(512*descr2) ;
tic ; 
list=dir('matdata');
mysize=size(list);
maxx=1;
maxxx=0;
maxname='';
maxxname='';
I=zeros(480,1280,1);
[names xx, yy] = textread('Nose.txt','%s %d %d');
for i=3:mysize(1)
 %   if(strcmp(list(i).name(1:length(list(i).name)-10),file_name(1:length(file_name)-6)))
  %       continue;end;
    load(strcat('matdata5/',list(i).name),'descr2','I2','frames2');
     matches=siftmatch( descr1, descr2 ) ;
     fprintf('Matched in %.3f s\n', toc) ;

     ynose1=yy(find(strcmp(names,file_name)));
     xnose1=xx(find(strcmp(names,file_name)));
     ynose2=yy(find(strcmp(names,list(i).name(1:length(list(i).name)-4))));
     xnose2=xx(find(strcmp(names,list(i).name(1:length(list(i).name)-4))));
     out=plotmatches(ynose1,xnose1,ynose2,xnose2,I1,I2,frames1(1:2,:),frames2(1:2,:),matches) ;
     if(out>maxxx & out<maxx) maxxx=out; maxxxname=list(i).name;end;
     if(out>=maxx) maxxx=out; maxxx=maxx;maxx=out;maxxxname=maxxname;maxxname=list(i).name;end;
     progressbar(i/mysize(1));
%drawnow ;
end

%figure,imshow(imread(strcat('bindata/',maxxname)));
if(maxxx<20)
msgbox('Khong Tim Thay Doi Tuong');
end
axes(handles.axes3)
imshow(imread(strcat('bindata/',maxxxname(1:length(maxxxname)-4))));
axes(handles.axes4)
imshow(imread(strcat('colordata/',maxxxname(1:length(maxxxname)-4))));
