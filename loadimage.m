function image1=loadimage

[file_name file_path] = uigetfile ('*.abs;*.bmp;*.BMP;*.tif;*.TIF;*.jpg','Chon anh kiem tra ','data');
        if file_path ~= 0
            if(strcmp(file_name(length(file_name)-2:length(file_name)),'abs')==1)
                file_name=strcat(file_name(1:length(file_name)-3),'jpg');
                file_path='data/bindata/';
            end
            image1 = imreadbw ([file_path,file_name]);
        end
      figure,  imshow(image1);



