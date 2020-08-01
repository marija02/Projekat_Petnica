
close all
clear all

N = 70;
folder = 'baza4\';

for i = 1:N
    % Ucitaj modifikovanu sliku 
    nameo = ['o' num2str(i,'%d')]; 
    imgo = (imread(strcat(folder,nameo), 'jpg'));
    
    names = ['s' num2str(i,'%d')];
    imgs = (imread(strcat(folder,names), 'jpg'));
    
    nameos = ['os' num2str(i,'%d')]; 
    if(ndims(imgo)>2)
        imgo = rgb2gray(imgo);
    end
    if(ndims(imgs)>2)
        imgs = rgb2gray(imgs);
    end
    img_diff = imbinarize(abs(imgo-imgs));
    
    imwrite(img_diff, strcat(folder, nameos, '.jpg'));
  
end