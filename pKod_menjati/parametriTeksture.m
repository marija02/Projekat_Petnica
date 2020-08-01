function ret = parametriTeksture(img, L1, L2)
%TextureParam prosledjuje kolonu koja odgovara slici, L1 i L2 blokovi
ret = zeros(12, L1*L2);
[b_h b_w] = size(img);
b_h = b_h / L1; %b_h da se vidi koliko je visok blok
b_w = b_w / L2;

for i = 1:L1
    for j=1:L2
        mask = zeros(size(img));
        mask((i-1)*b_h+1:i*b_h, (j-1)*b_w+1:j*b_w) = 1; % pravljenje maske
        [SRE,LRE,GLN,RP,RLN,LGRE,HGRE] = glrlm(img,16,mask);
        
        BW = edge(img,'canny');
        BW(~mask) = 0;

        gustina_piksela = sum(BW(:))/L1*L2;
        blok_img = img((i-1)*b_h+1:i*b_h, (j-1)*b_w+1:j*b_w);
        blok_img = im2double(blok_img(:));
        v2 = var(blok_img(:)); %vraca variance svih piksela
        m = mean(blok_img(:));
        %mean = sum(x)/length(x)
        %variance = sum((x - mean(x)).^2)/(length(x) - 1);
        k = kurtosis(double(blok_img(:)));
        s = skewness(blok_img(:));
        
        ret_v = [SRE,LRE,GLN,RP,RLN,LGRE,HGRE, gustina_piksela, v2, m, k, s];
        ret(:, i*j) = (ret_v');
    end
end

ret = mean(ret,2);
