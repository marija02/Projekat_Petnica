function ret = zaokruziRazliku(bw, L1,L2, p)
% Binarizuje sliku po svakoj komponenti po blokovima
% vraca 0 i 1 u Y

ret = imdilate(bw, strel('rectangle', [20 20]));
%imadilate spaja ivice, rupice
imshow(ret);

[b_h b_w] = size(bw);
b_h = b_h / L1;
b_w = b_w / L2;

ret = zeros(1,L1*L2);

for i = 1:L1
    for j = 1:L2
        suma = sum(sum(bw((i-1)*b_h + 1:i*b_h, (j-1)*b_w + 1:j*b_w)));
        if (suma / (b_w * b_h) > p)
            ret(i*j) = 1;
        end
    end
end

