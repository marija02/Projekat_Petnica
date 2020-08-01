function ret = binarizujRazliku(img_o, img_m, L, p)
% Binarizuje sliku po svakoj komponenti po blokovima

[b_h b_w] = size(img_o); % Block height and block width
b_h = b_h / L;
b_w = b_w /L;

% Ako je p procenata piksela u bloku 1 obelezi ceo blok kao 1
ret = zeros(L, L, 3);

for i=1:3
    bw = imbinarize(abs(img_m(:, :, i)-img_o(:, :, i)), 0);
    for j = 1:L
        for k = 1:L
            num_of_ones = sum(sum(bw((j - 1)*b_h + 1:j*b_h , (i - 1)*b_w + 1:i*b_w)));
            if (num_of_ones / (b_h*b_w) >= p)
                ret(j, k, i) = 1;
            end
        end
    end
end

