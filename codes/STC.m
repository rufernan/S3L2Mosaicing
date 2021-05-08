%% Función de arbol STC-S3
function [result, index] = STC(first_point, second_point, first_product, second_product)

    bin_min = zeros(1,20);            
    first_bin = flip(dec2bin(first_point)-'0'); % Data to binary compounds
    first_bin(length(bin_min)) = 0;
    second_bin = flip(dec2bin(second_point)-'0'); % Data to binary compounds
    second_bin(length(bin_min)) = 0;
    
    if (first_bin(3) == 1 || second_bin(3) == 1) % Any pixel is LAND?
        if (first_bin(3) == 1 && second_bin(3) == 1) % Both pixels are LAND?
            if (first_bin(18) ~= 1 || first_bin(17) ~= 1 || second_bin(18) ~= 1 || second_bin(17) ~= 1) % Pixels not classified as water,cloud or ice in OGVI
                if (first_bin(18) ~= 1 && first_bin(17) ~= 1 && second_bin(18) ~= 1 && second_bin(17) ~= 1)
                    if (first_product >= second_product) %OGVI1 >= OGVI2
                        result = first_product;
                        index = 2;
                    else
                        result = second_product;
                        index = 1;
                    end
                end
            elseif (first_bin(18) ~= 1 || first_bin(17) ~= 1) % First is not water, cloud or ice in OGVI
                result = first_product;
                index = 2;
            elseif (second_bin(18) ~= 1 || second_bin(17) ~= 1) % Second is not water, cloud or ice in OGVI
                result = second_product;
                index = 1;
            elseif (first_product >= second_product) % OGVI1 >= OGVI2
                result = first_product;
                index = 2;
            else
                result = second_product;
                index = 1;
            end
        end
        if (first_bin(3) == 1 || first_bin(5) == 1) % First is water or ice
            result = second_product;
            index = 1;
        else 
            result = first_product;
            index = 2;
        end
    end
    if (first_bin(5) == 1 || second_bin(5) == 1)
        if (first_bin(5) == 1 && second_bin(5) == 1)
            if (first_product >= second_product)
                result = first_product;
                index = 2;
            else
                result = second_product;
                index = 1;
            end
        end
        if (first_bin(2) == 1)
            result = second_product;
            index = 1;
        else
            result = first_product;
            index = 2;
        end
    end
    if (first_bin(2) == 1 || second_bin(2) == 1)
        if (first_bin(18) == 1 && second_bin(18) ~= 1)
            result = first_product;
            index = 2;
        elseif (first_bin(18) ~= 1 && second_bin(18) == 1)
            result = second_product;
            index = 1;
        elseif (first_product >= second_product)
            result = first_product;
            index = 2;
        else
            result = second_product;
            index = 1;
        end
        if (first_product >= second_product)
            result = first_product;
            index = 2;
        else 
            result = second_product;
            index = 1;
        end
    end
end
    