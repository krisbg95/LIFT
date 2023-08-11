function [f] = clean_duplicates(x) 
for i=1:size(x,1)
    y{i} = strcat(x{i,1},x{i,2});
end 
for i = 1:size(x,1)
    comb = strcat(x{i,2},x{i,1}); 
    altcomb = find([y{1,:}] == comb);  
    if size(altcomb,2)>0 
        x{altcomb,1:2} = [0,0];
    end 
end 
z = find(x{:,1}=="0");
x(z',:) = []; 
f = x;
end