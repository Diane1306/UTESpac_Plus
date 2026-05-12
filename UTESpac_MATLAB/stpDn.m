function output = stpDn(data,rows)

ouput = nan(floor(size(data,1)/2),size(data,2));
for i = rows:rows:length(data)
    output(floor(i/rows),:) = mean(data(i-rows+1:i,:), 'omitmissing');
end

end

