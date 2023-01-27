function data_array = load_data(zeros_array, filename, precision)

if isfile(filename)
    dimensions = size(zeros_array);
    total_elements = numel(zeros_array);
    
    file = fopen(filename);
    data = fread(file, total_elements, precision);
    fclose(file);
    
    data_array = reshape(data, dimensions);
else
    data_array = zeros_array;
end

end