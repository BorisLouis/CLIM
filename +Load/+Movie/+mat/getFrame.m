function Im = getFrame(filePath,frames)

    data = load(filePath);
    field = fieldnames(data);
    
    data = data.(field{1});
    
    Im = data(:,:,frames);




end