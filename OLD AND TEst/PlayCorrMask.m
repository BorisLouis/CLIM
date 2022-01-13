% play corrmask


while true
   figure(1) 
   for i = 1:size(allCorrMask,3)
       
       imagesc(allCorrMask(:,:,i))
       title(['Threshold = ' num2str(threshold(i))])
       axis image
       colormap('jet')
       pause(0.200)
       
   end
end

       