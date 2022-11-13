install.packages("magick")
library(magick)


image <- image_read("Table3.jpg")
res = dim(image)[2:1]
tiff("Table_3.tiff", 
     units="in", width=7, height=7, res=310, compression = 'lzw')
plot(image)
dev.off()

image1 <- image_read("Table4.jpg")
tiff("Table_4.tiff", 
     units="in", width=7, height=7, res=310, compression = 'lzw')
plot(image1)
dev.off()