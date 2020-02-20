from PIL import Image

nFiles = 26

for i in range(0,nFiles,1):
    fileName=str(i)+'_weno.png'
    im=Image.open(fileName)
    #print im.size
    #print im.getbbox()
    #im2 = im.crop(im.getbbox())
    sizeX = 6600
    sizeY = 2085
    cropX=(im.size[0]-sizeX)*0.5
    cropY=(im.size[1]-sizeY)*0.5

    offsetX = 150

    #(x0,y0,xN,yN)
    cropSize=(cropX-offsetX, cropY, cropX-offsetX+sizeX, cropY+sizeY)
        
    im2=im.crop(cropSize)
    fileName=str(i)+'_weno_crop.png'
    im2.save(fileName)
