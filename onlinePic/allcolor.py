from PIL import Image
import numpy as np

colorlist = []
for t in range(16):
    beg = t * 16
    end = beg + 16
    for k in range(256):
        ls = [ [i,j,k] for i in range(beg,end) for j in range(256) ]
        colorlist.append(ls)

im = Image.fromarray(np.uint8(colorlist))
im.save('allcolor.bmp')

print('done')

