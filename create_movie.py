import cv2
import numpy as np
import glob
import os
import sys

img_array = []
found=True
i=1
dir=os.path.abspath('./refdata/MDX1342/Images/')
while found:
    filename=os.path.join(dir,'Figure_%05d.png'%i)
    if os.path.exists(filename):
        img = cv2.imread(filename)
        height, width, layers = img.shape
        size = (width, height)
        img_array.append(img)
        i+=1
    else:
        print('%s not found'%filename)
        break

out = cv2.VideoWriter(os.path.join(dir,'movie.avi'), cv2.VideoWriter_fourcc(*'DIVX'), 3, size)

for i in range(len(img_array)):
    out.write(img_array[i])
out.release()