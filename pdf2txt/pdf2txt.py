#!/usr/bin/env python3

import pdf2image
import pytesseract as pat
import re


file = 'test.pdf'

p = re.compile('\n\n *')

# read file
imgs = pdf2image.convert_from_path(file)

file = ''
for n,im in enumerate(imgs):
    txt = pat.image_to_string(im,lang='chi_sim')

    txt = txt.replace('\n\x0c','\n\n')
    txt = re.sub(p,'###',txt)
    txt = re.sub('###*','###',txt)

    txt = txt.replace('\n','').replace('###','\n\n')
    file += '@Page-' + str(n+1) + '\n' + txt + '\n\n\n'


with open('xiang-si-wei.txt','wt') as f:
    f.write(file)


