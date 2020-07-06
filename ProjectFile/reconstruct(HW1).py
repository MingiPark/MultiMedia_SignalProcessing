# -*- coding:utf-8 -*-
from PIL import Image
import numpy as np
from ctypes import c_ubyte
import math


img = Image.open('original_image.jpeg')

dct_blockList = []

# 8x8 block 생성 -> 1024개
for i in range(32):
    for j in range(32):
        box = (j * 8, i * 8, (j + 1) * 8, (i + 1) * 8)
        block = img.crop(box)
        blocks = np.array(block)
        dct_blockList.append(blocks)

dct_blockLists = np.zeros((1024,8,8))

# 각 block에 대하여 dct 계산
for l in range(1024):
    for u in range(8):
        for v in range(8):
            sum=0.0
            for i in range(8):
                for j in range(8):
                    sum += dct_blockList[l][i][j]*np.cos(((2*i+1)*u*math.pi)/16)*np.cos(((2*j+1)*v*math.pi)/16)
                    if u == 0:
                        Cu = 1 / math.sqrt(2)
                    else:
                        Cu = 1
                    if v == 0:
                        Cv = 1 / math.sqrt(2)
                    else:
                        Cv = 1
            dct_blockLists[l][u][v] = (1/4)*(Cu*Cv*sum)
            dct_blockLists = dct_blockLists.astype(int)

idct = np.zeros((1024,8,8))
for l in range(1024):
    for i in range(8):
        for j in range(8):
            sum = 0
            for u in range(8):
                for v in range(8):
                    if u == 0:
                        Cu = 1 / math.sqrt(2)
                    else:
                        Cu = 1
                    if v == 0:
                        Cv = 1 / math.sqrt(2)
                    else:
                        Cv = 1
                    sum += dct_blockLists[l][u][v]*np.cos(((2*i+1)*u*math.pi)/16)*np.cos(((2*j+1)*v*math.pi)/16)*(1/4)*Cu*Cv
            idct[l][i][j] = sum
            idct= np.clip(idct, 0, 255)
            idct = idct.astype(c_ubyte)

idct_image = Image.new('L',(256,256))
n = 0
for i in range(0,256,8):
    for j in range(0,256,8):
        image = Image.fromarray(idct[n])
        area = (j,i,j+8,i+8)
        idct_image.paste(image,area)
        n += 1
idct_image.save("idct_image.jpg")
