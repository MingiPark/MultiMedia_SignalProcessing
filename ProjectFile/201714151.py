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
            '''
            f1 = open("dct_result.txt", 'a')
            sys.stdout = f1
            print(dct_blockLists[l])
sys.stdout = sys.__stdout__
f1.close()
'''
# 각 block에 대하여 Quantiztion 실행
quantization_table= [[16,11,10,16,24,40,51,61],
                     [12,12,14,19,26,58,60,55],
                     [14,13,16,24,40,57,69,56],
                     [14,17,22,29,51,87,80,62],
                     [18,22,37,56,68,109,103,77],
                     [24,35,55,64,81,104,113,92],
                     [49,64,78,87,103,121,120,101],
                     [72,92,95,98,112,100,103,99]]

quantization_blockList = np.zeros((1024,8,8))

for i in range(1024):
    for j in range(8):
        for k in range(8):
            quantization_blockList[i][j][k] = np.around(dct_blockLists[i][j][k] / quantization_table[j][k])
    '''
    f2 = open("quantized_block_result.txt", 'a')
    sys.stdout = f2
    print(quantization_blockList[i])

sys.stdout = sys.__stdout__
f2.close()
'''
# 각 block에 대하여 inverse quantized DCT coefficients 실행
inverse_qunatized_dct_blockList = np.zeros((1024,8,8))
for i in range(1024):
    for j in range(8):
        for k in range(8):
            inverse_qunatized_dct_blockList[i][j][k] = quantization_blockList[i][j][k] * quantization_table[j][k]
    '''
    f3 = open("inverse_quantized_dct_block_result.txt", 'a')
    sys.stdout =f3
    print(inverse_qunatized_dct_blockList[i])

sys.stdout = sys.__stdout__
f3.close()
'''

# 3-(a)에 대한 idct
a_idct = np.zeros((1024,8,8))
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
                    sum += inverse_qunatized_dct_blockList[l][u][v]*np.cos(((2*i+1)*u*math.pi)/16)*np.cos(((2*j+1)*v*math.pi)/16)*(1/4)*Cu*Cv
            a_idct[l][i][j] = sum
            a_idct= np.clip(a_idct, 0, 255)
            a_idct = a_idct.astype(c_ubyte)
'''
# 3-(a)에 대한 reconstruct image
a_reconstruct_image = Image.new('L',(256,256))
n = 0
for i in range(0,256,8):
    for j in range(0,256,8):
        b_image = Image.fromarray(a_idct[n])
        area = (j,i,j+8,i+8)
        a_reconstruct_image.paste(b_image,area)
        n += 1
a_reconstruct_image.save("3(a)_reconstruct_image.jpg")
'''
# 3-(b) F(0,0)유지 나머지 모두 truncate
b = np.zeros((1024,8,8))
for i in range(1024):
    for j in range(8):
        for k in range(8):
            if j is 0 and k is 0:
                b[i][j][k] = inverse_qunatized_dct_blockList[i][j][k]
            else:
                b[i][j][k] = 0
    '''
    f4 = open("truncate_a.txt", 'a')
    sys.stdout = f4
    print(a[i])

sys.stdout = sys.__stdout__
f4.close()
'''
# 3-(c) F(0,0), F(0,1), F(1,0) 유지 나머지 모두 truncate
c = np.zeros((1024,8,8))
for i in range(1024):
    for j in range(8):
        for k in range(8):
            if (j is 0 and k is 0) or (j is 0 and k is 1) or (j is 1 and k is 0):
                c[i][j][k] = inverse_qunatized_dct_blockList[i][j][k]
            else:
                c[i][j][k] = 0
    '''
    f5 = open("truncate_b.txt", 'a')
    sys.stdout = f5
    print(b[i])

sys.stdout = sys.__stdout__
f5.close()
'''

# 3-(d) F(0,0), F(1,0), F(0,1), F(1,1), F(2,0), F(0,2) 유지 나머지 모두 truncate
d = np.zeros((1024,8,8))
for i in range(1024):
    for j in range(8):
        for k in range(8):
            if (j is 0 and k is 0) or (j is 1 and k is 0) or \
                    (j is 0 and k is 1) or (j is 1 and k is 1) or (j is 2 and k is 0) or (j is 0 and k is 2):
                d[i][j][k] = inverse_qunatized_dct_blockList[i][j][k]
            else:
                d[i][j][k] = 0
    '''
    f6 = open("truncate_c.txt", 'a')
    sys.stdout = f6
    print(c[i])

sys.stdout = sys.__stdout__
f6.close()
'''

# 3-(d) F(0,0), F(1,0), F(0,1), F(1,1), F(2,0), F(0,2), F(0,3), F(1,2), F(2,1), F(3,0) 유지 나머지 모두 truncate
e = np.zeros((1024,8,8))
for i in range(1024):
    for j in range(8):
        for k in range(8):
            if (j is 0 and k is 0) or (j is 1 and k is 0) or\
                    (j is 0 and k is 1) or (j is 1 and k is 1) or (j is 2 and k is 0) or\
                    (j is 0 and k is 0) or (j is 1 and k is 2) or (j is 2 and k is 1) or\
                    (j is 3 and k is 0) or (j is 0 and k is 3):
                e[i][j][k] = inverse_qunatized_dct_blockList[i][j][k]
            else:
                e[i][j][k] = 0
    '''
    f7 = open("truncate_d.txt", 'a')
    sys.stdout = f7
    print(d[i])

sys.stdout = sys.__stdout__
f7.close()
'''

# 3-(b)에 대한 idct
b_idct = np.zeros((1024,8,8))
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
                    sum += b[l][u][v]*np.cos(((2*i+1)*u*math.pi)/16)*np.cos(((2*j+1)*v*math.pi)/16)*(1/4)*Cu*Cv
            b_idct[l][i][j] = sum
            b_idct= np.clip(b_idct, 0, 255)
            b_idct = b_idct.astype(c_ubyte)
'''
# 3-(b)에 대한 reconstruct image
b_reconstruct_image = Image.new('L',(256,256))
n = 0
for i in range(0,256,8):
    for j in range(0,256,8):
        b_image = Image.fromarray(b_idct[n])
        area = (j,i,j+8,i+8)
        b_reconstruct_image.paste(b_image,area)
        n += 1
b_reconstruct_image.save("3(b)_reconstruct_image.jpg")
'''
# 3-(c)에 대한 idct
c_idct = np.zeros((1024,8,8))
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
                    sum += c[l][u][v]*np.cos(((2*i+1)*u*math.pi)/16)*np.cos(((2*j+1)*v*math.pi)/16)*(1/4)*Cu*Cv
            c_idct[l][i][j] = sum
            c_idct = np.clip(c_idct, 0, 255)
            c_idct = c_idct.astype(c_ubyte)
'''
# 3-(c)에 대한 reconstruct image
c_reconstruct_image = Image.new('L',(256,256))
n = 0
for i in range(0,256,8):
    for j in range(0,256,8):
        c_image = Image.fromarray(c_idct[n])
        area = (j,i,j+8,i+8)
        c_reconstruct_image.paste(c_image,area)
        n += 1
c_reconstruct_image.save("3(c)_reconstruct_image.jpg")
'''
# 3-(d)에 대한 idct
d_idct = np.zeros((1024,8,8))
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
                    sum += d[l][u][v]*np.cos(((2*i+1)*u*math.pi)/16)*np.cos(((2*j+1)*v*math.pi)/16)*(1/4)*Cu*Cv
            d_idct[l][i][j] = sum
            d_idct = np.clip(d_idct, 0, 255)
            d_idct = d_idct.astype(c_ubyte)
'''
# 3-(d)에 대한 reconstruct image
d_reconstruct_image = Image.new('L',(256,256))
n = 0
for i in range(0,256,8):
    for j in range(0,256,8):
        d_image = Image.fromarray(d_idct[n])
        area = (j,i,j+8,i+8)
        d_reconstruct_image.paste(d_image,area)
        n += 1
d_reconstruct_image.save("3(d)_reconstruct_image.jpg")
'''
# 3-(e)에 대한 idct
e_idct = np.zeros((1024,8,8))
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
                    sum += e[l][u][v]*np.cos(((2*i+1)*u*math.pi)/16)*np.cos(((2*j+1)*v*math.pi)/16)*(1/4)*Cu*Cv
            e_idct[l][i][j] = sum
            e_idct = np.clip(e_idct, 0, 255)
            e_idct = e_idct.astype(c_ubyte)
'''
# 3-(e)에 대한 reconstruct image
e_reconstruct_image = Image.new('L',(256,256))
n = 0
for i in range(0,256,8):
    for j in range(0,256,8):
        e_image = Image.fromarray(e_idct[n])
        area = (j,i,j+8,i+8)
        e_reconstruct_image.paste(e_image,area)
        n += 1
e_reconstruct_image.save("3(e)_reconstruct_image.jpg")
'''
#3-(a) MSE 계산
MSE_A = 0
MSE_A = np.square(np.subtract(dct_blockList,a_idct)).mean(axis=None)
print("3-(a) MSE reult : " + str(MSE_A))

#3-(b) MSE 계산
MSE_B = 0
MSE_B = np.square(np.subtract(dct_blockList,b_idct)).mean(axis=None)
print("3-(b) MSE reult : " + str(MSE_B))

#3-(c) MSE 계산
MSE_C = 0
MSE_C = np.square(np.subtract(dct_blockList,c_idct)).mean(axis=None)
print("3-(c) MSE reult : " + str(MSE_C))

#3-(d) MSE 계산
MSE_D = 0
MSE_D = np.square(np.subtract(dct_blockList,d_idct)).mean(axis=None)
print("3-(d) MSE reult : " + str(MSE_D))

#3-(e) MSE 계산
MSE_E = 0
MSE_E = np.square(np.subtract(dct_blockList,e_idct)).mean(axis=None)
print("3-(e) MSE reult : " + str(MSE_E))

#3-(a) PSNR 계산
PSNR_A = 0
PSNR_A = 10 * math.log10(255**2 / MSE_A)
print("3-(A) PSNR result : " + str(PSNR_A) + "(dB)")

#3-(b) PSNR 계산
PSNR_B = 0
PSNR_B = 10 * math.log10(255**2 / MSE_B)
print("3-(B) PSNR result : " + str(PSNR_B) + "(dB)")

#3-(c) PSNR 계산
PSNR_C = 0
PSNR_C = 10 * math.log10(255**2 / MSE_C)
print("3-(C) PSNR result : " + str(PSNR_C) + "(dB)")

#3-(d) PSNR 계산
PSNR_D = 0
PSNR_D = 10 * math.log10(255**2 / MSE_D)
print("3-(D) PSNR result : " + str(PSNR_D) + "(dB)")


#3-(e) PSNR 계산
PSNR_E = 0
PSNR_E = 10 * math.log10(255**2 / MSE_E)
print("3-(E) PSNR result : " + str(PSNR_E) + "(dB)")