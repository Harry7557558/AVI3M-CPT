import os
from PIL import Image

files = [s for s in os.listdir() if s[-4:]=='.png']

for file in files:
    img = Image.open(file).convert('RGB')
    img.save(file.replace('.png','.jpg'))
