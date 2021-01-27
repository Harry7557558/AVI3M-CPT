# see if this hack works

from skimage.restoration import (denoise_tv_chambolle, denoise_bilateral,
                                 denoise_wavelet, estimate_sigma)
from skimage import io, img_as_float, img_as_ubyte

image = img_as_float(io.imread("test.png"))

print(image)

io.imsave("test1.png", img_as_ubyte(image))

image = denoise_tv_chambolle(image, weight=0.1, multichannel=True)

print(image)
io.imsave("test2.png", img_as_ubyte(image))
