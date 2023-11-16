r"""°°°
<span style="color:yellow">two kernel based extract random tile from wsi</span>
°°°"""
# |%%--%%| <jMbkpUz6xD|CVkE2fNOmq>

# run based "cellprofiler" kernel
import os
from histolab.slide import Slide
from histolab.tiler import RandomTiler

# |%%--%%| <CVkE2fNOmq|vko5qaGiOO>

os.chdir('/data/dk/crc/raw')
path = '/data/dk/crc/raw'
processed_path = '/home/ma/dk/crc/randomtile/'
file_list = os.listdir(path)
for file_name in file_list:
    if file_name != 'FS230104002-317_2023-07-02_11_38_48.svs':
        print(file_name)
        file_path = os.path.join(path, file_name)
        sample  = Slide(file_path, processed_path)
        print(f"Slide name: {sample.name}")
        print(f"Levels: {sample.levels}")
        print(f"Dimensions at level 0: {sample.dimensions}")
        print(f"Dimensions at level 1: {sample.level_dimensions(level=1)}")
        random_tiles_extractor = RandomTiler(
        tile_size=(512, 512),
        n_tiles=20,
        level=1,
        seed=42,
        check_tissue=True, # default
        tissue_percent=40.0, # default
        prefix=os.path.splitext(file_name)[0]+'_', # save tiles in the "random" subdirectory of slide's processed_path
        suffix=".png")# default
        random_tiles_extractor.extract(sample)
        random_tiles_extractor.locate_tiles(
        slide=sample,
        linewidth=2,
        scale_factor=24,  # default
        alpha=128,  # default
        outline="#40E0D0" # default
        )


# |%%--%%| <vko5qaGiOO|jXyMXk9awR>

a = random_tiles_extractor.locate_tiles(
slide=sample,
linewidth=4,
scale_factor=24,  # default
alpha=128,  # default
outline="#40E0D0" # default
)
type(a)
a.save('a.png')
# |%%--%%| <jXyMXk9awR|KMBbwzGHa6>
r"""°°°
<span style="color:yellow">'vahadane' kernel based normalize tiles</span>
°°°"""
# |%%--%%| <KMBbwzGHa6|MeOQRw3wNr>

import cv2
import matplotlib.pyplot as plt
import torchstain
import torch
from torchvision import transforms
import time
from PIL import Image
import numpy as np



# |%%--%%| <MeOQRw3wNr|ImQEW0huCK>

import os 
os.chdir("/home/ma/torchstain/")
root_dir = '/data/dk/crc/'

# |%%--%%| <ImQEW0huCK|cUKc4oIF46>

size = 1024
# target = cv2.resize(cv2.cvtColor(cv2.imread("/data/dk/ucec/randomtile/TCGA-A5-A0GA-01Z-00-DX1.D691BFC9-05DF-4815-9DCE-84B343C7015A_tile_8_level1_29636-35340-31684-37388.png"), cv2.COLOR_BGR2RGB), (size, size))
target = cv2.resize(cv2.cvtColor(cv2.imread("./data/target.png"), cv2.COLOR_BGR2RGB), (size, size))
normalizer = torchstain.normalizers.MacenkoNormalizer(backend='numpy')
normalizer.fit(target)
for file in os.listdir('/data/dk/crc/randomtile/'):
    print(file)
    try:
        to_transform = cv2.resize(cv2.cvtColor(cv2.imread(os.path.join(root_dir,'randomtile/',file)), cv2.COLOR_BGR2RGB), (size, size))
        T = transforms.Compose([
        transforms.ToTensor(),
        transforms.Lambda(lambda x: x*255)
        ])
        torch_normalizer = torchstain.normalizers.MacenkoNormalizer(backend='torch')
        torch_normalizer.fit(T(target))
        norm, H, E = normalizer.normalize(I=to_transform, stains=True)
        filename = os.path.join('/data/dk/crc/normalizedtile/',file)
        norm2  = cv2.cvtColor(norm, cv2.COLOR_RGB2BGR)
        cv2.imwrite(filename, norm2)
    except Exception as e:
        print(f"delete {file} :{str(e)}")
        remove_path = os.path.join('/data/dk/crc/randomtile/',file)
        os.remove(remove_path)


# |%%--%%| <cUKc4oIF46|XCYPN45SsU>

img = Image.open('/home/ma/dk/crc/CLAM/output_image.png').convert('RGB')
target = cv2.resize(cv2.cvtColor(cv2.imread("/home/ma/torchstain/data/target.png"), cv2.COLOR_BGR2RGB), (1024, 1024))
to_transform = cv2.resize(np.array(img), (1024, 1024))
normalizer = torchstain.normalizers.MacenkoNormalizer(backend='numpy')
normalizer.fit(target)

T = transforms.Compose([
transforms.ToTensor(),
transforms.Lambda(lambda x: x*255)
])
torch_normalizer = torchstain.normalizers.ReinhardNormalizer(backend='torch')
torch_normalizer.fit(T(target))
norm, H, E = normalizer.normalize(I=to_transform, stains=True)
plt.imshow(norm)
plt.savefig('/data/dk/crc/test/target.png')

# |%%--%%| <XCYPN45SsU|8grrZEkRzm>

img = "/home/ma/torchstain/data/target.png"
def is_empty(img):
   # Reading Image
   image = cv2.imread(img, 0)
   np.reshape(image, (-1,1))
   u, count_unique = np.unique(image, return_counts =True)
   
   if count_unique.size< 10:
      return "Image is empty"
   else:
      return "Image is not empty"

print(is_empty('whh.jpg'))

# |%%--%%| <8grrZEkRzm|mpn2lam652>

image = cv2.imread('/home/ma/dk/crc/CLAM/output_image.png')
gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
threshold = 200 
white_pixel_count = np.sum(gray_image > threshold)
total_pixels = gray_image.size
white_percentage = (white_pixel_count / total_pixels) * 100

white_threshold_percentage = 90  # Adjust this value as needed
if white_percentage >= white_threshold_percentage:
    print("The image is mostly white.")
else:
    print("The image is not mostly white.")

# |%%--%%| <mpn2lam652|yCVwqOJa9U>
r"""°°°
<span style="font-size: 26px;color:red">Tile Feature Extract</span><br>
<span style="color:yellow">after extract the random tile and normalize the image , run cellprofiler as a python package to extract morphology of wsi</span><br>
kernel : 'cf'
°°°"""
# |%%--%%| <yCVwqOJa9U|dlAYuNLXGq>

import cellprofiler_core.pipeline
import cellprofiler_core.preferences
import cellprofiler_core.utilities.java
import pathlib
import os
cellprofiler_core.preferences.set_headless()
cellprofiler_core.utilities.java.start_java()

# |%%--%%| <dlAYuNLXGq|kuVNvFFE4F>

os.chdir('/data/dk/crc/')

# |%%--%%| <kuVNvFFE4F|1EQNEtsVYX>


pipeline = cellprofiler_core.pipeline.Pipeline()
pipeline.load("/data/dk/cf/pathomic.cppipe")
current_dir = pathlib.Path().absolute()
cellprofiler_core.preferences.set_default_output_directory(f"{current_dir}/cf_out")
file_list = list(pathlib.Path('.').absolute().glob('normalizedtile/*.png'))
files = [file.as_uri() for file in file_list]
pipeline.read_file_list(files)

# |%%--%%| <1EQNEtsVYX|YUgmjbwY21>

output_measurements = pipeline.run()


# |%%--%%| <YUgmjbwY21|KTjQbpDieJ>

cellprofiler_core.utilities.java.stop_java()

# |%%--%%| <KTjQbpDieJ|neNhGlLnD8>


