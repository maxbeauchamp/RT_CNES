
import glob
from PIL import Image


# filepaths
fp_in = "./png/ssh_reconstructions_2017-*.png"
fp_out = "./png/animation_ssh_reconstructions.gif"
# filepaths
#fp_in = "./png/velocity_reconstructions_2017-*.png"
#fp_out = "./png/animation_velocity_reconstructions.gif"
# filepaths
#fp_in = "./png/vorticity_reconstructions_*.png"
#fp_out = "./png/animation_vorticity_reconstructions.gif"
#fp_in = "./png/ssh_sst_reconstructions_*.png"
#fp_out = "./png/animation_ssh_sst_reconstructions.gif"
print(sorted(glob.glob(fp_in)))
# https://pillow.readthedocs.io/en/stable/handbook/image-file-formats.html#gif
img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
img.save(fp=fp_out, format='GIF', append_images=imgs,
         save_all=True, duration=200, loop=0)