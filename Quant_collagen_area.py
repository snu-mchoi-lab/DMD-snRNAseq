import os
from PIL import Image, ImageDraw
import numpy as np

def count_color_pixels(image_path):
    img = Image.open(image_path)
    img_np = np.array(img)

    red_pixels_count = 0
    blue_pixels_count = 0
    red_pixels_coords = []
    blue_pixels_coords = []

    if len(img_np.shape) == 3:  # Check if the image is RGB
        red_channel = img_np[:, :, 0]
        green_channel = img_np[:, :, 1]
        blue_channel = img_np[:, :, 2]

        for i in range(img_np.shape[0]):
            for j in range(img_np.shape[1]):
                if (red_channel[i, j] < 120 and blue_channel[i, j] > 50) and not (30 < red_channel[i, j] < 100 and blue_channel[i, j] < 50):
                    blue_pixels_count += 1
                    blue_pixels_coords.append((j, i))
                elif blue_channel[i, j] < 210:
                    red_pixels_count += 1
                    red_pixels_coords.append((j, i))

    ratio_blue_to_red = blue_pixels_count / red_pixels_count if red_pixels_count > 0 else 0
    return red_pixels_count, blue_pixels_count, ratio_blue_to_red, red_pixels_coords, blue_pixels_coords

def process_images_in_folder(folder_path, output_folder_path):
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    for filename in os.listdir(folder_path):
        if filename.lower().endswith(('.tiff', '.tif', '.jpg')):
            image_path = os.path.join(folder_path, filename)
            red_count, blue_count, ratio, red_coords, blue_coords = count_color_pixels(image_path)

            img = Image.open(image_path)
            draw = ImageDraw.Draw(img)
            for coord in red_coords:
                draw.point(coord, fill='red')
            for coord in blue_coords:
                draw.point(coord, fill='blue')

            output_image_path = os.path.join(output_folder_path, filename)
            img.save(output_image_path)

# Example usage
folder_path = 'your folder path'
output_folder_path = 'your output folder path'
process_images_in_folder(folder_path, output_folder_path)
