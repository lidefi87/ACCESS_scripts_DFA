import os
import imageio

#Set file path where images are located
filein = r'C:\Users\ldfierro\Documents\ACCESS\ACCESS_scripts_DFA\Outputs\Figures\Maps\SeasonalSeaIce\IncObservations'
#Set the file path for the video to be save, include the name of the video and its extension
fileout = os.path.join(filein, 'VideoSeasonalSeaIce2000-17.gif')

#Get list of filepaths for all figures to be used in movie
filenames = os.listdir(filein)

#Create an empty list to save all read images
images = []

#Loop to read all images
for file in filenames:
    images.append(imageio.imread(os.path.join(filein, file)))

#Save video using path specified at the beginning. Each image to be shown for 0.5 seconds
imageio.mimsave(fileout, images, duration = 0.5)