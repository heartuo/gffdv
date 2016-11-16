# Gaussian File Format Data Visualizer

This project applies [yt project](http://yt-project.org/) to visualize the process of a hydrogen atom passing through aluminum foil, specifically, the changing of spatial electronic densities. Read the comments in the script, if you wish to re-use part of the code to apply to more visualizations in the future.

### Usage
The file "volume_render_ss.py" reads in a frame number and a source data file (.cub), which is in [Gaussian file format](http://paulbourke.net/dataformats/cube/), and generates one spherical image. The command is in the format of the following line:
```sh
$ python volume_render_ss.py <data_file_path> <frame_num>
```
Example usage:
```sh
$ python volume_render_ss.py density00000.cub 0
```
The process of generating a series of frames could be automated by a shell script. 

### Compiling into A 3-D Animation and Uploading to YouTube

1. Compile the sterero-spherical images into a video file (e.g. .mp4, .mov), using any tool of your choice (e.g. ImageJ, Adobe Premiere, iMovie).

2. Add the metadata required by YouTube to the file. Use [this](https://github.com/google/spatial-media/blob/master/spatialmedia/README.md) Python script to process the file. Read the official YouTube [support](https://support.google.com/youtube/answer/6178631?hl=en) on this if you wish. 

3. Upload the video file. 

###### This software was developed with support from NCSA through the Students Pushing Innovation (SPIN) program and this material is based upon work supported by the National Science Foundation under Grant No. DMR-1555153.
