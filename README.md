# **Diffusion-Model**
**Animated 3D Particle Diffusion**

This model uses a random walk approach to model and plot the diffusion of particles within a rectangular channel. Using principles from Fickian diffusion, each particle's path is randomly generated and then animated over time. Longitudinal flow patterns can be applied as well, with velocity profiles calculated using the log-law equation. Particles can be released either all at once or staggered over time. Coding was done in Anaconda's Spyder environment using the Matplotlib and Numpy libraries. Saving animations requires them to be written using FFMpeg.

### Path Generation 
Particle paths are created by first drawing step distances from a Gaussian Distribution with a standard deviation of sqrt(2Dt) as described by Fischer et al. (1979). A polar angle (phi) and azimuthal angle (theta) are then selected from a uniform random distribution between 0 and 2pi. Using the step distance as the radial distance in a spherical coordinate system and the particle's current position as the origin, the above can then be converted into cartesian coordinates relative to its own position. These x,y, and z-directional steps are found after each time step for every particle, and summarized in an array. To convert these steps into cartesian coordinates, the steps are cumulatively added over time. 
Once this initial path is found, the boundary conditions are checked. If particle's are outside of the channel dimensions, then they are simply reflected across the boundaries until they are back inside the channel. Since step conditions are random, no difference was observed whether these checks occurred each time step, or all checks were done at the end. As such, boundary checks were done at the end of step generation to reduce the computational load.

Using the z-directional information, a corresponding velocity array is calculated according the the log-law equation. The data initially used in the model are typical values for mean longitudinal and shear velocity for the Fraser River. The velocity array is then multiplied by a unit time-step to find the distance traveled along the x-axis caused by the flow, and the particle location array is displaced accordingly. The coordinates in the pathing array are plotted in each time step, and Matplotlib's FuncAnimation class is used to animate these plots over time. 

If a staggered release is set to true, then a specified number of particles is repeatedly released after a certain number of time steps. As the initial step array is randomly generated, this staggering effect is simply achieved by writing an array of zeros to step values for each particle until they are released. The same is done for the velocity array, and the two staggered arrays are combined in the manner as above to generate the final coordinate array. For staggered release, each released particle group is drawn with a different color for ease of viewing. 

### Usage
Simply run the code to display the animation! As a warning, saving the animation will significantly increase the run time. Saved animations will not be exactly the same as ones shown in the console due to the nature of FFMpeg's writing process. The frame interval and repeat delay arguments in FuncAnimation will not affect the mp4 saved by FFMpeg. To change the frame interval, adjust the fps argument in the writer variable accordingly. To freeze the final position on screen, copy the final frame of the video according to the instructions here or simply pause the video. If the save function is not working, verify that it is installed in your system, then check the file path used. To change the file path, paste the following code snippet somewhere below the matplotlib import declarations: 

plt.rcParams['animation.ffmpeg_path'] =r'Replace this string with the application location for FFMpeg on your computer'

Hope you enjoy, feel free to suggeset any improvements! 

