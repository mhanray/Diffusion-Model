import numpy as np 
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import axes3d as p3
from mpl_toolkits.mplot3d import art3d as art3d
import sys

def animate_scatters(iteration, data, scatters, data_ref, scatter_ref):
    for i in range(particle_n):
        scatters[i]._offsets3d=(data[iteration,0:1,i], data[iteration,1:2,i], data[iteration,2:,i])
    scatter_ref[0]._offsets3d=(data_ref[iteration,0:1,0],data_ref[iteration,1:2,0],data_ref[iteration,2:,0])
    
    return [scatters,scatter_ref]

def generate_steps():
    origin=np.zeros((1,3,particle_n))
    origin[:,1,:]+=y_initial-width/2
    origin[:,2,:]+=z_initial-depth/2
    
    steps=np.zeros((step_n,3,particle_n))
    rho=np.random.normal(loc=0,scale=sd,size=(step_n,particle_n))
    theta=np.random.uniform(size=(step_n,particle_n))*2*np.pi
    phi=np.random.uniform(size=(step_n,particle_n))*2*np.pi
    steps[:,0,:]=rho*np.sin(phi)*np.cos(theta)
    steps[:,1,:]=rho*np.sin(phi)*np.sin(theta)
    steps[:,2,:]=rho*np.cos(phi)
    
    if stagger_release==False:
        path=np.concatenate((origin,steps),0).cumsum(0)
        
    else:
        steps=stagger_path(steps)
        path=np.concatenate((origin,steps),0).cumsum(0)
    
    return path

def generate_displacement(path):
    if stagger_release==False:
        displacement=((u_mean+u_shear/vk*(np.log((path[1:,2,:]+depth/2)/depth)+1))*dt).cumsum(0)
        
    else:
        velocity=((u_mean+u_shear/vk*(np.log((path[1:,2,:]+depth/2)/depth)+1))*dt)
        displacement=stagger_velocity(velocity).cumsum(0)
    return displacement

def check_boundaries(path):
    for k in range(particle_n):
        for i in range(step_n+1):
            while np.abs(path[i,1,k])>width/2:
                if path[i,1,k]>width/2:
                    path[i,1,k]-= (path[i,1,k]-width/2)*2
                if path[i,1,k]<-width/2:
                    path[i,1,k]-= (path[i,1,k]+width/2)*2 
                    
            while np.abs(path[i,2,k])>depth/2:     
                if path[i,2,k]>depth/2:
                    path[i,2,k]-= (path[i,2,k]-depth/2)*2
                if path[i,2,k]<-depth/2:
                    path[i,2,k]-= (path[i,2,k]+depth/2)*2 
    return path

def stagger_path(data):
    staggered=data.copy()
    if dt<1: step_delay=int(release_delay/dt)
    else: step_delay=release_delay
    
    for i in range(int(step_n/step_delay)):
        for count in range(step_delay*(i), step_delay*(i+1)):
            for k in range(int(release_n*(i+1)),int(particle_n)):
                staggered[count,:,k]=0
    return staggered

def stagger_velocity(velocity):
    staggered=velocity.copy()
    if dt<1: step_delay=int(release_delay/dt)
    else: step_delay=release_delay
    for count in range(int(step_n/step_delay)):
        for i in range(step_delay*(count), step_delay*(count+1)):
            for k in range(int(release_n*(count+1)),int(particle_n)):
                staggered[i,k]=0
    return staggered

def generate_reference():
    reference=np.zeros((step_n+1,3,1))
    reference[:,1,:]+=y_initial-width/2
    reference[:,2,:]+=depth/2
    
    u_origin=u_mean+u_shear/vk*(np.log((depth)/depth)+1)
    for i in range(step_n+1):
        reference[i,0,0]+=u_origin*dt*(i)
    return reference

def generate_streamline(data):
    streamline=np.empty((step_n+1,3))
    for i in range(step_n+1):
        streamline[i,0]=(u_mean+u_shear/vk*(np.log((z_initial+depth/2)/depth)+1))*dt*i
        streamline[i,1]=np.average(data[i,1,:])
        streamline[i,2]=np.average(data[i,2,:])
    return streamline

def draw_water(length,width,depth,color):
    water_surface=Rectangle((0,-width/2),length,width,alpha=0.1,color=color)
    water_bottom=Rectangle((0,-width/2),length,width,alpha=0.1,color=color)
    wall_left=Rectangle((0,-depth/2),length,depth,alpha=0.1,color=color)
    wall_right=Rectangle((0,-depth/2),length,depth,alpha=0.1,color=color)
    upstream=Rectangle((-width/2,-depth/2),width,depth,alpha=0.1,color=color)
    downstream=Rectangle((-width/2,-depth/2),width,depth,alpha=0.1,color=color)
    ax.add_patch(water_surface)
    ax.add_patch(water_bottom)
    ax.add_patch(wall_left)
    ax.add_patch(wall_right)
    ax.add_patch(upstream)
    ax.add_patch(downstream)
    art3d.pathpatch_2d_to_3d(water_surface, z=depth/2, zdir="z")
    art3d.pathpatch_2d_to_3d(water_bottom, z=-depth/2, zdir="z")
    art3d.pathpatch_2d_to_3d(wall_left, z=-width/2, zdir="y")
    art3d.pathpatch_2d_to_3d(wall_right, z=width/2, zdir="y")
    art3d.pathpatch_2d_to_3d(upstream, z=-0, zdir="x")
    art3d.pathpatch_2d_to_3d(downstream, z=length, zdir="x")

def check_model():
    if step_n<=0 or diff<=0 or dt<0 or particle_n<=0 or vk<=0 or depth<=0 or width<=0 or z_initial>depth or z_initial<0 or y_initial<0 or y_initial>width:
        print('Please check model parameters.')
        sys.exit()
    elif time%dt!=0:
        print('Simulation time must be a multiple of the time step.')
        sys.exit() 
    elif stagger_release==True and (step_n)%(release_delay)!= 0:
        print('Total number of steps must be a multiple of the release delay.')
        sys.exit()
    else: return
    
#set up model parameters 
time=10     #seconds, length of time to simulate
dt=0.25   #seconds, time step
diff=7**(-3)    #m^3/s, diffusive coefficient
particle_n=200   #total number of particles in simulation
step_n=int(time/dt)  #total number of time steps 

u_mean=1.4     #m/s, mean longitudinal velocity 
u_shear=0.1     #m/s, shear velocity
vk=0.41     #von Karman constant
sd=np.sqrt(2*diff*dt)   #standard deviation of Gaussian distribution

width=2     #m, width in y-axis
depth=2     #m, depth in z-axis
y_initial=1     #m, initial y location
z_initial=0.2    #m, initial z location

stagger_release=False    #Repeatedly release set number of particles after a specified delay
release_delay=2    #seconds, delay period 
release_n=50     #number of particles to be released 

save=False   #save animation using FFmpeg   
streamline_visible=False    #draw streamline when plotting 
water_visible=True    #draw water outline in plot

#check model parameters
check_model()

#generate particle origin with channel centre of channel at y=0,z=0, generate step array by sampling from Gaussian distribution, then combine as particle coordinates
path=generate_steps()

#check particle coordinates, reflect back in if out of bounds
path=check_boundaries(path)

#generate velocity array, then displacement due to velocity.
displacement=generate_displacement(path)

#transform coordinates using displacement array
for k in range(particle_n):
    for i in range(step_n):
        path[i+1,0,k]+=displacement[i,k]
        
#add a reference particle traveling on water surface
reference=generate_reference()
if streamline_visible==True:streamline=generate_streamline(path)

#set up figure 
fig = plt.figure(figsize=(10,5))
ax = p3.Axes3D(fig)

if u_mean==0: ax.set_xlim3d([-diff*step_n*2,diff*step_n*2])
else: ax.set_xlim3d([-1,step_n*u_mean*dt*1.5])
ax.set_ylim3d([-width/2-1,width/2+1])
ax.set_zlim3d([-depth/2-1,depth/2+1])


ax.set_xlabel('X Distance (m)')
ax.set_ylabel('Y Distance (m)')
ax.set_zlabel('Z Distance (m)')

if stagger_release==False: plot_title= '3D Diffusion of %d Particles Over %.2f Seconds' % (particle_n, time)
else: plot_title= '3D Diffusion of %d Particles Released Every %.1f Seconds over %.1f Seconds' % (release_n, release_delay, time)
plt.suptitle(plot_title,fontsize=12)
plt.title('Mean Longitudinal Flow Velocity: %.2f m/s' % (u_mean),y=0.96, fontsize=10)

#visual components, draw water if selected and vary particle colors
if stagger_release==False: release_n=particle_n
if water_visible==True and u_mean>0: draw_water(step_n*u_mean*dt*1.5,width,depth,'cyan')
if streamline_visible==True: ax.plot(streamline[:,0], streamline[:,1], streamline[:,2], color='red')
coloring=cm.winter(np.linspace(0,1,int(particle_n/release_n)))

#initialize scatter plots
scatters=[ax.scatter(path[0,0:1,i], path[0,1:2,i], path[0,2:,i], color=coloring[int(i/(release_n))], alpha=0.5,edgecolors='face',marker=u'o') for i in range(particle_n)]    
scatter_ref=[ax.scatter(reference[0,0:1,0], reference[0,1:2,0], reference[0,2:,0], color='magenta', alpha=1,edgecolors='face',marker='o', s=50)]  
iterations=step_n+1


#animate and save video as mp4 if selected
anim=animation.FuncAnimation(fig,animate_scatters, iterations,fargs=(path, scatters, reference, scatter_ref), interval= 250, blit=False, repeat=True, repeat_delay=5000)
plt.show()
if save==True:
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, extra_args=['-vcodec', 'libx264'])
    anim.save('3D_Diffusion.mp4', writer=writer)
