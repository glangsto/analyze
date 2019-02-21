import numpy as np
from numba import jit
from matplotlib import pyplot as plt
from matplotlib import colors 
#from matplotlib import inline

#matplotlib inline

def mandelbrot(c,maxiter):
    """
    Low level function to compute the mandelbrot set for a single complex location
    """
    z = c
    for n in range(maxiter):
        if abs(z) > 2:
            return n
        z = z*z + c
    return 0


def mandelbrot_set(xmin,xmax,ymin,ymax,width,height,maxiter):
    """
    compute the mandelbrot set for a whole region
    """
    r1 = np.linspace(xmin, xmax, width)
    r2 = np.linspace(ymin, ymax, height)
    n3 = np.empty((width,height))
    for i in range(width):
        for j in range(height):
            n3[i,j] = mandelbrot(r1[i] + 1j*r2[j],maxiter)
        if int(i/50)*50 == i:  # print every so often
            print "Row: ", i
    return (r1,r2,n3)

image_counter = 32

def save_image(fig):
    global image_counter
    filename = "mandelbrodt_%d.png" % image_counter
    image_counter += 1
    fig.savefig(filename)

def mandelbrot_image(xmin,xmax,ymin,ymax,width=15,height=15,maxiter=100):
    dpi = 72
    img_width = dpi * width
    img_height = dpi * height
    x,y,z = mandelbrot_set(xmin,xmax,ymin,ymax,img_width,img_height,maxiter)
    
    fig, ax = plt.subplots(figsize=(width, height),dpi=72)
    ticks = np.arange(0,img_width,3*dpi)
    x_ticks = xmin + (xmax-xmin)*ticks/img_width
    plt.xticks(ticks, x_ticks)
    y_ticks = ymin + (ymax-ymin)*ticks/img_width
    plt.yticks(ticks, y_ticks)
    
    ax.imshow(z.T,origin='lower') 
    
    save_image(fig)

def mandelbrot_ii(xmin,xmax,ymin,ymax,width=10,height=10,maxiter=80,cmap='jet'):
    dpi = 72
    img_width = dpi * width
    img_height = dpi * height
    x,y,z = mandelbrot_set(xmin,xmax,ymin,ymax,img_width,img_height,maxiter)
    
    fig, ax = plt.subplots(figsize=(width, height),dpi=72)
    ticks = np.arange(0,img_width,3*dpi)
    x_ticks = xmin + (xmax-xmin)*ticks/img_width
    plt.xticks(ticks, x_ticks)
    y_ticks = ymin + (ymax-ymin)*ticks/img_width
    plt.yticks(ticks, y_ticks)
    ax.set_title(cmap)
    
    ax.imshow(z.T,cmap=cmap,origin='lower') 
    
    save_image(fig)

#mandelbrot_image(-2.0,0.5,-1.25,1.25)
#mandelbrot_image(-1.8,-1.0,0.0,.8)
mandelbrot_image(-1.3,-1.15,0.1,.25)
