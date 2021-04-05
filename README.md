# Simulation platform for Acoustic Levitation Traps (SALT)
#### A configurable implementation of a matrix method to approximate the pressure field in acoustic levitation devices

## How to use SALT:

### Installation
#### 1. Install python
SALT was tested in a virtual environment using python 3.7.4, its only dependencies are ```matplotlib``` and ```numpy```.

#### 2. Clone this repo
Use github desktop and copy-paste the URL of this repo or just download the zip-file manually.

#### 3. Launch the application
Launch the application by running the file ```main.py```, simply by double clicking the file or navigating to the folder where it is located in and typing ```python main.py``` in your command prompt or terminal.

### Usage

**Choose components:** The top component can be either an array or a transducer while the bottom component can be either an array or a reflector.
The difference between arrays/transducers and reflectors is that they do not simulate excitation of the acoustic field.  

<a href="https://ibb.co/Mn6ZfW4"><img src="https://i.ibb.co/PxQFM3n/component.png" alt="component" border="0"></a>

**Convaity:** Checked by default, if uncheck results in a flat array/transducer/reflector and thus the radius of curvature slider has no effect.

<a href="https://ibb.co/xhJDjbG"><img src="https://i.ibb.co/1Z8LqNK/concavity.png" alt="concavity" border="0"></a>

**Vertical offset:** Determines the distance between the center of the array/transducer/reflector from the origin. E.g. the total distance between two arrays with this value set to 64mm on both sides is 128mm.

<a href="https://ibb.co/Yj4R8Jw"><img src="https://i.ibb.co/C0p1v4D/vertical.png" alt="vertical" border="0"></a>

**Radius of curvature:** The radius of a sphere which the surface of the transducers will lie on.

<a href="https://ibb.co/y58K8Jt"><img src="https://i.ibb.co/Xpz6zT1/roc.png" alt="roc" border="0"></a>

**Maximum size of the socket:** The size of the socket determines how many rows of transducers can be put in an array, experimenting with the relation between this and radius of curvature is encouraged to understand how to define a geometry that works. 

<a href="https://ibb.co/B4V9PKL"><img src="https://i.ibb.co/XyZq4kX/socket.png" alt="socket" border="0"></a>

**Phase shift:** The phase component of the wave.

<a href="https://ibb.co/7XjtdXQ"><img src="https://i.ibb.co/t4QpT4z/phase.png" alt="phase" border="0"></a>

**Adjust frequency:** Fine tune the frequency at which you wish to excite the acoustic field.

<a href="https://ibb.co/09WVnCp"><img src="https://i.ibb.co/n1V8w7q/frequency.png" alt="frequency" border="0"></a>

**Configure array:** Almost any combination of number of rows and density on each row is available, experiment!

<a href="https://ibb.co/Wx03Tzj"><img src="https://i.ibb.co/ZLSJZg7/transducers.png" alt="transducers" border="0"></a>

**Options:** To familiarize oneself with how to configure the acoustic trap it is recommended to use the "display geometry" option, before running the simulation. When the simulation is run, the terminal window from where the GUI was launched will print messages on the progress of the simulation. When the simulation is done two plots will appear, one is a contour/surface plot of the acoustic radiation pressure, the other is the profile of the acoustic radiation pressure. 

<a href="https://ibb.co/X3q0hFb"><img src="https://i.ibb.co/hZSQP2y/options.png" alt="options" border="0"></a>

When the simulation is done three text files will be created containing comma-seperated values, one will contain the points on the contour-plot called acoustic_radiation.txt, one contains the values on the x-axis of the profile of the radiation pressure (xPlot.txt) and one profile.txt contains the values on the y-axis of the acoustic radiation profile.
 
#### 
## What's left?
Among many things that can be added are:

Convex shape arrays
Square/rectangular arrays with novel arrangement patterns
Mathematics coded in C for faster execution, better scaleability.
