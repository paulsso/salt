import os
from MatrixMethodTest import MatrixMethod, CreateGeometry
try:
    import numpy as np
    import matplotlib
    from tkinter import *
    import tkinter as tk
    from tkinter import ttk
    from matplotlib import pyplot as plt
    matplotlib.use("TkAgg")
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
    from mpl_toolkits.mplot3d import axes3d, Axes3D
except ModuleNotFoundError as e:
    print(e)
    os.system("pip install matplotlib")
    os.system("pip install numpy")
    from matplotlib import pyplot as plt
    import numpy as np
    from tkinter import *
    import tkinter as tk
    from tkinter import ttk
    import matplotlib
    matplotlib.use("TkAgg")
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
    from mpl_toolkits.mplot3d import axes3d, Axes3D

class masterframe(tk.Tk):

    def __init__(self, *args, **kwargs):

        def window1(self, frame):
            """Frame for top component being an array"""

            self.vdisp1_lab = Label(frame, text="Vertical displacement (mm)")
            self.vdisp1_lab.pack(pady=5)
            self.vdisp1_var = DoubleVar() # Vertical displacement
            self.vdisp1_scale = Scale(frame, from_=0, to=100, tickinterval=50,
            orient=HORIZONTAL, variable = self.vdisp1_var, length=250)
            self.vdisp1_scale.pack(padx=10, pady=2)

            self.roc1_lab = Label(frame, text="Radius of Curvature (mm)")
            self.roc1_lab.pack(pady=5)
            self.rcurv1_var = DoubleVar() # Radius of curvature
            self.rcurv1_scale = Scale(frame, from_=0, to=100, tickinterval=50,
            orient=HORIZONTAL, variable = self.rcurv1_var, length=250)
            self.rcurv1_scale.pack(padx=10, pady=2)

            self.layers1_lab = Label(frame, text="Number of layers")
            self.layers1_lab.pack(pady=5)
            self.layers1_var = DoubleVar() # Number of layers
            self.layers1_scale = Scale(frame, from_=1, to=8, tickinterval=1,
            orient=HORIZONTAL, variable = self.layers1_var, length=250)
            self.layers1_scale.pack(padx=10, pady=2)

            self.sock1_lab = Label(frame, text="Socket radius (mm)")
            self.sock1_lab.pack(pady=5)
            self.radius1_var = DoubleVar() # Disc radius
            self.radius1_scale = Scale(frame, from_=0, to=100, tickinterval=50,
            orient=HORIZONTAL, variable = self.radius1_var, length=250)
            self.radius1_scale.pack(padx=10, pady=2)

            self.phase1_lab = Label(frame, text = "Phase shift (degrees)")
            self.phase1_lab.pack(pady=5)
            self.phase1_var = DoubleVar()
            self.phase1_scale = Scale(frame, from_=0, to=179, tickinterval = 90,
            orient=HORIZONTAL, variable = self.phase1_var, length=250)
            self.phase1_scale.pack(padx=10, pady=2)

            self.transfreq1_lab = Label(frame, text="Transducer frequency (kHz):")
            self.transfreq1_lab.pack(pady=5)
            self.transfreq1_var = DoubleVar()
            self.transfreq1_ent = Entry(frame, width = 8, textvariable = self.transfreq1_var)
            self.transfreq1_ent.pack(padx=10, pady=2)

            self.depth1_lab = Label(frame, text="Number of transducers per layer")
            self.depth1_lab.pack(pady=5)
            self.depth1_var = [0]*8
            self.depth1_ent = [0]*8

            for i in range(0,8):
                self.depth1_var[i] = IntVar()
                self.depth1_ent[i] = Entry(frame, width=2, textvariable = self.depth1_var[i])
                self.depth1_ent[i].pack(side=LEFT, padx=10, pady=10)

        def window2(self, frame):
            """Frame for top component being a reflector"""
            self.depth1_var = [0]*8
            self.depth1_ent = [0]*8
            self.transfreq1_var = DoubleVar()
            self.transfreq1_ent = Entry(frame, width = 8, textvariable=self.transfreq1_var)

            for i in range(0,8):
                self.depth1_var[i] = IntVar()
                self.depth1_ent[i] = Entry(frame, width=2, textvariable = self.depth1_var[i])

            self.layers1_var = DoubleVar() # Number of layers
            self.layers1_scale = Scale(frame, from_=1, to=8, tickinterval=1,
            orient=HORIZONTAL, variable = self.layers1_var, length=250)

            self.vdisp1_lab = Label(frame, text="Vertical displacement (mm)")
            self.vdisp1_lab.pack(pady=5)
            self.vdisp1_var = DoubleVar() # Vertical displacement
            self.vdisp1_scale = Scale(frame, from_=0, to=100, tickinterval=50,
            orient=HORIZONTAL, variable = self.vdisp1_var, length=250)
            self.vdisp1_scale.pack(padx=10, pady=2)

            self.roc1_lab = Label(frame, text="Radius of curvature (mm)")
            self.roc1_lab.pack(pady=5)
            self.rcurv1_var = DoubleVar() # Radius of curvature
            self.rcurv1_scale = Scale(frame, from_=0, to=100, tickinterval=50,
            orient=HORIZONTAL, variable = self.rcurv1_var, length=250)
            self.rcurv1_scale.pack(padx=10, pady=2)

            self.sock1_lab = Label(frame, text="Socket radius (mm)")
            self.sock1_lab.pack(pady=5)
            self.radius1_var = DoubleVar() # Disc radius
            self.radius1_scale = Scale(frame, from_=0, to=100, tickinterval=50,
            orient=HORIZONTAL, variable = self.radius1_var, length=250)
            self.radius1_scale.pack(padx=10, pady=2)

        def window3(self, frame):
            """Frame for bottom component being an array"""

            self.vdisp2_lab = Label(frame, text="Vertical displacement (mm)")
            self.vdisp2_lab.pack(pady=5)
            self.vdisp2_var = DoubleVar() # Vertical displacement
            self.vdisp2_scale = Scale(frame, from_=0, to=100, tickinterval=50,
            orient=HORIZONTAL, variable = self.vdisp2_var, length=250)
            self.vdisp2_scale.pack(padx=10, pady=2)

            self.roc2_lab = Label(frame, text="Radius of curvature (mm)")
            self.roc2_lab.pack(pady=5)
            self.rcurv2_var = DoubleVar() # Radius of curvature
            self.rcurv2_scale = Scale(frame, from_=0, to=100, tickinterval=50,
            orient=HORIZONTAL, variable = self.rcurv2_var, length=250)
            self.rcurv2_scale.pack(padx=10, pady=2)

            self.layers2_lab = Label(frame, text="Number of layers")
            self.layers2_lab.pack(pady=5)
            self.layers2_var = DoubleVar() # Number of layers
            self.layers2_scale = Scale(frame, from_=1, to=8, tickinterval=1,
            orient=HORIZONTAL, variable = self.layers2_var, length=250)
            self.layers2_scale.pack(padx=10, pady=2)

            self.sock2_lab = Label(frame, text="Socket radius (mm)")
            self.sock2_lab.pack(pady=5)
            self.radius2_var = DoubleVar() # Disc radius
            self.radius2_scale = Scale(frame, from_=0, to=100, tickinterval=50,
            orient=HORIZONTAL, variable = self.radius2_var, length=250)
            self.radius2_scale.pack(padx=10, pady=2)

            self.phase2_lab = Label(frame, text = "Phase shift (degrees)")
            self.phase2_lab.pack()
            self.phase2_var = DoubleVar()
            self.phase2_scale = Scale(frame, from_=0, to=179, tickinterval = 90,
            orient=HORIZONTAL, variable = self.phase2_var, length=250)
            self.phase2_scale.pack(padx=10, pady=2)

            self.transfreq2_lab = Label(frame, text="Transducer frequency (kHz)")
            self.transfreq2_lab.pack(pady=5)
            self.transfreq2_var = DoubleVar()
            self.transfreq2_ent = Entry(frame, width = 8, textvariable = self.transfreq2_var)
            self.transfreq2_ent.pack(padx=10, pady=2)

            self.depth2_lab = Label(frame, text="Number of transducers per layer")
            self.depth2_lab.pack(pady=5)
            self.depth2_var = [0]*8
            self.depth2_ent = [0]*8

            for i in range(0,8):
                self.depth2_var[i] = IntVar()
                self.depth2_ent[i] = Entry(frame, width=2, textvariable = self.depth2_var[i])
                self.depth2_ent[i].pack(side=LEFT, padx=10, pady=10)

        def window4(self, frame):
            """Frame for bottom component being a reflector"""

            self.depth2_var = [0]*8
            self.depth2_ent = [0]*8
            self.transfreq2_var = DoubleVar()
            self.transfreq2_ent = Entry(frame, width = 8, textvariable = self.transfreq2_var)

            for i in range(0,8):
                self.depth2_var[i] = IntVar()
                self.depth2_ent[i] = Entry(frame, width=2, textvariable = self.depth2_var[i])

            self.layers2_var = DoubleVar() # Number of layers
            self.layers2_scale = Scale(frame, from_=1, to=8, tickinterval=1,
            orient=HORIZONTAL, variable = self.layers2_var, length=250)

            self.vdisp2_lab = Label(frame, text="Vertical displacement")
            self.vdisp2_lab.pack(pady=5)
            self.vdisp2_var = DoubleVar() # Vertical displacement
            self.vdisp2_scale = Scale(frame, from_=0, to=100, tickinterval=50,
            orient=HORIZONTAL, variable = self.vdisp2_var, length=250)
            self.vdisp2_scale.pack()

            self.roc2_lab = Label(frame, text="Radius of curvature")
            self.roc2_lab.pack(pady=5)
            self.rcurv2_var = DoubleVar() # Radius of curvature
            self.rcurv2_scale = Scale(frame, from_=0, to=100, tickinterval=50,
            orient=HORIZONTAL, variable = self.rcurv2_var, length=250)
            self.rcurv2_scale.pack()

            self.sock2_lab = Label(frame, text="Socket radius")
            self.sock2_lab.pack(pady=5)
            self.radius2_var = DoubleVar() # Disc radius
            self.radius2_scale = Scale(frame, from_=0, to=100, tickinterval=50,
            orient=HORIZONTAL, variable = self.radius2_var, length=250)
            self.radius2_scale.pack()

        def window5(self,frame):
            """Top single transducer"""

            self.vdisp1_lab = Label(frame, text="Vertical displacement (mm)")
            self.vdisp1_lab.pack(pady=5)
            self.vdisp1_var = DoubleVar() # Vertical displacement
            self.vdisp1_scale = Scale(frame, from_=0, to=100, tickinterval=50,
            orient=HORIZONTAL, variable = self.vdisp1_var, length=250)
            self.vdisp1_scale.pack(padx=10, pady=2)

            self.roc1_lab = Label(frame, text="Radius of Curvature (mm)")
            self.roc1_lab.pack(pady=5)
            self.rcurv1_var = DoubleVar() # Radius of curvature
            self.rcurv1_scale = Scale(frame, from_=0, to=100, tickinterval=50,
            orient=HORIZONTAL, variable = self.rcurv1_var, length=250)
            self.rcurv1_scale.pack(padx=10, pady=2)

            self.sock1_lab = Label(frame, text="Socket radius (mm)")
            self.sock1_lab.pack(pady=5)
            self.radius1_var = DoubleVar() # Disc radius
            self.radius1_scale = Scale(frame, from_=0, to=100, tickinterval=50,
            orient=HORIZONTAL, variable = self.radius1_var, length=250)
            self.radius1_scale.pack(padx=10, pady=2)

            self.phase1_lab = Label(frame, text = "Phase shift (degrees)")
            self.phase1_lab.pack(pady=5)
            self.phase1_var = DoubleVar()
            self.phase1_scale = Scale(frame, from_=0, to=179, tickinterval = 90,
            orient=HORIZONTAL, variable = self.phase1_var, length=250)
            self.phase1_scale.pack(padx=10, pady=2)

            self.transfreq1_lab = Label(frame, text="Transducer frequency (kHz):")
            self.transfreq1_lab.pack(pady=5)
            self.transfreq1_var = DoubleVar()
            self.transfreq1_ent = Entry(frame, width = 8, textvariable = self.transfreq1_var)
            self.transfreq1_ent.pack(padx=10, pady=2)

        def show1(none=0):
            """function which destroys old frame widgets and calls function to open new frame"""
            if hasattr(self, 'vdisp1_scale'):
                self.vdisp1_scale.destroy()
                self.rcurv1_scale.destroy()
                self.vdisp1_lab.destroy()
                self.roc1_lab.destroy()

            if hasattr(self, 'layers1_scale'):
                self.layers1_scale.destroy()
            if hasattr(self, 'layers1_lab'):
                self.layers1_lab.destroy()

            if hasattr(self, 'radius1_scale'):
                self.radius1_scale.destroy()
            if hasattr(self, 'sock1_lab'):
                self.sock1_lab.destroy()

            if hasattr(self, 'depth1_lab'):
                self.depth1_lab.destroy()
                for i in range(0,8):
                    self.depth1_ent[i].destroy()

            if hasattr(self, 'transfreq1_ent'):
                self.transfreq1_ent.destroy()
            if hasattr(self, 'transfreq1_lab'):
                self.transfreq1_lab.destroy()

            if hasattr(self, 'phase1_lab'):
                self.phase1_lab.destroy()
            if hasattr(self, 'phase1_scale'):
                self.phase1_scale.destroy()

            checkThisString = self.clicked1.get()
            if checkThisString == "Array":
                window1(self,frame)
            elif checkThisString == "Reflector":
                window2(self,frame)
            elif checkThisString == "Transducer":
                window5(self,frame)

        def show2(none=0):
            """function which destroys old frame widgets and calls function to open new frame"""

            if hasattr(self, 'vdisp2_scale'):
                self.vdisp2_scale.destroy()
                self.rcurv2_scale.destroy()
                self.vdisp2_lab.destroy()
                self.roc2_lab.destroy()

            if hasattr(self, 'layers2_scale'):
                self.layers2_scale.destroy()
            if hasattr(self, 'layers2_lab'):
                self.layers2_lab.destroy()

            if hasattr(self, 'radius2_scale'):
                self.radius2_scale.destroy()
            if hasattr(self, 'sock2_lab'):
                self.sock2_lab.destroy()

            if hasattr(self, 'depth2_lab'):
                self.depth2_lab.destroy()
                for i in range(0,8):
                    self.depth2_ent[i].destroy()

            if hasattr(self, 'transfreq2_ent'):
                self.transfreq2_ent.destroy()
            if hasattr(self, 'transfreq2_lab'):
                self.transfreq2_lab.destroy()

            if hasattr(self, 'phase2_lab'):
                self.phase2_lab.destroy()
            if hasattr(self, 'phase2_scale'):
                self.phase2_scale.destroy()

            checkThisString = self.clicked2.get()
            if checkThisString == "Array":
                window3(self,frame2)
            elif checkThisString == "Reflector":
                window4(self,frame2)

        def render(self):
            """Function organizes input from widgets and plots the geometry"""
            def properties_fix(properties):
                if properties["Type"] == "Reflector":
                    properties["Layers"] = 0
                    properties["Depth"] = [0, 0, 0, 0, 0, 0, 0, 0]
                    properties["TransFreq"] = 0
                    properties["Amplitude"] = 0
                    properties["TransRadius"] = 0

                properties_fixed = properties
                return properties_fixed

            if self.clicked1.get() == "Array":
                zPosProperties = {
                                "Orientation": 1,
                                "Type": self.clicked1.get(),
                                "Phase": self.phase1_scale.get(),
                                "RadiusCurvature": self.rcurv1_scale.get()*1e-3,
                                "Layers": self.layers1_scale.get(),
                                "Depth": [int(self.depth1_ent[0].get()), int(self.depth1_ent[1].get()), int(self.depth1_ent[2].get()), int(self.depth1_ent[3].get()), int(self.depth1_ent[4].get()), int(self.depth1_ent[5].get()),
                                int(self.depth1_ent[6].get()), int(self.depth1_ent[7].get())],
                                "Radius": self.radius1_scale.get()*1e-3,
                                "TransFreq": float(self.transfreq1_ent.get())*1e3,
                                "Amplitude": 1e-6,
                                "Displacement": self.vdisp1_scale.get()*1e-3,
                                "TransRadius": 4.5e-3
                }
            elif self.clicked1.get() == "Transducer":
                zPosProperties = {
                                "Orientation": 1,
                                "Type": self.clicked1.get(),
                                "Phase": self.phase1_scale.get(),
                                "RadiusCurvature": self.rcurv1_scale.get()*1e-3,
                                "Layers": 0,
                                "Depth": [0,0,0,0,0,0,0,0],
                                "Radius": self.radius1_scale.get()*1e-3,
                                "TransFreq": float(self.transfreq1_ent.get())*1e3,
                                "Amplitude": 1e-6,
                                "Displacement": self.vdisp1_scale.get()*1e-3,
                                "TransRadius": 4.5e-3
                }
            if self.clicked2.get() == "Array":
                zNegProperties = {
                                "Orientation": -1,
                                "Type": self.clicked2.get(),
                                "Phase": self.phase2_scale.get(),
                                "RadiusCurvature": self.rcurv2_scale.get()*1e-3,
                                "Layers": self.layers2_scale.get(),
                                "Depth": [int(self.depth2_ent[0].get()), int(self.depth2_ent[1].get()), int(self.depth2_ent[2].get()), int(self.depth2_ent[3].get()), int(self.depth2_ent[4].get()), int(self.depth2_ent[5].get()),
                                int(self.depth2_ent[6].get()), int(self.depth2_ent[7].get())],
                                "Radius": self.radius2_scale.get()*1e-3,
                                "TransFreq": float(self.transfreq2_ent.get())*1e3,
                                "Amplitude": 1e-6,
                                "Displacement": self.vdisp2_scale.get()*1e-3,
                                "TransRadius": 4.5e-3
                }
            elif self.clicked2.get() == "Reflector":
                zNegProperties = {
                                "Orientation": -1,
                                "Type": self.clicked2.get(),
                                "Phase": 0,
                                "RadiusCurvature": self.rcurv2_scale.get()*1e-3,
                                "Layers": 0,
                                "Depth": [0,0,0,0,0,0,0,0],
                                "Radius": self.radius2_scale.get()*1e-3,
                                "TransFreq": float(self.transfreq2_ent.get())*1e3,
                                "Amplitude": 1e-6,
                                "Displacement": self.vdisp2_scale.get()*1e-3,
                                "TransRadius": 4.5e-3
                }

            Ux, Uy, Uz, Vx, Vy, Vz = CreateGeometry(zPosProperties, zNegProperties)

            fig = plt.figure()
            ax = Axes3D(fig)
            ax.scatter3D(Vx, Vy, Vz, marker='.')
            ax.scatter3D(Ux, Uy, Uz, marker='.')
            ax.set_xlim([-0.05,0.05])
            ax.set_ylim([-0.05,0.05])
            ax.set_zlim([-0.065,0.065])
            plt.show()

        def compute(self):
            """Calls function to run computations and plots result"""

            mediumProperties = {"Density": 1.2, "SpeedOfSound": 343}
            if self.clicked1.get() == "Array":
                zPosProperties = {
                                "Orientation": 1,
                                "Type": self.clicked1.get(),
                                "Phase": self.phase1_scale.get(),
                                "RadiusCurvature": self.rcurv1_scale.get()*1e-3,
                                "Layers": self.layers1_scale.get(),
                                "Depth": [int(self.depth1_ent[0].get()), int(self.depth1_ent[1].get()), int(self.depth1_ent[2].get()), int(self.depth1_ent[3].get()), int(self.depth1_ent[4].get()), int(self.depth1_ent[5].get()),
                                int(self.depth1_ent[6].get()), int(self.depth1_ent[7].get())],
                                "Radius": self.radius1_scale.get()*1e-3,
                                "TransFreq": float(self.transfreq1_ent.get())*1e3,
                                "Amplitude": 1e-6,
                                "Displacement": self.vdisp1_scale.get()*1e-3,
                                "TransRadius": 4.5e-3
                }
            elif self.clicked1.get() == "Transducer":
                zPosProperties = {
                                "Orientation": 1,
                                "Type": self.clicked1.get(),
                                "Phase": self.phase1_scale.get(),
                                "RadiusCurvature": self.rcurv1_scale.get()*1e-3,
                                "Layers": 0,
                                "Depth": [0,0,0,0,0,0,0,0],
                                "Radius": self.radius1_scale.get()*1e-3,
                                "TransFreq": float(self.transfreq1_ent.get())*1e3,
                                "Amplitude": 1e-6,
                                "Displacement": self.vdisp1_scale.get()*1e-3,
                                "TransRadius": 4.5e-3
                }
            if self.clicked2.get() == "Array":
                zNegProperties = {
                            "Orientation": -1,
                            "Type": self.clicked2.get(),
                            "Phase": self.phase2_scale.get(),
                            "RadiusCurvature": self.rcurv2_scale.get()*1e-3,
                            "Layers": self.layers2_scale.get(),
                            "Depth": [int(self.depth2_ent[0].get()), int(self.depth2_ent[1].get()), int(self.depth2_ent[2].get()), int(self.depth2_ent[3].get()), int(self.depth2_ent[4].get()), int(self.depth2_ent[5].get()),
                            int(self.depth2_ent[6].get()), int(self.depth2_ent[7].get())],
                            "Radius": self.radius2_scale.get()*1e-3,
                            "TransFreq": float(self.transfreq2_ent.get())*1e3,
                            "Amplitude": 1e-6,
                            "Displacement": self.vdisp2_scale.get()*1e-3,
                            "TransRadius": 4.5e-3}
            elif self.clicked2.get() == "Reflector":
                zNegProperties = {
                            "Orientation": -1,
                            "Type": self.clicked2.get(),
                            "Phase": 0,
                            "RadiusCurvature": self.rcurv2_scale.get()*1e-3,
                            "Layers": 0,
                            "Depth": [0,0,0,0,0,0,0,0],
                            "Radius": self.radius2_scale.get()*1e-3,
                            "TransFreq": float(self.transfreq2_ent.get())*1e3,
                            "Amplitude": 1e-6,
                            "Displacement": self.vdisp2_scale.get()*1e-3,
                            "TransRadius": 4.5e-3}

            start = time.time()
            relative_potential, pressure, x_span, z_span = MatrixMethod(mediumProperties,zPosProperties,zNegProperties)
            end = time.time()
            diff = end - start
            print("Total time elapsed was %.4f" % diff, "seconds")
            x = len(x_span)
            z = len(z_span)
            _relative_potential = np.real(relative_potential).reshape([z, x])
            _pressure = np.real(pressure).reshape([z, x])
            maxPotential = np.max(_relative_potential)
            minPotential = np.min(_relative_potential)
            maxPressure = np.max(_relative_potential)
            minPressure = np.min(_relative_potential)
            fig=plt.figure(figsize=(8, 6), dpi= 80, facecolor='w', edgecolor='k')
            cs = plt.contourf(_relative_potential, np.linspace(minPotential,maxPotential,1000), cmap='jet', extend='both')
            plt.colorbar()
            plt.show()

            fig=plt.figure(figsize=(8, 6), dpi= 80, facecolor='w', edgecolor='k')
            cs = plt.contourf(_pressure, np.linspace(minPressure,maxPressure,1000), cmap='jet', extend='both')
            plt.colorbar()
            plt.show()


        tk.Tk.__init__(self, *args, **kwargs)
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.geometry("790x800")
        self.iconbitmap(".\icon.ico")
        self.title("Simulation Platform for Acoustic Levitation Traps (SALT)")
        self.frames = {}

        self.clicked1 = StringVar()
        self.clicked1.set("Window 1")
        frame = LabelFrame(self, text="Top component", padx=10, pady=10)
        frame.pack(padx=10, pady=10)
        frame.place(x=50, y=0)
        drop = OptionMenu(frame, self.clicked1, "Array", "Transducer", command=show1)
        drop.pack()

        self.clicked2 = StringVar()
        self.clicked2.set("Window 2")
        frame2 = LabelFrame(self, text="Bottom components", padx=10, pady=10)
        frame2.pack(padx=10, pady=10)
        frame2.place(x=420, y=0)
        drop2 = OptionMenu(frame2, self.clicked2, "Array", "Reflector", command=show2)
        drop2.pack()

        btn1 = Button(self, text="Display Geometry", command=lambda:render(self))
        btn1.pack(side=BOTTOM,pady=10)
        btn2 = Button(self, text="Compute Relative Acoustic Potential", command=lambda:compute(self))
        btn2.pack(side=BOTTOM,pady=10)

root = masterframe()
root.mainloop()
