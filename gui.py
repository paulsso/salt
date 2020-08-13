import tkinter as tk
import time, sys, os
try:
    import numpy as np
except:
    os.system("pip install numpy")
    import numpy as np

try:
    from matplotlib import pyplot as plt
except:
    os.system("pip install matplotlib")
    from matplotlib import pyplot as plt

from MatrixMethod import MatrixMethod, CreateGeometry
from mpl_toolkits.mplot3d import axes3d, Axes3D

def definedicts(page_inst, masterframe_inst):
    c = float(masterframe_inst.medium_c.get())
    rho = float(masterframe_inst.medium_density.get())

    medium_properties = {"Density": rho, "SpeedOfSound": c}

    top_properties = {"Orientation":1,
    "Type": masterframe_inst.top_selection.get(),
    "Concave":page_inst.properties_page1["Concave"].get()}

    bot_properties = {"Orientation":-1,
    "Type":masterframe_inst.bot_selection.get(),
    "Concave":page_inst.properties_page2["Concave"].get()}

    for key in page_inst.properties_page1:
        if key != "Depth":
            top_properties.update({key : page_inst.properties_page1[key].get()})
            if key in ["Displacement", "RadiusCurvature", "Radius"]:
                top_properties[key] = float(top_properties[key])*1e-3
            if key == "TransFreq":
                top_properties[key] = float(top_properties[key])*1e3
            if key == "Layers":
                top_properties[key] = int(top_properties[key])

        else:
            v = []
            for item in page_inst.properties_page1[key]:
                v = np.append(v, int(item.get()))
            top_properties.update({key : v})

    for key in page_inst.properties_page2:
        if key != "Depth":
            bot_properties.update({key : page_inst.properties_page2[key].get()})
            if key in ["Displacement", "RadiusCurvature", "Radius"]:
                bot_properties[key] = float(bot_properties[key])*1e-3
            if key == "TransFreq":
                bot_properties[key] = float(bot_properties[key])*1e3
            if key == "Layers":
                bot_properties[key] = int(bot_properties[key])
        else:
            v = []
            for item in page_inst.properties_page2[key]:
                v = np.append(v, int(item.get()))
            bot_properties.update({key : v})

    return bot_properties, top_properties, medium_properties

class types:
    def array(self, frame):
        """ Method adds widgets relevant to an array of transducers """
        labels, scales, entries = [], [], []
        properties = {}

        var_checkbox = tk.BooleanVar()
        checkbox = tk.Checkbutton(frame, text="Concave", variable=var_checkbox, onvalue=True, offvalue=False)
        checkbox.pack(pady=10)

        properties.update({"Concave" : var_checkbox})

        labels = np.append(labels, tk.Label(frame, text="Vertical Displacement (mm)"))
        labels[0].pack(pady=0)
        scales = np.append(scales, tk.Scale(frame, from_=0, to=100, tickinterval=50,
        orient=tk.HORIZONTAL, variable = tk.DoubleVar(), length=250))
        scales[0].pack(padx=10, pady=2)
        scales[0].set(68)

        properties.update({"Displacement" : scales[0]})

        labels = np.append(labels, tk.Label(frame, text="Radius of Curvature (mm)"))
        labels[1].pack(pady=0)
        scales = np.append(scales, tk.Scale(frame, from_=0, to=100, tickinterval=50,
        orient=tk.HORIZONTAL, variable = tk.DoubleVar(), length=250))
        scales[1].pack(padx=10, pady=2)
        scales[1].set(48)

        properties.update({"RadiusCurvature" : scales[1]})

        labels = np.append(labels, tk.Label(frame, text="Socket Radius (mm)"))
        labels[2].pack(pady=0)
        scales = np.append(scales, tk.Scale(frame, from_=0, to=100, tickinterval=50,
        orient=tk.HORIZONTAL, variable = tk.DoubleVar(), length=250))
        scales[2].pack(padx=10, pady=2)
        scales[2].set(50)

        properties.update({"Radius" : scales[2]})

        labels = np.append(labels, tk.Label(frame, text="Phase shift (degrees)"))
        labels[3].pack(pady=0)
        scales = np.append(scales, tk.Scale(frame, from_=0, to=179, tickinterval=90,
        orient=tk.HORIZONTAL, variable = tk.DoubleVar(), length=250))
        scales[3].pack(padx=10, pady=2)

        properties.update({"Phase" : scales[3]})

        labels = np.append(labels, tk.Label(frame, text="Transducer frequency (kHz):"))
        labels[4].pack(pady=0)
        entries = np.append(entries, tk.Entry(frame, width = 8))
        entries[0].pack(padx=10, pady=2)
        entries[0].insert(0, 40.0)

        properties.update({"TransFreq" : entries[0]})

        labels = np.append(labels, tk.Label(frame, text="Number of Layers:"))
        labels[5].pack(pady=0)
        scales = np.append(scales, tk.Scale(frame, from_=1, to=8, tickinterval=1,
        orient=tk.HORIZONTAL, variable = tk.IntVar(), length=250))
        scales[4].pack(padx=10, pady=2)
        scales[4].set(3)

        properties.update({"Layers" : scales[4]})

        labels = np.append(labels, tk.Label(frame, text = "Number of transducer per layer:"))
        labels[6].pack(pady=0)
        v = [6, 12, 18, 24, 32, 36, 42, 48]
        for i in range(0,8):
            entries = np.append(entries, tk.Entry(frame, width=2))
            entries[1+i].pack(side=tk.LEFT, padx=10, pady=10)
            entries[1+i].insert(0,v[i])

        properties.update({"Depth" : entries[1:]})

        return properties, labels, scales, entries, checkbox

    def transducer(self, frame):
        """ Method adds widgets relevant to a transducer """
        labels, scales, entries = [], [], []
        properties = {}

        var_checkbox = tk.BooleanVar()
        checkbox = tk.Checkbutton(frame, text="Concave", variable=var_checkbox, onvalue=True, offvalue=False)
        checkbox.pack(pady=10)

        properties.update({"Concave" : var_checkbox})

        labels = np.append(labels, tk.Label(frame, text="Vertical Displacement (mm)"))
        labels[0].pack(pady=0)
        scales = np.append(scales, tk.Scale(frame, from_=0, to=100, tickinterval=50,
        orient=tk.HORIZONTAL, variable = tk.DoubleVar(), length=250))
        scales[0].pack(padx=10, pady=2)
        scales[0].set(68)

        properties.update({"Displacement" : scales[0]})

        labels = np.append(labels, tk.Label(frame, text="Radius of Curvature (mm)"))
        labels[1].pack(pady=0)
        scales = np.append(scales, tk.Scale(frame, from_=0, to=100, tickinterval=50,
        orient=tk.HORIZONTAL, variable = tk.DoubleVar(), length=250))
        scales[1].pack(padx=10, pady=2)
        scales[1].set(33)

        properties.update({"RadiusCurvature" : scales[1]})

        labels = np.append(labels, tk.Label(frame, text="Socket Radius (mm)"))
        labels[2].pack(pady=0)
        scales = np.append(scales, tk.Scale(frame, from_=0, to=100, tickinterval=50,
        orient=tk.HORIZONTAL, variable = tk.DoubleVar(), length=250))
        scales[2].pack(padx=10, pady=2)
        scales[2].set(35)

        properties.update({"Radius" : scales[2]})

        labels = np.append(labels, tk.Label(frame, text="Phase shift (degrees)"))
        labels[3].pack(pady=0)
        scales = np.append(scales, tk.Scale(frame, from_=0, to=179, tickinterval=90,
        orient=tk.HORIZONTAL, variable = tk.IntVar(), length=250))
        scales[3].pack(padx=10, pady=2)

        properties.update({"Phase" : scales[3]})

        labels = np.append(labels, tk.Label(frame, text="Transducer frequency (kHz):"))
        labels[4].pack(pady=0)
        entries = np.append(entries, tk.Entry(frame, width = 8))
        entries[0].pack(padx=10, pady=2)
        entries[0].insert(0,40.0)
        properties.update({"TransFreq" : entries[0]})

        return properties, labels, scales, entries, checkbox

    def reflector(self, frame):
        """ Method adds widgets relevant to a reflector """
        labels, scales, entries = [], [], []
        properties = {}

        var_checkbox = tk.BooleanVar()
        checkbox = tk.Checkbutton(frame, text="Concave", variable=var_checkbox, onvalue=True, offvalue=False)
        checkbox.pack(pady=10)

        properties.update({"Concave" : var_checkbox})

        labels = np.append(labels, tk.Label(frame, text="Vertical Displacement (mm)"))
        labels[0].pack(pady=0)
        scales = np.append(scales, tk.Scale(frame, from_=0, to=100, tickinterval=50,
        orient=tk.HORIZONTAL, variable = tk.DoubleVar(), length=250))
        scales[0].pack(padx=10, pady=2)
        scales[0].set(68)

        properties.update({"Displacement" : scales[0]})

        labels = np.append(labels, tk.Label(frame, text="Radius of Curvature (mm)"))
        labels[1].pack(pady=0)
        scales = np.append(scales, tk.Scale(frame, from_=0, to=100, tickinterval=50,
        orient=tk.HORIZONTAL, variable = tk.DoubleVar(), length=250))
        scales[1].pack(padx=10, pady=2)
        scales[1].set(33)

        properties.update({"RadiusCurvature" : scales[1]})

        labels = np.append(labels, tk.Label(frame, text="Socket Radius (mm)"))
        labels[2].pack(pady=0)
        scales = np.append(scales, tk.Scale(frame, from_=0, to=100, tickinterval=50,
        orient=tk.HORIZONTAL, variable = tk.DoubleVar(), length=250))
        scales[2].pack(padx=10, pady=2)
        scales[2].set(30)

        properties.update({"Radius" : scales[2]})

        return properties, labels, scales, entries, checkbox

class Pager:
    def __init__(self):
        self.labels_page1 = []
        self.labels_page2 = []
        self.properties_page1 = []
        self.properties_page2 = []
        self.scales_page1 = []
        self.scales_page2 = []
        self.entries_page1 = []
        self.entries_page2 = []

    def toggle_page1(self, selection, frame):
        type = types()

        if selection.get() == "Transducer":
            if len(self.labels_page1) > 0:
                for i in range(len(self.labels_page1)):
                    self.labels_page1[i].destroy()
                for i in range(len(self.scales_page1)):
                    self.scales_page1[i].destroy()
                for i in range(len(self.entries_page1)):
                    self.entries_page1[i].destroy()
                self.checkbox_page1.destroy()

            self.properties_page1, self.labels_page1, self.scales_page1, self.entries_page1, self.checkbox_page1 = type.transducer(frame)

        if selection.get() == "Array":
            if len(self.labels_page1) > 0:
                for i in range(len(self.labels_page1)):
                    self.labels_page1[i].destroy()
                for i in range(len(self.scales_page1)):
                    self.scales_page1[i].destroy()
                for i in range(len(self.entries_page1)):
                    self.entries_page1[i].destroy()
                self.checkbox_page1.destroy()

            self.properties_page1, self.labels_page1, self.scales_page1, self.entries_page1, self.checkbox_page1 = type.array(frame)

    def toggle_page2(self, selection, frame):
        type = types()
        if selection.get() == "Reflector":
            if len(self.labels_page2) > 0:
                for i in range(len(self.labels_page2)):
                    self.labels_page2[i].destroy()
                for i in range(len(self.scales_page2)):
                    self.scales_page2[i].destroy()
                for i in range(len(self.entries_page2)):
                    self.entries_page2[i].destroy()
                self.checkbox_page2.destroy()

            self.properties_page2, self.labels_page2, self.scales_page2, self.entries_page2, self.checkbox_page2 = type.reflector(frame)

        if selection.get() == "Array":
            if len(self.labels_page2) > 0:
                for i in range(len(self.labels_page2)):
                    self.labels_page2[i].destroy()
                for i in range(len(self.scales_page2)):
                    self.scales_page2[i].destroy()
                for i in range(len(self.entries_page2)):
                    self.entries_page2[i].destroy()
                self.checkbox_page2.destroy()

            self.properties_page2, self.labels_page2, self.scales_page2, self.entries_page2, self.checkbox_page2 = type.array(frame)


    def compute_potential(self, masterframe_inst):
        top_selection = masterframe_inst.top_selection
        bot_selection = masterframe_inst.bot_selection

        if top_selection.get() != "Select Type" or bot_selection.get() != "Select Type":
            bot_properties, top_properties, medium_properties = definedicts(self, masterframe_inst)

            start = time.time()
            relative_potential, pressure, x_span, z_span = MatrixMethod(medium_properties,top_properties,bot_properties)
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
            fig = plt.figure(figsize=(8, 6), dpi= 80, facecolor='w', edgecolor='k')
            cs = plt.contourf(_relative_potential, np.linspace(minPotential,maxPotential,1000), cmap='bone', extend='both')
            plt.colorbar()
            plt.show()
        else:
            print("Please define system geometry")

    def render(self, masterframe_inst):

        if masterframe_inst.top_selection.get() != "Select Type" or masterframe_inst.bot_selection.get() != "Select Type":

            bot_properties, top_properties, medium_properties = definedicts(self, masterframe_inst)

            Ux, Uy, Uz, Vx, Vy, Vz = CreateGeometry(top_properties, bot_properties)

            fig = plt.figure()
            ax = Axes3D(fig)
            ax.scatter3D(Vx, Vy, Vz, marker='.')
            ax.scatter3D(Ux, Uy, Uz, marker='.')
            ax.set_xlim([-0.05,0.05])
            ax.set_ylim([-0.05,0.05])
            ax.set_zlim([-0.065,0.065])
            plt.show()
        else:
            print("Please define system geometry")

    def compute_pressure(self, masterframe_inst):
        top_selection = masterframe_inst.top_selection
        bot_selection = masterframe_inst.bot_selection

        if top_selection.get() != "Select Type" or bot_selection.get() != "Select Type":
            bot_properties, top_properties, medium_properties = definedicts(self, masterframe_inst)

            start = time.time()
            relative_potential, pressure, x_span, z_span = MatrixMethod(medium_properties,top_properties,bot_properties)
            end = time.time()
            diff = end - start
            print("Total time elapsed was %.4f" % diff, "seconds")
            x = len(x_span)
            z = len(z_span)
            _pressure = np.real(pressure).reshape([z, x])
            maxPressure = np.max(_pressure)
            minPressure = np.min(_pressure)
            fig = plt.figure(figsize=(8, 6), dpi= 80, facecolor='w', edgecolor='k')
            cs = plt.contourf(_pressure, np.linspace(minPressure,maxPressure,1000), cmap='bone', extend='both')
            plt.colorbar()
            plt.show()
        else:
            print("Please define system geometry")

class masterframe(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        self.frames = {}
        self.drops = {}
        self.btns = {}
        self.drop_selections = {}

    def setouterproperties(self):
        self.geometry("750x850")
        self.iconbitmap(".\icon.ico")
        self.title("Simulation Platform for Acoustic Levitation Traps (SALT)")

    def placemainwidgets(self):
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        pager = Pager()

        self.top_selection = tk.StringVar()
        self.top_selection.set("Select Type")
        self.frames[0] = tk.LabelFrame(self, text="Top component", padx=10, pady=10)
        self.frames[0].pack(padx=10, pady=10)
        self.frames[0].place(x=20, y=0)
        self.drops[0] = tk.OptionMenu(self.frames[0], self.top_selection,
        "Array", "Transducer", command=lambda x:pager.toggle_page1(self.top_selection, self.frames[0]))
        self.drops[0].pack()

        self.bot_selection = tk.StringVar()
        self.bot_selection.set("Select Type")
        self.frames[1] = tk.LabelFrame(self, text="Bottom components", padx=10, pady=10)
        self.frames[1].pack(padx=10, pady=10)
        self.frames[1].place(x=420, y=0)
        self.drops[1] = tk.OptionMenu(self.frames[1], self.bot_selection,
        "Array", "Reflector", command=lambda x:pager.toggle_page2(self.bot_selection, self.frames[1]))
        self.drops[1].pack()

        self.frames[2] = tk.LabelFrame(self, text="Medium Properties")
        self.frames[2].pack(padx=5,pady=5,anchor=tk.CENTER)

        sub_frame = tk.LabelFrame(self.frames[2], text="Speed of sound (m/s)")

        self.medium_c = tk.Entry(sub_frame, width=8,justify=tk.CENTER)
        self.medium_c.insert(0,343)
        self.medium_c.pack(padx=5,pady=5,side=tk.BOTTOM)
        sub_frame.pack(side=tk.LEFT)

        sub_frame2 = tk.LabelFrame(self.frames[2], text="Zero-density (kg/m^3)")

        self.medium_density = tk.Entry(sub_frame2, width=8, justify=tk.CENTER)
        self.medium_density.insert(0,1.2)
        self.medium_density.pack(padx=5,pady=5,side=tk.BOTTOM)
        sub_frame2.pack(side=tk.RIGHT)

        self.frames[3] = tk.LabelFrame(self, text="Options", padx=10, pady=10)
        self.frames[3].pack(padx=10,pady=10,side=tk.BOTTOM)


        self.btns[0] = tk.Button(self.frames[3], text="Display Geometry",
        command=lambda:pager.render(self))
        self.btns[0].pack(pady=10, padx =5, side=tk.LEFT)

        self.btns[1] = tk.Button(self.frames[3], text="Compute Relative Acoustic Potential",
        command=lambda:pager.compute_potential(self))
        self.btns[1].pack(pady=10, padx = 5, side=tk.LEFT)

        self.btns[2] = tk.Button(self.frames[3], text="Compute Pressure",
        command=lambda:pager.compute_pressure(self))
        self.btns[2].pack(pady=10, padx = 5, side=tk.LEFT)

        # self.btns[3] = tk.Button(self.frames[3], text="Animate Phase Shift",
        # command=lambda:pager.animate(self))
        # self.btns[3].pack(pady=10, padx = 5, side=tk.LEFT)

root = masterframe()
root.setouterproperties()
root.placemainwidgets()
root.mainloop()
