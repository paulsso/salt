import tkinter as tk

class masterframe(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        mainwidgets = placemainwidgets(masterframe).placemainwidgets()

class placemainwidgets(masterframe):
    def placemainwidgets():
        container = tk.Frame(masterframe)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        masterframe.geometry("790x800")
        masterframe.iconbitmap(".\icon.ico")
        masterframe.title("Simulation Platform for Acoustic Levitation Traps (SALT)")
        masterframe.frames = {}

        masterframe.clicked1 = StringVar()
        masterframe.clicked1.set("Window 1")
        frame = LabelFrame(masterframe, text="Top component", padx=10, pady=10)
        frame.pack(padx=10, pady=10)
        frame.place(x=50, y=0)
        drop = OptionMenu(frame, masterframe.clicked1, "Array", "Transducer", command=show1)
        drop.pack()

        masterframe.clicked2 = StringVar()
        masterframe.clicked2.set("Window 2")
        frame2 = LabelFrame(masterframe, text="Bottom components", padx=10, pady=10)
        frame2.pack(padx=10, pady=10)
        frame2.place(x=420, y=0)
        drop2 = OptionMenu(frame2, masterframe.clicked2, "Array", "Reflector", command=show2)
        drop2.pack()

        btn1 = Button(masterframe, text="Display Geometry", command=lambda:render(masterframe))
        btn1.pack(side=BOTTOM,pady=10)
        btn2 = Button(masterframe, text="Compute Relative Acoustic Potential", command=lambda:compute(masterframe))
        btn2.pack(side=BOTTOM,pady=10)

root = masterframe()
root.mainloop()
