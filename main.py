#!/usr/bin/env python

from gui import masterframe

if __name__ == "__main__":
    root = masterframe()
    root.setouterproperties()
    root.placemainwidgets()
    root.mainloop()
