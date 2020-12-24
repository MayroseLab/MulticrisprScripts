from tkinter import *
import call_MULTICRISPR_on_families
import re
#from tkinter import tkMessageBox

class App(object):
    def __init__(self):
        self.root = Tk()
        self.root.geometry("700x140")
        self.root.wm_title("MultiCRISPR 1.0")
        self.label = Label (self.root, text= "Enter folder path")
        self.label.pack()


        self.entrytext = StringVar()
        Entry(self.root, textvariable=self.entrytext, width = 70).pack()

        self.buttontext = StringVar()
        self.buttontext.set("Find")
        Button(self.root, textvariable=self.buttontext, command=self.clicked1).pack()

        self.label = Label(self.root, text="")
        self.label.pack()

        self.root.mainloop()


    def clicked1(self):
        the_input = self.entrytext.get()
        ##call "call multyCRISPR"
        #result = int(input)*2
        if "," in the_input:
#            path, Omega = the_input.split(", ") ##old version
            path, parameters = the_input.split(", ")
            if "G" in parameters:
                start_with_G = True
            else:
                start_with_G = False
            if " 1 " in parameters:
                Omega = 1
            else:
                OmegaCompiled = re.compile("0\.\d+")
                OmegaInLst = re.findall(OmegaCompiled, parameters)
                if len(OmegaInLst) != 0:
                    Omega = float(OmegaInLst[0])
                else:
                    Omega = 0.011
            lengthCompiled = re.compile("\d\d-\d\d")
            lengthInLst = re.findall(lengthCompiled, parameters)
            if len(lengthInLst) != 0:
                min_and_max_length = lengthInLst[0].split("-")
                min_length, max_length = int(min_and_max_length[0]), int(min_and_max_length[1])
            else:
                min_length, max_length = 20, 20
        else:
            path = the_input
            Omega = 0.011
            start_with_G = False
            min_length, max_length = 20, 20
        result = call_MULTICRISPR_on_families.call_using_CasSites(path, Omega, min_length, max_length,start_with_G)
        self.label.configure(text=result)

    def button_click(self, e):
        pass

App()
