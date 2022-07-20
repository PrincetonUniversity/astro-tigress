import sys
sys.path.append("./module/") #add path for import
import tigress_read #scripts for reading data

#master directory of where the data is stored
dir_master = "./data/" #this is default
model_id = "R8_4pc" #name of the simulation model

model = tigress_read.Model(model_id, dir_master=dir_master) #reading the model information

#one can download a data snapshot using `download` method
#dataset can be one of ["MHD", "chem", "CO_lines", "history", "input", "all"]
model.download(300,dataset="history")
model.download(300,dataset="input")
model.download(300,dataset="MHD")


