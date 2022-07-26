import sys
sys.path.insert(0, "./module/") #add path for import
import tigress_read #scripts for reading data

#master directory of where the data is stored
dir_master = "./data/" #this is default

#one can download a data snapshot using `download` method
#dataset can be one of ["MHD", "chem", "CO_lines", "history", "input", "all"]
model_id_4pc = "R8_4pc" #name of the simulation model
model_4pc = tigress_read.Model(model_id_4pc, dir_master=dir_master) #reading the model information
model_4pc.download(300,dataset="history")
model_4pc.download(300,dataset="input")
model_4pc.download(300,dataset="MHD")

model_id_2pc = "R8_2pc" #name of the simulation model
model_2pc = tigress_read.Model(model_id_2pc, dir_master=dir_master) #reading the model information
model_2pc.download(300,dataset="history")
model_2pc.download(300,dataset="input")
model_2pc.download(300,dataset="MHD")
#model_2pc.download(300,dataset="chem") #large file 15GB
model_2pc.download(300, dataset="CO_lines", iline=1) #CO(J=1-0)
model_2pc.download(300, dataset="CO_lines", iline=2) #CO(J=2-1)
