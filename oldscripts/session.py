import BootTest as bt
from array import array
ar = array('d',[1,2,3,4])
rd = [ar,ar+1,ar+2,ar+3]
rd = [ar,ar+[1,1,1,1],ar+[2,2,2,2],ar+[3,3,3,3]]
rd = [ar,ar+array('d',[1,1,1,1]),ar+array('d',[2,2,2,2]),ar+array('d',[3,3,3,3])]
rd
ar = [array('d',[1,2,3]),array('d',[4,5,6]),array('d',[7,8,9])]
ar
bt.CreateBoot(ar,10,.9)
data = bt.CreateBoot(ar,10,.9)
data[1].avg
data[1].Avg
data[0].Avg
data[2].Avg
import readline
readline.write_history_file('./session.py')
