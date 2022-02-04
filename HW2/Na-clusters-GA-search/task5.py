from ase.db import connect
db6 = connect('./6/gadb.db')
db7 = connect('./7/gadb.db')
db8 = connect('./8/gadb.db')
atoms6 = db6.get('id=1').toatoms()
atoms7 = db7.get('id=1').toatoms()
atoms8 = db8.get('id=1').toatoms()
#print(atoms6)
