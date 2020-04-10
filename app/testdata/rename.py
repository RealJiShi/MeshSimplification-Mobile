import os
import os.path

files = os.listdir("./")
ext = ".off"
index = 46

for filename in files:
	portion = os.path.splitext(filename)

	if portion[1] == ext:
	    newname = "mesh" + str(index) + ".off"
	    os.chdir("./")
	    os.rename(filename, newname)
	    index += 1
