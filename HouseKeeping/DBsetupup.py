import csv

def genHECSV():
	with open("HighEnergyGroups.csv", "w") as DB:
		DB.write("N1C=NN=N1 |c:1,3|\n")
		DB.write("N1C=NN=N1 |c:1,3|\n")
		DB.write("N1C=NN=N1 |c:1,3|\n")
		DB.write("N1C=NN=N1 |c:1,3|\n")

genHECSV()