file = open("momentum5000.000000_xPosition0.098000_yPosition0.098000_theta0.300000_phi0.600000_NO_wParaMeter_test_Pion.txt","r")

lineList = file.readlines()

ll = []
for temp in lineList:
    tempSplit = temp.split()
    ll.append(tempSplit[len(tempSplit)-1])


string = lineList[len(lineList)-1]

print string
stringList = string.split()

i = 0

v5andErr = []
vWar = []
test5Parameters = []
test5Warnings = []
start = 0
for temp in stringList:
    if (i>= 0+start and i <10+start and i <= 10000):
        test5Parameters.append(temp)
    elif (i >=10+start and i <= 10000):
        test5Parameters.append(temp)
        v5andErr.append(" ".join(test5Parameters))
        print i," ".join(test5Parameters)
        test5Parameters = []
        #i = i - 1
        start=10+start
    elif (i >10000 and i >=0+start and i <4 + start):
        test5Warnings.append(temp)
    elif (i > 10000 and i >= 4+start):
        test5Warnings.append(temp)
        vWar.append(" ".join(test5Warnings))
        print i, " ".join(test5Warnings)
        test5Warnings = []
        start = 4 + start
    i=i+1

#print stringList[len(stringList)-1]

text_file = open("Output.txt", "w")
for i in range(0,len(vWar)):
    text_file.write("%d %s %s %s \n" %(i, v5andErr[i], vWar[i], ll[i]))

text_file.close();

