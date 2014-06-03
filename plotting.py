
#R9   = 9 + cv.normal(N, ncols = 1, mean = 0.0, std = 0.01)
#RX9  = 9 + cv.normal(N, ncols = 1, mean = 0.0, std = 0.01)

#R14   = 14 + cv.normal(N, ncols = 1, mean = 0.0, std = 0.01)
#RX14  = 14 + cv.normal(N, ncols = 1, mean = 0.0, std = 0.01)

#R30  = 30 + cv.normal(N, ncols = 1, mean = 0.0, std = 0.01)
#RX30 = 30 + cv.normal(N, ncols = 1, mean = 0.0, std = 0.01)

#R39  = 39 + cv.normal(N, ncols = 1, mean = 0.0, std = 0.01)
#RX39 = 39 + cv.normal(N, ncols = 1, mean = 0.0, std = 0.01)

#### Absolute Time
#f = plt.figure()
#ax = f.add_subplot(1,1,1)
#ax.set_yscale('log')

#plt.plot(R9,time9[::,0],'ro')
#plt.plot(RX9,elapsedX9,'bo')
#plt.plot(R9,time9[::,3],'mx')
#plt.plot(RX9,elsX9,'cx')

#plt.plot(R14,time14[::,0],'ro')
#plt.plot(RX14,elapsedX14,'bo')
#plt.plot(R14,time14[::,3],'mx')
#plt.plot(RX14,elsX14,'cx')

#plt.plot(R30,time30[::,0],'ro')
#plt.plot(RX30,elapsedX30,'bo')
#plt.plot(R30,time30[::,3],'mx')
#plt.plot(RX30,elsX30,'cx')

#plt.plot(R39,time39[::,0],'ro')
#plt.plot(RX39,elapsedX39,'bo')
#plt.plot(R39,time39[::,3],'mx')
#plt.plot(RX39,elsX39,'cx')


#plt.show()

#vec = [0,1,2,3]

### Speed up
#f = plt.figure()
#ax = f.add_subplot(1,1,1)
#ax.set_yscale('log')
#plt.plot(vec,[Speedup9,Speedup14,Speedup30, Speedup39])
#plt.show()

### mean of time/ Bar plots
#mte9 = sum(time9[::,0])/time9.size[0]
#mteX9 = sum(elapsedX9)/elapsedX9.size[0]
#Speedup9 = mteX9/mte9

#mte14 = sum(time14[::,0])/time14.size[0]
#mteX14 = sum(elapsedX14)/elapsedX14.size[0]
#Speedup14 = mteX14/mte14

#mte30 = sum(time30[::,0])/time30.size[0]
#mteX30 = sum(elapsedX30)/elapsedX30.size[0]
#Speedup30 = mteX30/mte30

#mte39 = sum(time39[::,0])/time39.size[0]
#mteX39 = sum(elapsedX39)/elapsedX39.size[0]
#Speedup39 = mteX39/mte39

#m9 = cv.matrix(np.mean(time9,axis=0))
#m14 = cv.matrix(np.mean(time14,axis=0))
#m30 = cv.matrix(np.mean(time30,axis=0))
#m39 = cv.matrix(np.mean(time39,axis=0))

#M = cv.matrix([[m9/m9[0,0]],[m14/m14[0,0]],[m30/m30[0,0]], [m39/m39[0,0]]])
#tot = (M[1,::] + M[2,::] + M[3,::])
#tot_ = cv.matrix(1,(1,4)) - tot

#width = 0.35
#ind = np.arange(4)
#p1 = plt.bar(vec,M[1,::], width, color = 'r')
#p2 = plt.bar(vec,M[2,::], width, color = 'g',bottom=M[1,::])
#p3 = plt.bar(vec,M[3,::], width, color = 'y',bottom=M[2,::])
#p4 = plt.bar(vec,tot_, width, color = 'b',bottom=tot)


#plt.ylabel('time [s]')
#plt.title('network size')
#plt.xticks(ind+width/2., ('case9', 'case14', 'case30','case39') )
#plt.yticks(np.arange(0,1.1,0.1))
#plt.legend( (p1[0], p2[0], p3[0], p4[0]), ('Init', 'decomposition', 'solver', 'other'),loc=10)

#plt.show()

#### growth original vs. decomposed (realtive time)
#Growth = [mte9, mte14, mte30, mte39]
#Growth = [i/mte9 for i in Growth]

#GrowthX = [mteX9, mteX14, mteX30, mteX39]
#GrowthX = [i/mteX9 for i in GrowthX]

#f = plt.figure()
#ax = f.add_subplot(1,1,1)
#ax.set_yscale('log')

#plt.plot(vec, Growth, 'ro')
#plt.plot(vec, GrowthX, 'bo')
#plt.xlim((-0.5,3.5))
#plt.show()

#f = open('DataGen/timeMeasurements9', 'w')
#np.save(f,np.array(time9))

#f = open('DataGen/timeMeasurements14', 'w')
#np.save(f,np.array(time14))

#f = open('DataGen/timeMeasurements30', 'w')
#np.save(f,np.array(time30))

#f = open('DataGen/timeMeasurements39', 'w')
#np.save(f,np.array(time39))
