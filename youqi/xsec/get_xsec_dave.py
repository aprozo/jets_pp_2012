# drupal.star.bnl.gov/STAR/blog/djs232/pythia6-xsections-pp-embeeding-pau-200-gev-2015-collisions 

data = dict()
weighted = dict()

with open("Xsec_log_lines",'r') as f_in:
    for L in f_in:
        if ("Binary" in L):
            continue
        L = L.split('pt')[1]
        L = L.replace("infty","99")
        A = L.split('_')
        #print(A)
        if (L[0][0] == "_"): # some are like pt_2_3
            key = (int(A[1]),int(A[2]))
        else: # others are like pt3_4
            key = (int(A[0]),int(A[1]))
        A = L.split()
        nevents = int(A[7])
        Xsection = float(A[10].replace("D","E"))

        if key in data:
            data[key].append((nevents,Xsection))
            weighted[key][0] += nevents
            weighted[key][1] += nevents*Xsection
        else:
            data[key] = [(nevents,Xsection),]
            weighted[key] = [nevents, nevents*Xsection]

# print to intermediate file:
with open("Xsec-nEvents",'w') as f_out:
    keys = list(data.keys())
    keys.sort()
    for K in keys:
        f_out.write('pthatrange: %-4i - %-i\n'%(K[0],K[1]))
        for val in data[K]:
            f_out.write("%i %g\n"%(val[0],val[1]))
# print weighted values:
with open('Xsec-pythia.txt','w') as f_out:
    keys = list(weighted.keys())
    keys.sort()
    f_out.write("(pthatrange)  nevents  weighted-Xsection\n")
    for K in keys:
        f_out.write(" %2i,%2i,      %8i,  %16g\n" %(K[0],K[1],weighted[K][0], weighted[K][1]/weighted[K][0]))
