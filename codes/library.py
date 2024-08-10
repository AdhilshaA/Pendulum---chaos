import math as m

def RK4(dydxlist,x0,y0list,x1,dx):
    # Using simple RK4 to solve coupled ODE equations.
    # solves dydxlist list of ODE fns. with intial y0 values at x0.
    # return graph data

    if len(dydxlist) != len(y0list):
        print('Boundary values doesnt match with the coupled ODEs')
        return None

    x1 -= dx/2 #for proper stopping of iterations

    # setting the common k1,k2,k3,k4 values for each variable
    n = len(y0list) #no. of variables
    k1 = [0]*n
    k2 = [0]*n
    k3 = [0]*n
    k4 = [0]*n
    tempy0list = [0]*n #intermediate y values needed

    # setting the solution lists
    datY = []
    for i in range(n):
        datY.append([y0list[i]])
    datX = [x0]

    # iterations till you reach x1
    while x0 < x1:

        #calculating k1,k2,k3,k4
        for i in range(n):
            k1[i] = dx*dydxlist[i](y0list,x0)
        
        for i in range(n):
            tempy0list[i] = y0list[i] + (k1[i] / 2)
        for i in range(n):
            k2[i] = dx*dydxlist[i](tempy0list, (x0 + (dx/2)))

        for i in range(n):
            tempy0list[i] = y0list[i] + (k2[i] / 2)
        for i in range(n):
            k3[i] = dx*dydxlist[i](tempy0list, (x0 + (dx/2)))
            
        for i in range(n):
            tempy0list[i] = y0list[i] + k3[i]
        for i in range(n):
            k4[i] = dx*dydxlist[i](tempy0list, (x0 + dx))

        #getting the next set of y values
        for i in range(n):
            y0list[i] += ((k1[i] + (2 * k2[i]) + (2 * k3[i]) + (k4[i])) / 6)
        x0 += dx
        
        for i in range(n):
            datY[i].append(y0list[i])
        datX.append(x0)
    return datX, datY

def RK4_singlestep(dydxlist,x0,y0list,htry):

   

    if len(dydxlist) != len(y0list):
        print('current values doesnt match with the coupled ODEs')
        return None

    # setting the common k1,k2,k3,k4 values for each variable
    n = len(y0list) #no. of variables
    k1 = [0]*n
    k2 = [0]*n
    k3 = [0]*n
    k4 = [0]*n
    tempy0list = [0]*n #intermediate y values needed  

    #calculating k1,k2,k3,k4 for single step of size h
    for i in range(n):
        k1[i] = htry*dydxlist[i](y0list,x0)
    
    for i in range(n):
        tempy0list[i] = y0list[i] + (k1[i] / 2)
    for i in range(n):
        k2[i] = htry*dydxlist[i](tempy0list, (x0 + (htry/2)))

    for i in range(n):
        tempy0list[i] = y0list[i] + (k2[i] / 2)
    for i in range(n):
        k3[i] = htry*dydxlist[i](tempy0list, (x0 + (htry/2)))
        
    for i in range(n):
        tempy0list[i] = y0list[i] + k3[i]
    for i in range(n):
        k4[i] = htry*dydxlist[i](tempy0list, (x0 + htry))

    #getting the next set of y values
    for i in range(n):
        y0list[i] += ((k1[i] + (2 * k2[i]) + (2 * k3[i]) + (k4[i])) / 6)
    x0 += htry
    
    return x0, y0list[:]

def ASRK4(dydxlist,x0,y0list,xf,h,tolerance,errprioirity = None):

    if errprioirity is None:
        errprioirity = [i for i in range(len(y0list))]

    n = len(y0list) #no. of variables

    # print(x0,xf)

    #error for each variable
    yerr = [0]*n
    counts = [0]

    datX = [x0]
    datY = []
    for j in range(n):
        datY.append([y0list[j]])

    k = 0
    while x0 < xf:
        k += 1
        # print('{}th step, '.format(k),end = '')
        rho = 0
        count = -1
        while rho < 1:
            # print("\nInside loop start")
            count += 1
            dx = 2*h
            # one double step
            y1list = RK4_singlestep(dydxlist,x0,y0list[:],2*h)[1]
            # two single steps
            x2, y2list = RK4_singlestep(dydxlist,x0,y0list[:],h)
            # print(y2list,"(y2list)")
            # print("x values",x2,x0+h)
            x2, y2list = RK4_singlestep(dydxlist,x2,y2list,h)
            # print('solution by single step = {}'.format(y1list))
            # print('solution by double step = {}'.format(y2list))
            
            # yerr = abs(y1list[0] - y2list[0]) / 30

            yerr = 0
            for i in errprioirity:
                yerr += ((y1list[i] - y2list[i])/30)**2
            yerr = m.sqrt(yerr)
            # print(yerr)
            #maximum error 
            # errmax = max(yerr)

            try: 
                rho =  (tolerance * h) / yerr
                # print('rho is',rho)
                # print("temp h is",h*(rho**0.25))
                h = min( h*(rho**(1/4)), 2*h )
            except: 
                h = 2*h
                break
            # print("reduced to {}".format(h_new))
        # if count != 0: print("{} reductions".format(count))
        counts.append(count)
        x0 += dx
        y0list = y2list
        datX.append(x0)
        for j in range(n):
            datY[j].append(y0list[j])
        # print('current x0 is',x0)
    # print(datX,datY)

    return datX,datY,counts

def print_coltable(data):
    # prints tables with heading as keys and columns as values from the dictionary data given
    # first column is formatted as for serial number  (integers)
    # all data is assumed to be numbers (floats to be specific, except from the first col, which is int)

    cols = len(data)
    key = list(data.keys())

    rows = 0
    for k in key:
        if len(data[k]) > rows:
            rows = len(data[k])
    
    maxchar = []
    for k in key:
        max = len(k)
        for value in data[k]:
            if int(value)//1 == 0 and value < 0 and max < 8:
                max = 8
                continue
            if (len(str(int(value)))+7) > max:
                max = len(str(int(value)))
        maxchar.append(max + 2)
    # print(maxchar)

    line = '|'.join(['-' * spaces for spaces in maxchar])
    line = ''.join(['|',line,'|'])
    i = 0
    print(line)
    print('|',end='',sep ='')
    for k in key:
        max = maxchar[i]
        spaces = (max - len(k))
        start =  spaces // 2
        stop = spaces - start
        # if k == 'N': print('spaces',max,len(k),max - len(k),spaces)
        print(' '*start,k,' '*stop,'|',end = '',sep='')
        i += 1
    print('')
    print(line)

    for i in range(rows):
        print('|',end='',sep='')
        spaces = maxchar[0] - len(str(data[key[0]][i]))
        start = (spaces) // 2
        stop = spaces - (start)
        print(' '*start,data[key[0]][i],' '*stop,'|',end = '',sep='')
                
        for j in range(1,cols):
            try: spaces = maxchar[j] - (len(str(int(data[key[j]][i]))) + 7)
            except:
                start = maxchar[j] // 2
                stop = maxchar[j] - (start + 1)
                print(' '*start,'-',' '*stop,'|',end = '',sep='')
                continue
            if int(data[key[j]][i])//1 == 0 and value < 0:
                spaces -= 1
            start = spaces // 2
            stop = spaces - start
            print(' '*start,end='',sep = '')
            print('{:f}'.format(data[key[j]][i]),end='',sep = '')
            print(' '*stop,'|',end='',sep = '')
        print('')   
    print(line)




    
    
    



