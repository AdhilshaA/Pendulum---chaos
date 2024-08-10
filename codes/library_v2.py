import math as m
import numpy as np


def parse_trajectory_data(file_path, offset = 3, updates = False):
    # Parse the trajectory data from the txt file into numpy matrix
    # offset is the number of lines to skip from the top of the file 
    #       (usually first 3 lines are the header, so offset = 3)
    #       (in first line of data, we add the lines to skip manually, this is detected using 'auto')
    file_name = file_path.split('\\')[-1]
    if updates: print(f"Opening {file_path} and reading data !")
    with open(file_path, 'r') as f:
        data = []
        lines = f.readlines()
        if offset == 'auto':
            try: offset = int(lines[0].split(':')[-1])
            except: offset = 3
        skipped = []
        for i in range(offset, len(lines)):
            try: 
                row = np.float64(lines[i].split('\t'))
            except:
                print(f"Data not in proper format for {file_name} at line {i+1} !")
                row = lines[i].split('\t')
                if row[1] == '': # mass1 is missing
                    j = 1
                    if row[3] != '': k = 3 #mass2 known
                    else : k = 5 #mass2 also missing
                else: # mass1 is known and mass2 is missing
                    j = 3
                    k = 5
                temp = row[k:5]
                row = row[:j]
                while j < k:
                    row.append(0)
                    skipped.append((i,j))
                    j += 1
                if temp:
                    row = row + temp
                row = [np.float64(val) for val in row]
            data.append(row)
    if updates: print(f"Data read successfully !")
    data = np.array(data)
    
    if skipped: #handling skipped data
        print(f"interpolating {len(skipped)} coordinates!")
        
        # removing offset from skipped as the data starts from 0th row
        skipped = [(i-offset,j) for i,j in skipped]
        
        for i, j in skipped:
            if (i-1,j) in skipped: # must have been processed in the longers skips in previous iterations
                continue
            if (i+1,j) in skipped: # longer skips handling
                start = (i-1,j)
                end = (i+2,j)
                while end in skipped:
                    end = (end[0]+1,end[1])
                diff_each = (data[end] - data[start]) / ((end[0] - start[0]))
                for k in range(i,end[0]):
                    if data[k,j] == 0: data[k,j] = data[k-1,j] + diff_each
            else:
                data[i,j] = (data[i-1,j] + data[i+1,j]) / 2
    
    # printing corrected ones
    for i,j in skipped:
        print(f'line {i+1}, column {j+1}:  {data[i,j]} [{data[i-1,j],data[i+1,j]}]' )
     
    return np.array(data)

def RK4(dydxlist,x0,y0list,x1,dx,print_updates = False):
    # Using simple RK4 to solve coupled ODE equations.
    # solves dydxlist list of ODE fns. with intial y0 values at x0.
    # return graph data

    # x0 = np.float64(x0)
    # x1 = np.float64(x1)
    # dx = np.float64(dx)
    # y0list = [np.float64(i) for i in y0list]
    
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
    
    if print_updates:
        checkpoints = [i for i in range(int(x0),int(x1)+2,1)]
        curr_idx = 1
    
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
         
        if print_updates and int(x0) == checkpoints[curr_idx]:
            print(f"completing {int(x0)}/{int(x1)}...")
            curr_idx += 1
        
        for i in range(n):
            datY[i].append(y0list[i])
        datX.append(x0)
        
    del k1,k2,k3,k4,tempy0list
    
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
    
    del k1,k2,k3,k4,tempy0list
    
    return x0, y0list[:]

def ASRK4(dydxlist,x0,y0list,xf,h,tolerance,errprioirity = None,print_updates = False):

    # x0 = np.float64(x0)
    # xf = np.float64(xf)
    # h = np.float64(h)
    # y0list = [np.float64(i) for i in y0list]

    if errprioirity is None:
        errprioirity = [i for i in range(len(y0list))]

    n = len(y0list) #no. of variables

    # print(x0,xf)

    #error for each variable
    # yerr = [0]*n # unwanted ig
    counts = [0]

    datX = [x0]
    datY = []
    for j in range(n):
        datY.append([y0list[j]])
    
    if print_updates:    
        checkpoints = [i for i in range(int(x0),int(xf)+2,1)]
        curr_idx = 1
    
    k = 0
    while x0 < xf:
        k += 1
        # print('{}th step, '.format(k),end = '')
        rho = np.float64(0)
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
        
        if print_updates  and x0 >= checkpoints[curr_idx]:
            print(f"completing {int(x0)}/{int(xf)}...")
            curr_idx += 1
        
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

import os

def save_matrix_as_txt(path, filename, matrix):
    
    # Save the matrix as txt file
    with open(f"{path}\{filename}", 'w') as f:
        for row in matrix:
            row_str = [str(val) for val in row]
            f.write('\t'.join(row_str) + '\n')


def read_matrix_from_txt(path, filename):
    # Create the file path
    
    # Read the matrix from txt file
    with open(f"{path}\{filename}", 'r') as f:
        matrix = []
        for line in f:
            row = [np.float64(val) for val in line.split('\t')]
            matrix.append(row)
    
    return matrix