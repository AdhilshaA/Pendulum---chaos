# A library updated to all functions till now

from math import sqrt
import math as m
import matplotlib.pyplot as plt
    
def sum_odd(n):
    #returns sum of n odd numbers
    sum = 0
    for i in range(1,2 * n,2):
        sum += i
    return sum

def fact(n):
    #returns n factorial
    factorial = 1
    for i in range(n,1,-1):
        factorial = factorial * i
    return factorial
    
def sum_ap(start,cd,n):
    sum = 0
    term = start
    #returns the sum of first n terms in an AP specified by first term and common difference
    for i in range(n):
        sum += term
        term += cd
    return sum

def sum_hp(start,cd,n):
    sum = 0
    term = start
    #returns the sum of first n terms in an HP specified by first term and common difference of its corresponding AP
    for i in range(n):
        sum += 1 / term
        term += cd
    return sum

def sum_gp(start,cr,n):
    sum = 0
    term = start
    #returns the sum of first n terms in an GP specified by first term and common ratio
    for i in range(n):
        sum += term
        term = term * cr
    return sum

def mat_mult(A,B): 
    #finding AB matrix

    #verifiying matrix multiplication is possible (cols of A ?= rows of B)
    if len(A[0]) != len(B):
        return None

    #creating the sum matrix with zeroes
    sum = []
    for i in range(len(A)): #no. of row in sum is row of A
        sum.append(list(0 for i in range( len(B[0]) ))) #no. of col in sum is col of B

    #filling the sum matrix
    for row in range(len(sum)):
        for col in range(len(sum[0])):
            #finding each term in sum matrix
            temp_sum = 0
            for i in range(len(B)):
                temp_sum += A[row][i] * B[i][col]
                sum[row][col]=temp_sum

    return sum

def mat_scalar(A,k):
    for i in range(len(A)):
        for j in range(len(A[0])):
            A[i][j] = A[i][j]*k
    return A

def mat_add(A,B):
    if len(A) != len(B) or len(A[0]) == len(A[0]):
        return None
    for i in range(len(A)):
        for j in range(len(A[0])):
            A[i][j]= A[i][j] + B[i][j]
    return A

def print_mat2(A):
    #printing a given matrix A
    maxs = []
    for j in range(len(A[0])):  #finding maximum integer length in each column and saving to maxs list
        max = 0
        for i in range(len(A)):
            length = len(str(int(A[i][j])))
            if max < length:
                max = length
        maxs.append(max + 2)   #adding 2 for equal space left and right to the maximum long term
    
    #printing
    for i in range(len(A)):
        print('[ ',end='')
        for j in range(len(A[0])):
            val = A[i][j]
            length = len(str(int(val))) 
            if int(val) == 0 and val < 0: length += 1
            spaces = maxs[j] - length
            stop = (spaces) // 2       #finding how many spaces before and after each term to align centrally in each column
            start = spaces - stop
            print(' ' * start,'{:.2f}'.format(A[i][j]),' ' * stop, end = '')
        print(' ]')
    print('')

def print_mat4(A):
    #prints matrix with floats upto 8 decimal points

    if A is None:
        print('None')
        return
    #printing a given matrix A
    maxs = []
    for j in range(len(A[0])):  #finding maximum integer length in each column and saving to maxs list
        max = 0
        for i in range(len(A)):
            length = len(str(int(A[i][j])))
            if max < length:
                max = length
        maxs.append(max + 2)   #adding 2 for equal space left and right to the maximum long term
    
    #printing
    for i in range(len(A)):
        print('[ ',end='')
        for j in range(len(A[0])):
            val = A[i][j]
            length = len(str(int(val))) 
            if int(val) == 0 and val < 0: length += 1
            spaces = maxs[j] - length
            stop = (spaces) // 2       #finding how many spaces before and after each term to align centrally in each column
            start = spaces - stop
            print(' ' * start,'{:.4f}'.format(A[i][j]),' ' * stop, end = '')
        print(' ]')
    print('')  



def print_mat8(A):
    #prints matrix with floats upto 8 decimal points

    if A is None:
        print('None')
        return
    #printing a given matrix A
    maxs = []
    for j in range(len(A[0])):  #finding maximum integer length in each column and saving to maxs list
        max = 0
        for i in range(len(A)):
            length = len(str(int(A[i][j])))
            if max < length:
                max = length
        maxs.append(max + 2)   #adding 2 for equal space left and right to the maximum long term
    
    #printing
    for i in range(len(A)):
        print('[ ',end='')
        for j in range(len(A[0])):
            val = A[i][j]
            length = len(str(int(val))) 
            if int(val) == 0 and val < 0: length += 1
            spaces = maxs[j] - length
            stop = (spaces) // 2       #finding how many spaces before and after each term to align centrally in each column
            start = spaces - stop
            print(' ' * start,'{:.8f}'.format(A[i][j]),' ' * stop, end = '')
        print(' ]')
    print('')  

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


def mat_dot(A,B):
    #returns the dot product of two column matrices
    if len(A) != len(B) or len(A[0]) != 1 or len(B[0]) != 1: #works only if column matrice sof same length
        return None
    dotprdct = 0
    for row in range(len(A)):
        dotprdct += (A[row][0] * B[row][0])
    return dotprdct

def mat_copy(A):
    #returns a copy of the matrix without any links in memory
    A1 = []
    for i in range(len(A)):
        A1.append(A[i][:]) #splicing done to avoid changes linking
    return A1

def check_symmetry(A):
    #Function to check symmetry of A, return True or False
    if len(A) != len(A[0]): #not square matix
        return False

    for i in range(len(A)):
        for j in range(i + 1, len(A[0])): #checking only off-diagonal elements

            if A[i][j] != A[j][i]: #symmetry condition False
                return False
    return True 


class mat:
    
    def __check(x):
        if type(x) != mat:
            return False

    def copy(A):
        A1 = mat([[]])
        A1.data = []
        for i in range(len(A.data)):
            A1.data.append(A.data[i][:]) #splicing done to avoid changes linking
        return A1

    def __init__(self,A):
        if type(A) != list:
            print('Not a matrix !')
            return
        n = len(A[0])
        flag = 1
        for i in range(1,len(A)):
            if type(A[i]) != list or n != len(A[i]):
                flag = 0
                print('Not a matrix !')
                return
        if flag == 1:
            self.data = A

    def __add__(self,B):
        if mat.__check(B) == False:
            print ('matrix addition impossible !')
            return None

        sum = mat([[]])
        sum.data = []
        for i in range(len(self.data)):
            sum.data.append([])
            for j in range(len(self.data[0])):
                sum.data[i].append(self.data[i][j] + B.data[i][j])
        return sum
        
    def __sub__(self,B):
        if mat.__check(B) == False:
            print ('matrix addition impossible !')
            return None

        sum = mat([[]])
        sum.data = []
        for i in range(len(self.data)):
            sum.data.append([])
            for j in range(len(self.data[0])):
                sum.data[i].append(self.data[i][j] - B.data[i][j])
        return sum

    def __str__(self):
        return str(self.data)

    def table(self,decimals = 2):
        if decimals == 8:
            print_mat8(self.data)
        elif decimals == 4:
            print_mat4(self.data)
        else:
            print_mat2(self.data)
        

    def __mul__(self,B):
        if mat.__check(B) == False:
            if type(B) != int and type(B) != float:
                print ('cannot perform multiplication with matrix !')
                return None

        
        if type(B) == int or type(B) == float:
            sol = mat.copy(self)
            for i in range(len(self.data)):
                for j in range(len(self.data[0])):
                    sol.data[i][j] *= B
            return sol
        else:
            #finding AB matrix

            #verifiying matrix multiplication is possible (cols of A ?= rows of B)
            if len(self.data[0]) != len(B.data):
                print("data mismatch")
                return None

            #creating the sum matrix with zeroes
            sum = mat([[]])
            sum.data = []
            for i in range(len(self.data)): #no. of row in sum is row of A
                sum.data.append(list(0 for i in range( len(B.data[0]) ))) #no. of col in sum is col of B

            #filling the sum matrix
            for row in range(len(sum.data)):
                for col in range(len(sum.data[0])):
                    #finding each term in sum matrix
                    temp_sum = 0
                    for i in range(len(B.data)):
                        temp_sum += self.data[row][i] * B.data[i][col]
                        sum.data[row][col]=temp_sum

            return sum 

    def dot(self,B):
        #returns the dot product of two column matrices
        if len(self.data) != len(B.data) or len(self.data[0]) != 1 or len(B.data[0]) != 1: #works only if column matrice sof same length
            return None
        dotprdct = 0
        for row in range(len(self.data)):
            dotprdct += (self.data[row][0] * B.data[row][0])
        return dotprdct

    def __iadd__(self,B):
        if mat.__check(B) == False:
            print ('matrix addition impossible !')
            return None

        for i in range(len(self.data)):
            for j in range(len(self.data[0])):
                self.data[i][j] += B.data[i][j]
        return self
        
    def __isub__(self,B):
        if mat.__check(B) == False:
            print ('matrix addition impossible !')
            return None

        for i in range(len(self.data)):
            for j in range(len(self.data[0])):
                self.data[i][j] -= B.data[i][j]
        return self

    
    def __imul__(self,B):
        if mat.__check(B) == False:
            if type(B) != int and type(B) != float:
                print ('cannot perform multiplication with matrix !')

        
        if type(B) == int or type(B) == float:
            for i in range(len(self.data)):
                for j in range(len(self.data[0])):
                    self.data[i][j] *= B
            return self
        else:
            #finding AB matrix

            #verifiying matrix multiplication is possible (cols of A ?= rows of B)
            if len(self.data[0]) != len(B.data):
                return None

            #creating the sum matrix with zeroes
            sum = mat([[]])
            sum.data = []
            for i in range(len(self.data)): #no. of row in sum is row of A
                sum.data.append(list(0 for i in range( len(B.data[0]) ))) #no. of col in sum is col of B

            #filling the sum matrix
            for row in range(len(sum.data)):
                for col in range(len(sum.data[0])):
                    #finding each term in sum matrix
                    temp_sum = 0
                    for i in range(len(B.data)):
                        temp_sum += self.data[row][i] * B.data[i][col]
                        sum.data[row][col]=temp_sum
            self = mat.copy(sum)
            return self
    
    def __abs__(self):
        #only for vectors!

        sum = 0
        for i in range(len(self.data)):
            sum += self.data[i][0]**2
        return m.sqrt(sum)

    def check_symmetry(self):
        #Function to check symmetry of A, return True or False
        if len(self.data) != len(self.data[0]): #not square matix
            return False

        for i in range(len(self.data)):
            for j in range(i + 1, len(self.data[0])): #checking only off-diagonal elements

                if self.data[i][j] != self.data[j][i]: #symmetry condition False
                    return False
        return True 


class myComplex:

  def __init__(self, real, imag):
    self.real = real
    self.imag = imag

  def __str__(self): #for printing complex no. in (a + bi) format when called print on any class instance
    return '{} + ({})i'.format(self.real, self.imag)

  def __add__(self,z2): # Method to add it with z2 of same class by using + between the instances
    add_real = self.real + z2.real
    add_imag = self.imag + z2.imag
    z3 = myComplex(add_real, add_imag)
    return z3
  
  def __mul__(self,z2): # Method to multiply it with z2 of same class using * between instances
    mult_real = (self.real * z2.real) - (self.imag * z2.imag)
    mult_imag = (self.real * z2.imag) + (self.imag * z2.real)
    z3 = myComplex(mult_real, mult_imag)
    return z3
  
  def __abs__(self):  #method to find modulus by using abs() on the instance
    val = sqrt((self.real ** 2) + (self.imag ** 2))
    return val

def parse(file_name):
    # A function to parse data in a specific format. Read README file in the repository for more details
    with open(file_name) as file:
        lines = file.readlines()
        file.close()
    inputs = {}
    numlines = len(lines)
    line_index = 3
    while bool(lines[line_index].split()) is False:
        line_index += 1
    while line_index < numlines:
        line_index += 2
        type_name = lines[line_index - 1].split()
        if type_name[0] == 'int':
            inputs[type_name[1]] = int(lines[line_index].split()[0])
            line_index += 1
            continue
        
        elif type_name[0] == 'float':
            inputs[type_name[1]] = float(lines[line_index].split()[0])
            line_index += 1
            continue

        elif type_name[0] == 'str':
            inputs[type_name[1]] = lines[line_index][:-1]
            line_index += 1
            continue
        
        elif type_name[0] == 'complex':
            real_imag = list(map(float,lines[line_index].split()))
            inputs[type_name[1]] = myComplex(real_imag[0], real_imag[1])
            line_index += 1
            continue

        elif type_name[0] == 'int_list':
            inputs[type_name[1]] = list(map(int,lines[line_index].split()))
            line_index += 1
            continue

        elif type_name[0] == 'float_list':
            inputs[type_name[1]] = list(map(float,lines[line_index].split()))
            line_index += 1
            continue

        elif type_name[0] == 'str_list':
            inputs[type_name[1]] = []
            while lines[line_index][0] != '#':
                inputs[type_name[1]].append(lines[line_index][:-1])
                line_index += 1

        elif type_name[0] == 'int_mat':
            inputs[type_name[1]] = []
            while lines[line_index][0] != '#':
                inputs[type_name[1]].append(list(map(int,lines[line_index].split())))
                line_index += 1
        
        elif type_name[0] == 'float_mat':
            inputs[type_name[1]] = []
            while lines[line_index][0] != '#':
                inputs[type_name[1]].append(list(map(float,lines[line_index].split())))
                line_index += 1
            inputs[type_name[1]] = mat(inputs[type_name[1]])

        elif type_name[0] == 'str_mat':
            inputs[type_name[1]] = []
            while lines[line_index][0] != '#':
                inputs[type_name[1]].append(lines[line_index].split())
                line_index += 1
            inputs[type_name[1]] = mat(inputs[type_name[1]])

        elif type_name[0] == 'xydata':
            inputs[type_name[1]] = []
            inputs[type_name[2]] = []
            while lines[line_index][0] != '#':
                lines[line_index] = lines[line_index].split()
                try:
                    inputs[type_name[1]].append(float(lines[line_index][0]))
                    inputs[type_name[2]].append(float(lines[line_index][1]))
                except:
                    print('{} input inproper! check {}'.format(type_name[0],file_name))
                line_index += 1

        elif type_name[0] == 'xysigdata':
            inputs[type_name[1]] = []
            inputs[type_name[2]] = []
            inputs[type_name[3]] = []
            while lines[line_index][0] != '#':
                lines[line_index] = lines[line_index].split()
                try:
                    inputs[type_name[1]].append(float(lines[line_index][0]))
                    inputs[type_name[2]].append(float(lines[line_index][1]))
                    inputs[type_name[3]].append(float(lines[line_index][2]))
                except:
                    print('{} input inproper! check {}'.format(type_name[0],file_name))
                line_index += 1

        else:
            return inputs

class randgen():
    def __init__(self,seed, a = 1103515245, c = 12345 ,m = 32768,interval = (0,1), integer = False, decimals = None):
        #initiation of data input, seed , other LCG parameters, and the interval in which you need the random numbers in.

        self.term = seed
        self.a = a
        self.c = c
        self.m = m
        self.decimals = decimals
        self.integer = integer
        if interval[0] > interval[1]:
            self.interval = (interval[1],interval[0])
        elif interval[0] == interval[1]:
            print('Invalid interval for LCG')
        else:
            self.interval = interval

    def gen(self):
        #generates a random number in the range (0,1)
        self.term = (((self.a * self.term) + self.c) % self.m)
        value = (((self.interval[1]-self.interval[0])*(self.term / self.m)) + self.interval[0])
        if self.integer is True:
            value = int(value)
        if self.decimals is not None:
            value = round(value,self.decimals)
        return value

    def genlist(self,length):
        # returns a list of 'n' random numbers in the range (0,1) where 'n' is 'length'.
        RNs = []
        for i in range(length):
            self.term = (((self.a * self.term) + self.c) % self.m)
            value = (((self.interval[1]-self.interval[0])*(self.term / self.m)) + self.interval[0])
            if self.integer is True:
                value = int(value)
            if self.decimals is not None:
                value = round(value,self.decimals)  
                
            RNs.append(value)
        return RNs

def Randomwalk2D_sim(seed,steps,start = (0,0)):
    #REWRITE
    #simulates Random walk, return list of coordinates and prints a Random walk plot.

    random = randgen(seed)
    # Random_numbers = random.genlist(steps * 2) #need two coordinates
    points = [start] #list of coordinates visited, stored as tuples
    for i in range(steps):
        x = points[i][0] + (2 * random.gen()) - 1
        y = points[i][1] + (2 * random.gen()) - 1
        points.append((x,y))
    
    #graph formatting
    plt.title('Random walk with {} steps'.format(steps))
    plt.axhline(0,lw = 1,c = 'k')
    plt.axvline(0,lw = 1,c = 'k') 
    plt.xlabel('X axis')
    plt.ylabel('Y axis')
    plt.plot([start[0]],[start[1]],'go',ms = 4) #plotting thte start point as green
    plt.plot([points[i][0] for i in range(steps+1)],[points[i][1] for i in range(steps+1)],'-bo',lw = '1',ms = 1, mec = 'k', mfc = 'b' )
    #points[i][0] is the (i)th x cordinate and points[i][1] is the (i)th y cordinate, then a list of each
    plt.plot([points[-1][0]],[points[-1][1]],'ro',ms = 3) #plotting the end point as red
    plt.show()
    return points


def rms_walk(walk):
    #calculating rms distance from a 2D walk

    steps = len(walk) - 1  #n steps = n + 1 coordinates
    sumof_dsquared = 0
    for i in range(1,steps):
        sumof_dsquared += (((walk[i][0] - walk[i-1][0]) ** 2) + ((walk[i][1] - walk[i-1][1]) ** 2))
        # eqn.                  (  x2   -   x1  ) ^ 2         +      (  y2   -   y1  ) ^ 2 
    rms = sqrt(sumof_dsquared / steps)
    return rms


def netdisplace_walk(walk):
    #caluculate the net displacement of walk

    netdis = sqrt(((walk[-1][0] - walk[0][0]) ** 2) + ((walk[-1][1] - walk[0][1]) ** 2))
    #  eqn.          (last x    -   first x) ^ 2    +    (last y    -   first y) ^ 2
    return netdis

def solve_checklinearform(A,B):

    #checking if A is square matrix
    if len(A[0]) != len(A):
        print('A is not square matrix')
        return False

    #checking B is a column matrix
    for i in range(len(B)):
        if len(B[i]) != 1:
            print('B is not a column matrix')
            return False

    #checking if B has same number of elements as rows of A
    if len(A) != len(B):
        print('B is not compatible with A dimensions') 
        return False
    return True

def inv_mat_GJ(A):

    if len(A) != len(A[0]): #if not square matrix, exit
        return None

    n = len(A) #the dimension of matrix will be n*n
    I = []
    for row in range(n):
        I.append(list())
        for col in range(n):
            if col == row:
                I[row].append(1)
            else:
                I[row].append(0)
    
    
    for curr in range(n): #curr takes the index of each column we are processing
        #row swap if zero
        if A[curr][curr] == 0:
            max_row = curr
            for row in range(curr + 1,n):

                if abs(A[row][curr]) > abs(A[max_row][curr]):
                    max_row = row
            if max_row == curr: #if max elemnt is still zero, max_row is not changed; no unique solution
                return None
            A[curr],A[max_row] = A[max_row], A[curr]
            I[curr],I[max_row] = I[max_row], I[curr]
        #making the pivot element 1
        if A[curr][curr] != 1:
            pivot = A[curr][curr]
            for i in range(n):
                A[curr][i] = A[curr][i] / pivot
                I[curr][i] = I[curr][i] / pivot

        #making others zero
        for i in range(0,n):
            if i == curr: #skipping the pivot point
                continue
            if A[i][curr] != 0:
                lead = A[i][curr]
                for j in range(0,n):
                    A[i][j] = A[i][j] - (A[curr][j] * lead)
                    I[i][j] = I[i][j] - (I[curr][j] * lead)

    return I

def solve_GJ(A,B):
    #solves linear equations using Gauss-Jordon elimination method
    #takes the A and B matrix as input from the form A.X = B where X is the unknown matrix
    #returns solved X 

    #checking if matA and matB properly represent linear equation form
    if solve_checklinearform(A,B) == False:
        print('The matrices doesnt properly represent linear eqns form')
        return None
    
    n = len(A) #the dimension is a frequently used value, therefore stored.
    
    #convering matrix A into augmented matrix 
    for row in range(n):
        A[row].append(B[row][0])
    
    for curr in range(n): #curr takes the index of each column we are processing
        #row swap if zero
        if A[curr][curr] == 0:
            max_row = curr
            for row in range(curr + 1,n):

                if abs(A[row][curr]) > abs(A[max_row][curr]):
                    max_row = row
            if max_row == curr: #if max elemnt is still zero, max_row is not changed; no unique solution
                return None
            A[curr],A[max_row] = A[max_row], A[curr]


        #making the pivot element 1
        if A[curr][curr] != 1:
            pivot_term = A[curr][curr]
            for i in range(n + 1):
                A[curr][i] = A[curr][i] / pivot_term
     
        #making others zero
        for i in range(0,n):
            if i == curr: #skipping the pivot point
                continue
            if A[i][curr] != 0:
                main_term = A[i][curr] 
                for j in range(curr,n + 1): #A[curr][j] is zero for j<curr, so no changes to the elements to the left in this row.
                    A[i][j] = A[i][j] - (A[curr][j] * main_term)

    solution = []
    for i in range(n):
        solution.append([A[i][-1]]) #Taking last elements into a list to form column matrix
    return solution

def mat_det(A):
    #matrix determinant using row echelon form
    n = len(A)
    rowswaps = 0 #no. of row swaps done
    if n != len(A[0]):
        print('Not a square matrix')
        return None
    for curr in range(n): #curr takes the index of each column we are processing
        #row swap if zero
        if A[curr][curr] == 0:
            max_row = curr
            for row in range(curr + 1,n):

                if abs(A[row][curr]) > abs(A[max_row][curr]):
                    max_row = row
            if max_row == curr: #if max elemnt is still zero, max_row is not changed; no unique solution
                print('The matrix is singular!')
                return None
            A[curr],A[max_row] = A[max_row], A[curr]
            rowswaps += 1

        #making others elements below lead term zero
        for i in range(curr + 1,n):
            if A[i][curr] != 0:
                main_term = A[i][curr] / A[curr][curr]
                for j in range(curr,len(A[i])): #A[curr][j] is zero for j<curr, so no changes to the elements to the left of curr in this row.
                    A[i][j] = A[i][j] - (A[curr][j] * main_term) 
    prdct = 1
    if rowswaps % 2 == 1:
        prdct = -1
    for i in range(n):
        prdct *= A[i][i]
    return prdct
    
def Chol_dec(A):
    #function that performs the Cholesky decomposition on A and returns (L and T_transpose superimposed matrix where the diagonal of both are same).

    #IF NOT SYMMETRIC, EXIT.
    if check_symmetry(A) == False:
        print('Non symmetric!')
        return None

    n = len(A)

    #Performing the decomposition
    for row in range(n):
        for col in range(row,n): #I am finding the L_transpose and copying elements to the lower part to complete the decomposition
            
            if row == col:
                sum = 0
                for i in range(row):
                    sum += (A[row][i] ** 2)
                A[row][row] = sqrt(A[row][row] - sum)
            else:
                sum = 0
                for i in range(row):
                    sum += (A[row][i] * A[i][col])
                A[row][col] = (A[row][col] - sum) / A[row][row]
                A[col][row] = A[row][col] #As the matrix isresultant L-L_transpose matrix is symmetric

    return A

def forward_backward_cholesky(A,B): 
    # Performs forward substitution and backward substitution given an (L-Ltranspose superimposed) matrix A and and matrix B. returns X.
    if A is None:
        print('Cholesky decomposition failed!')
        return None
    
    #forward substition
    Y = []
    n = len(B)
    for i in range(n):
        sum = 0
        for j in range(i):
            sum += (A[i][j] * Y[j])
        Y.append((B[i][0]-sum) / A[i][i])

    #backward substitution
    X = Y # Here X and Y can be stored in the same matrix due to the way of calculation, for the sake understanding I have made a name X.
          # By assigning a matrix like this, only a another name for the same list is created.
    for i in range(n-1,-1,-1):
        sum = 0
        for j in range(i + 1, n):
            sum += (A[i][j] * X[j])
        X[i] = (Y[i] - sum) / A[i][i]

    for i in range(n):
        X[i] = [X[i]]
    
    return X

def solve_Cholesky(A,B): 
    #solves AX = B using cholesky Decomposition. Returns X.

    #checking if matA and matB properly represent linear equation form
    if solve_checklinearform(A,B) == False:
        print('The matrices doesnt properly represent linear eqns form')
        return None
    
    A = Chol_dec(A)
    if A is None:
        print('Cholesky Solve not possible!')
        return None

    return forward_backward_cholesky(A,B) #Calling forward backward substitution on Cholesky decomposed A and B and returning the solution

def LU_dec(A,B):
    # performs doolittle method LU decompostition of A and return a (L and U superimposed matrix where L_ii = 1 is understood).

    n = len(A) #matrix dimension
    if n != len(A[0]): #checking if the given matrix is a square matrix
        print('Not square!')
        return None
    
    # checking if the determinant of A matrix is 0, if so no unique solution exist
    A1 = mat_copy(A)
    if mat_det(A1) == 0:
        print('Singular matrix!')
        return None
    
    #making leading submatrices non-zero by pivoting
    for i in range(n-1): #last sub matrix is the matarix itself, so avoided.
        
        #constructing submatrix
        submat = []
        for j in range(i+1):
            submat.append(A[j][:i+1])
        
        if mat_det(submat) == 0: #if submat determinant is zero, the last row is replaced with the next rows until the det is non zero. 
            curr = i + 1
            flag = 0
            while curr != n:
                submat[i] = A[curr][:i+1]
                if mat_det(submat) != 0:
                    flag = 1
                    break
                else:
                    curr += 1
            if flag == 0:
                print('Reordering incomplete! submat of order',i+1,'has determinant zero while the matrix itself is non singular.') #this is not supposed to happen, still added in case of any unexpected exception.
                return None
            A[i],A[curr] = A[curr],A[i] #then the original A matrix is changed as per the non-singular submatrix swap.
            B[i],B[curr] = B[curr],B[i] #then the original B matrix is changed as per the non-singular submatrix swap.

    #performing LU decompostion
    for j in range(n):
        for i in range(1,n):
            if i <= j:
                sum = 0
                for k in range(0,i):
                    sum += (A[i][k] * A[k][j])
                A[i][j] = A[i][j] - sum
            else:
                sum = 0
                for k in range(0,j):
                    sum += (A[i][k] * A[k][j])
                A[i][j] = (A[i][j] - sum) / A[j][j]
    return A,B

    
def forward_backward_LU(A,B):

    if A is None: #checking if LU decomposition failed
        print('LU decomposition failed!')
        return None

    Y = [] # the intermidiate Y matrix where Y = U.X .thus Y is also be in equation L.Y = B, we solve this first
    n = len(B)
    for i in range(n):
        sum = 0
        for j in range(i):
            sum += (A[i][j] * Y[j])
        Y.append(B[i][0]-sum)
    
    #now we solve U.X = Y
    X = Y
    for i in range(n-1,-1,-1):
        sum = 0
        for j in range(i + 1, n):
            sum += (A[i][j] * X[j])
        X[i] = (Y[i] - sum) / A[i][i]

    #making X into column matrix
    for i in range(n):
        X[i] = [X[i]]

    return X

def LU_solve(A,B):
    #solves the set of linear eqn.s in the form of AX = B using LU decomposition method

    #checking if matA and matB properly represent linear equation form
    if solve_checklinearform(A,B) == False:
        print('The matrices doesnt properly represent linear eqns form')
        return None

    #PErforming the LU decomposition
    A,B = LU_dec(A,B)

    return forward_backward_LU(A,B) #Calling forward backward substitution of LU on LU decomposed A and B and returning the solution

def mat_makeSDD(A,B):
    #makes a strictly diagonally dominant matrix A after row swaps on A and the same row operations on B. returns both.

    n = len(A)
    

    for i in range(n):
        sum = 0
        for j in range(n):
            if j != i:
                sum += abs(A[i][j])
        if abs(A[i][i]) > sum: #checkis if leading element is strictly dominant over others in the row. if so, continue to the next row.
            # print('row',i,'is fine')
            continue
        else: #else, check for the same position dominance in next rows.
            curr = i + 1
            flag = 0
            while curr < n:
                sum = 0
                for j in range(n):
                    if j != i:
                        sum += abs(A[curr][j])
                # print('sum is',sum)
                if abs(A[curr][i]) > sum:
                    flag = 1
                    break
                else:
                    curr += 1
            if flag == 0: #if no dominant term was found in the rows below, it cannot be made strictly diagonally dominant matrix.
                # print('failed at',i)          #it is so, from the fact that a row/column cannot have 2 diagonally dominant terms, no row/column swap can separate those.
                return None,None
            else:
                A[i],A[curr] = A[curr],A[i]
                B[i],B[curr] = B[curr],B[i]
    return A,B

def solve_GS(A,B,tolerance):
    #solving linear equations (A.X=B) using Gauss-seidel method with tolerace, where X is the solution. return solution X and number of steps in which it is done
    
    #checking if matA and matB properly represent linear equation form
    if solve_checklinearform(A,B) == False:
        print('The matrices doesnt properly represent linear eqns form')
        return None    
    
    n = len(A)

    #trying to make strictly diagonally dominant , now DISABLED
    # A,B= mat_makeSDD(A,B)   
    # if A is None:
    #     print('matrix cannot be made strictly diagonally dominant\nGauss-Seidel method exited!')
    

    X = list([0] for i in range(n)) # initial guess
    
    for steps in range(100):
        flag = 1
        for i in range(n):
            sum = 0
            for j in range(i):
                sum += (A[i][j] * X[j][0])
            for j in range(i+1,n):
                sum += (A[i][j] * X[j][0])
            temp = (B[i][0] - sum) / (A[i][i])
            if abs((temp) - (X[i][0])) > tolerance: #checks the tolerance at each new value
                flag = 0
            X[i][0] = temp
        if flag == 1:
            return X,steps + 1
    print('Eqn not solved after 100 steps')
    return None,100

def solve_Jacobi(A,B,tolerance):

    #checking if matA and matB properly represent linear equation form
    if solve_checklinearform(A,B) == False:
        print('The matrices doesnt properly represent linear eqns form')
        return None

    #trying to make strictly diagonally dominant
    A,B= mat_makeSDD(A,B)
    if A is None:
        print('matrix cannot be made strictly diagonally dominant\nJacobi method exited!')

    n = len(A)

    currX = [0] * n # initial guess
    newX = currX[:] # new guess, to be processed

    steps = 0
    while steps < 150: #max steps kept high for 10^-6 precision
        
        #finding new X
        for i in range(n):
            sum = 0
            for j in range(n):
                if i != j: sum += (A[i][j] * currX[j])
            newX[i] = (B[i][0] - sum) / A[i][i]

        #comparing newX and currentX against tolerance
        flag = 1 #assumed true
        for i in range(n):
            if abs(newX[i]-currX[i]) > tolerance:
                flag = 0 #looking for exception
                steps += 1
                currX = newX[:]
                break
        
        if flag == 1:
            #converting the answer to column matrix
            for i in range(n):
                newX[i] = [newX[i]]
            return newX, steps
    #if not done in 150 steps, return None
    print('After more than 150 iterations, it is not solved!')
    return None, steps

def root_bracketing(f,a,b):

    #makes the interval 1.5 times bigger and shifts to closer to zero side for 12 times
    for i in range(12):
        if f(a) * f(b) < 0:
            return a,b
        if abs(f(a)) < abs(f(b)):
            temp = a
            a = a - (1.5*(b-a))
            b = temp
        else:
            temp = b
            b = b + (1.5*(b-a))
            a = temp

    #after 12 times, if not bracketed, the process restarts with 2.5 times bigger including the current interval
    if abs(f(a)) < abs(f(b)):
        return root_bracketing(a - (1.5*(b-a)),b)
    else:
        return root_bracketing(a,b + 1.5*(abs(b-a)))


def root_bisection(f,a,b,epsilon,delta,details = False):

    # bracketing the root if not already done
    a,b = root_bracketing(f,a,b)

    steps = 1

    #if details kept true, it will return a list of the values converging to the required root
    if details is True: 
        tables = []
        while steps<100:
            if (a-b) < epsilon:
                if abs(f(a)) < delta or abs(f(b)) < delta:
                    return tables
            c = (a + b) / 2
            tables.append(c)
            if f(a) * f(c) < 0:
                b = c
            else:
                a = c
            steps += 1
        print('Not converging after 100 steps, terminated')
        return None

    #if details kept False (default), it will return the root and the number of steps taken
    else: 
        while steps<100:
            if (a-b) < epsilon:
                if abs(f(a)) < delta or abs(f(b)) < delta:
                    return (a+b)/2,steps
            c = (a + b) / 2
            if f(a) * f(c) < 0:
                b = c
            else:
                a = c
            steps += 1
        print('Not converging after 100 steps, terminated')
        return None,None

def root_newtonraphson(f,df,x0,epsilon,delta,details = False):
    #returns the root for function f from guess x0 with epsilon tolerance for closeness to root and delta tolerance for f's closeness to zero at that root
    
    if details is True:
        tables = []
        steps = 1
        while steps < 151: #150 steps as max. limit

            #new guess
            x1 = x0 - (f(x0)/df(x0))
            tables.append(x1)
            #checking is guess is good enough
            if abs(x1 - x0) < epsilon and f(x1) < delta:
                return tables

            #preparing for next iteration
            x0 = x1
            steps += 1
        
        print('Not converging after 151 steps, terminated')
        return None
    else:
        steps = 1
        while steps < 151: #150 steps as max. limit

            #new guess
            x1 = x0 - (f(x0)/df(x0))

            #checking is guess is good enough
            if abs(x1 - x0) < epsilon and f(x1) < delta:
                return x1,steps

            #preparing for next iteration
            x0 = x1
            steps += 1
        
        print('Not converging after 151 steps, terminated')
        return None,None

def root_regulafalsi(f,a,b,epsilon,delta,details = False):
    
    #bracketing the root if not done
    a,b = root_bracketing(f,a,b)

    # doing the first run, as no comparison is done here
    c = b - (((b-a)*f(b))/(f(b) - f(a)))
    if f(c)*f(a) < 0:
        b = c
    else:
        a = c

    #if details kept true, it will return a list of the values converging to the required root
    if details is True:
        #adding the c from first step
        tables = []
        tables.append(c)

        steps = 2
        while steps < 100:
            oldc = c
            c = b - (((b-a)*f(b))/(f(b) - f(a)))
            tables.append(c)
            if f(c)*f(a) < 0:
                b = c
            else:
                a = c
            if abs(c-oldc) < epsilon and f(c) < delta:
                return tables
            oldc = c
            steps += 1
        print('Not converging after 100 steps, terminated')
        return None,None
    
    #if details kept False (default), it will return the root and the number of steps taken
    else:
        steps = 2
        while steps < 100:
            oldc = c
            c = b - (((b-a)*f(b))/(f(b) - f(a)))
            if f(c)*f(a) < 0:
                b = c
            else:
                a = c
            if abs(c-oldc) < epsilon and f(c) < delta:
                return c,steps
            oldc = c
            steps += 1
        print('Not converging after 100 steps, terminated')
        return None,None



def lagrange_intrapolate(X,Y,x1):
    N = len(X)
    if N != len(Y):
        print('data mismatch')
        return None
    sum = 0
    for i in range(N):
        prdct = 1
        for k in range(N):
            if k == i:
                continue
            prdct = prdct * ((x1 - X[k])/(X[i]-X[k]))
        sum += prdct * Y[i]
    return sum


def px_deflate(px, root):
    # here px is in the format a0,a1,a2,..an where n is the polynomial power.
    n = len(px)
    #Here, the passed 'root' is a verified root from main body with a certain tolerance, therefore checking not done.
    if n == 1:
        print('P(x) doesnt contain any x: deflation exited!')
        return None
    n -= 1
    if px[n] != 1:
        lead = px[n]
        for i in range(len(px)):
            px[i] = px[i] / lead
    n -= 1
    while n >= 0:
        px[n] = (px[n+1] * root) + px[n]
        n -= 1
    return px[1:]



def px_derivative(px):
        # here px is in the format a0,a1,a2,..an where n is the polynomial power.

    n = len(px)
    i = 1
    while i != n:
        px[i] = px[i] * (i)
        i += 1
    return px[1:]

def px_value(px,x):
        # here px is in the format a0,a1,a2,..an where n is the polynomial power.

    sum = 0
    i = 0
    n = len(px)
    while i != n:
        sum+= (px[i] * (x**i))
        i += 1
    return sum

def root_laguire(px,guess,tolerance):
        # here px is in the format a0,a1,a2,..an where n is the polynomial power.

    roots = []
    steps = 1
    N = len(px) - 1
    n = N
    while steps <= N:
        if abs(px_value(px,guess)) < tolerance:
            print('{}th root, {:.4f}, is found in 0 steps.'.format(steps,guess))
            steps += 1
            roots.append(guess)
            px = px_deflate(px,guess)
            n -= 1
            continue
        # print('Finding the {}th root for {}'.format(steps,px))
        dpx = px_derivative(px[:])
        ddpx = px_derivative(dpx[:])
        i = 1
        theguess = guess
        while True:
            g = px_value(dpx,theguess)/px_value(px,theguess)
            h = (g**2) - (px_value(ddpx,theguess)/px_value(px,theguess))
            if g < 0:
                a = (n / (g - m.sqrt((n-1)*((n*h)-(g**2)))))
            else:
                a = (n / (g + m.sqrt((n-1)*((n*h)-(g**2)))))
            # print('a is',a)
            newguess = theguess - a
            # print(px_value(px,theguess),'and',px_value(px, newguess))
            
            if i < 26 :
                # if abs(a) < tolerance and px_value(px,newguess) < tolerance:
                if px_value(px,newguess) < tolerance:
                    print('{}th root, {:.4f}, is found in {} steps.'.format(steps,newguess,i))
                    steps += 1
                    roots.append(newguess)
                    px = px_deflate(px,newguess)
                    n -= 1
                    break
                # else:
                #     print('Guess discarded.\n')

            else:
                print('The guess for {} th root was not found in 25 steps'.format(steps))
                return None
            theguess = newguess
            i += 1
    return roots

def fit_polyleastsq(dataX,dataY,order,datasigma=None):


    n = len(dataX) #no. of datasets

    if len(dataY) != n:
        print("Data mismatch! exited!")
        return None

    matX = [[0]*(order+1) for i in range(order+1)]

    #completing matX
    matX[0][0] = n
    for i in range(1,(2 * order) + 1):
        sum = 0
        for j in range(n):
            sum += (dataX[j]**i)
        if i <= order:
            startX = 0
            startY = i
            while startY >= 0:
                # print('p',startX,startY)

                matX[startX][startY] = sum
                startY -= 1
                startX += 1
        else:
            
            startX = i - order
            startY = order
            while startX <= order:
                # print('p',startX,startY)
                matX[startX][startY] = sum
                startY -= 1
                startX += 1

    matY = []
    sum = 0
    for i in range(n):
        sum += dataY[i]
    matY.append([sum])
    for i in range(1,order+1):
        sum = 0
        for j in range(n):
            sum += (dataX[j]**i)*dataY[j]
        matY.append([sum])

    # print(matX,matY)
    px = solve_GJ(matX,matY)
    for j in range(len(px)):
        px[j] = px[j][0]

    return px

def fit_linearleastsq(dataX,dataY,datasigma=None):
    #returns the linear equation in polynomial form where the line is a0+a1x

    n = len(dataX) #no. of datasets

    if len(dataY) != n:
        print("Data mismatch! exited!")
        return None
    Sxx = 0
    Syy = 0 #for pearson R square calc.
    Sxy = 0
    Sx = 0
    Sy = 0

    if datasigma is None:
        S = n
        for i in range(n):
            Sxx += dataX[i]**2
            Syy += dataY[i]**2
            Sxy += dataX[i] * dataY[i]
            Sx += dataX[i]
            Sy += dataY[i]

    else:
        S = 0
        for i in range(n):
            S += 1/(datasigma[i]**2)
            Sxx += (dataX[i]**2)/(datasigma[i]**2)
            Syy += (dataY[i]**2)/(datasigma[i]**2)
            Sxy += (dataX[i] * dataY[i])/(datasigma[i]**2)
            Sx += dataX[i]/(datasigma[i]**2)
            Sy += dataY[i]/(datasigma[i]**2)

    delta = (S*Sxx)-(Sx**2)
    a0 = ((Sxx*Sy) - (Sx*Sxy)) / delta
    a1 = ((Sxy*S) - (Sx*Sy)) / delta

    R2 = (((n*Sxy) - (Sx*Sy))**2)/(((n*Sxx)-(Sx**2)) * ((n*Syy)-(Sy**2)))

    return [a0,a1],R2


def px_graphdata(px,start,stop,number):
    # px in the format [x0,x1,x2...]
    
    step = (stop - start) / (number - 1)
    Xvalues = []
    Yvalues = []
    stop += step
    while start < stop:
        Xvalues.append(start)
        Yvalues.append(px_value(px,start))
        start += step
    return Xvalues, Yvalues

def integrate_midpoint(f,a,b,N):
    # integrates f over a to b by dividing to N parts in midpoint method
    #getting the equation 'h * [ sum(f(x)) for all x ~ middle of subsections ] '

    if N < 1:
        N = 10 #default value if N is not a viable value!
        print("N not viable, took default value 10.")

    if N % 1 != 0:
        N = N // 1 #to get integer if it is not

    step = (b - a) / N
    x = a + (step / 2) #reaching middle of a section and jumping steps
    
    sum = 0
    while x < b:
        sum += (f(x))
        x += step
    sum *= step

    return sum
            
def integrate_trapezoidal(f,a,b,N):
    #integrates f over a to b by dividing to N parts in trapezoidal method
    #getting the equation 'h * [ (f(a)+f(b))/2 + ( sum(f(x)) for all x subsections ) ] '
    
    if N < 1:
        N = 10 #default value if N is not a viable value!
        print("N not viable, took default value 10.")

    if N % 1 != 0:
        N = N // 1 #to get integer if it is not

    step = (b-a)/N
    sum = (f(a)+f(b))/2
    b -= (step/2) #to ensure that even if rounding error accumulates in using the step multiple times, we stop before the last value
    a += step
    while a < b:
        sum += f(a)
        a += step
    sum *= step
    return sum

def integrate_simpson(f,a,b,N):
    #integrates f over a to b by dividing to N parts in simpson method
    #getting the equation '(step/3) * [ f(a) + f(b) + ( 4*sum(f(x)) for all x ~ middle of subsections ) + ( 2*sum(f(x)) for all x subsections ) ] '

    if N < 1:
        N = 5 #default value if N is not a viable value!
        print("N not viable, took default value 5.")

    if N % 1 != 0:
        N = N // 1 #to get integer if it is not

    h = (b-a)/N
    step = h/2
    sum = f(a)+f(b)+(4*f(a+(step)))
    b -= step # for stopping before the last value if any rounding errors accumulate
    a += h

    while a < b:
        sum += ((2*f(a)) + (4*f(a + step)))
        a += h
    sum *= (step/3)
    return sum

def integrate_montecarlo(f,a,b,N):
    # does monte-carlo integration by taking N random numbers in the a to b interval and doing a sum over the function applied to those random numbers

    if N < 1:
        N = 10 #default value if N is not a viable value!
        print("N not viable, took default value 10.")

    if N % 1 != 0:
        N = N // 1 #to get integer if it is not.

    rand = randgen(seed = 1234,interval = (a,b))
    sum = 0
    for i in range(N):
        sum += f(rand.gen())
    return sum*((b-a)/N)

def fit_powerlaw(dataX,dataY):
    # fits the data given by dataX and dataY using (a*x^b) model
    # returns a,b and pearsons R square
    
    n = len(dataX) #no. of datasets

    #checking the data mismatch
    if len(dataY) != n:
        print("Data mismatch! exited!")
        return None

    #convering to linear form data
    logxdata = []
    logydata = []
    for i in range(n):
        logxdata.append(m.log(dataX[i]))
        logydata.append(m.log(dataY[i]))
    
    
        #using least square fit for linear form
    px,prs = fit_linearleastsq(logxdata,logydata)

    #calculating the a and b from the linear polynomial solution
    a = m.exp(px[0])
    b = px[1]

    return a,b,prs

def fit_exponential(dataX,dataY):
    # fits the data given by dataX and dataY using (a*e^bx) model
    # returns a,b and pearsons R square

    n = len(dataX) #no. of datasets

    #checking data mismatch
    if len(dataY) != n:
        print("Data mismatch! exited!")
        return None

    #converting to linear form
    logydata = []
    for i in range(n):
        logydata.append(m.log(dataY[i]))
    
    #using least square fit for linear form
    px,prs = fit_linearleastsq(dataX,logydata) #prs is pearson r square

    #calculating the a and b from the linear polynomial solution
    a = m.exp(px[0])
    b = px[1]

    return a,b,prs

def powerlaw_data(a,b,start,stop,number):
    # gives the graph data for y = a*x^b between start and stop with "number" no. of points
    step = (stop - start) / (number - 1)
    Xvalues = []
    Yvalues = []
    stop += step
    while start < stop:
        Xvalues.append(start)
        Yvalues.append(a*(start**b))
        start += step
    return Xvalues, Yvalues


def exp_graphdata(a,b,start,stop,number):
    # gives the graph data for y = a*e^bx between start and stop with "number" no. of points
    step = (stop - start) / (number - 1)
    Xvalues = []
    Yvalues = []
    stop += step
    while start < stop:
        Xvalues.append(start)
        Yvalues.append(a*(m.exp(start*b)))
        start += step
    return Xvalues, Yvalues

def fx_graphdata(f,a,b,number):
    # gives the graph data for f(x) from a to b and gives n number of sections.
    step = (b-a)/number
    datX = []
    datY = []
    for i in range(number+1):
        datX.append(a)
        datY.append(f(a))
        a += step
    return datX,datY

def diff_eulerforward(dydx,x0,y0,x1,dx):
    # integrates and give data points between x0 and x1 for the solution of a given dydx
    # x0 y0 is the given intial solution, dx is the delta X  used in each iteration 

    x0
    datX = [x0]
    datY = [y0]
    while x0 < x1:
        y0 = y0 + dx*(dydx(y0,x0))
        x0 += dx
        datX.append(x0)
        datY.append(y0)

    return datX, datY

def diff_predictorcorrector(dydx,x0,y0,x1,dx):
    # integrates and give data points between x0 and x1 for the solution of a given dydx
    # x0 y0 is the given intial solution, step is the delta X  used in each iteration 

    start = x0
    datX = [x0]
    datY = [y0]

    while x0 < x1:
        k1 = dx*dydx(y0, x0)
        k2 = dx*dydx(y0 + k1, x0 + dx)

        y0 = y0 + (k1+k2)/2
        x0 = x0 + dx
        datX.append(x0)
        datY.append(y0)

    return datX, datY

def diff_RK2(dydx,x0,y0,x1,dx):
    # integrates and give data points between x0 and x1 for the solution of a given dydx
    # x0 y0 is the given intial solution, step is the delta X  used in each iteration 

    datX = [x0]
    datY = [y0]
    while x0 < x1:
        k1 = dx * dydx(y0,x0)
        k2 = dx * dydx((y0 + (k1/2)),(x0 + (dx/2)))
        y0 = y0 + k2
        x0 += dx
        datX.append(x0)
        datY.append(y0)
        
    return datX, datY

def diff_RK4(dydx,x0,y0,x1,dx):
    # integrates and give data points between x0 and x1 for the solution of a given dydx
    # x0 y0 is the given intial solution, step is the delta X  used in each iteration 

    datX = [x0]
    datY = [y0]
    while x0 < x1:
        k1 = dx * dydx(y0,x0)
        k2 = dx * dydx((y0 + (k1/2)),(x0 + (dx/2)))
        k3 = dx * dydx((y0 + (k2/2)),(x0 + (dx/2)))
        k4 = dx * dydx((y0 + k3),(x0 + dx))
        y0 = y0 + ((k1 + k4 + (2*(k2 + k3)))/6)
        x0 += dx
        datX.append(x0)
        datY.append(y0)
        
    
    return datX, datY


def RK4_coupled(dydxlist,x0,y0list,x1,dx):
    # Using RK4 to solve coupled ODE equations.
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

def diff_shooting(flist,x0,y0,x1,y1,dy0guess1,tolerance,step):
    # does shooting method using RK4
    # return X values and Y values solution for graph
    
    # running initial dy0guess
    datX,datY = RK4_coupled(flist,x0,[y0,dy0guess1],x1,step)
    yval1 = datY[0][-1]
    if abs(yval1 - y1) < tolerance:
        return datX,datY
    if yval1 < y1:
        guess1side = -1
    else :
        guess1side = 1

    # getting second dy0guess
    dy0guess2 = dy0guess1 + 1   
    yval2 = RK4_coupled(flist,x0,[y0,dy0guess2],x1,step)[1][0][-1]
    if yval2 < y1:
        guess2side = -1
    else :
        guess2side = 1
    while guess1side * guess2side != -1:
        #checking if the second guess brackets the needed solution guess

        #appropriate movement of guess
        if abs(y1-yval2) > abs(y1-yval1):
            dy0guess2 = dy0guess1 - 1.5*abs(dy0guess2-dy0guess1)
        else:
            dy0guess2 += abs(dy0guess2-dy0guess1)
        
        yval2 = RK4_coupled(flist,x0,[y0,dy0guess2],x1,step)[1][0][-1]
        if yval2 < y1:
            guess2side = -1
        else :
            guess2side = 1
    
    #arranging so that guess1 is leftside one and guess2 is rightside
    if guess1side == 1:
        dy0guess1,dy0guess2 = dy0guess2,dy0guess1

    #doing lagrange interpolation and checking solution till the needed value is within tolerance
    i = 0
    while True:
        newguess = dy0guess1 + (((dy0guess2 - dy0guess1)/(yval1 - yval2))*(y1 - yval2))
        i += 1
        datY = RK4_coupled(flist,x0,[y0,newguess],x1,step)[1][0]
        yvalnew = datY[-1]

        if abs(yvalnew - y1) < tolerance:
            break
        if yvalnew < y1:
            dy0guess1 = newguess
            yval1 = yvalnew
        else:
            dy0guess2 = newguess
            yval2 = yvalnew
    return datX, datY

def graphsearch(Xlist,Ylist,x0 = None, y0 = None):
    # given the data points for a graph while assuming close enough points, we find x0 for a given y0 or vice versa
    # done using lagrange intrapolate by taking a pair of x values and y values between which the x0 or y0 lies

    # if y0 is given
    if x0 is None:
        if Ylist[0] <= y0:
            i = 1
            while Ylist[i] <= y0:
                i += 1
        else:
            i = 1
            while Ylist[i] > y0:
                i += 1
        return ( Xlist[i-1] + ((Xlist[i]-Xlist[i-1])*((y0-Ylist[i-1])/(Ylist[i] - Ylist[i-1]))) )

    # if x0 is given
    elif y0 is None:
        if Xlist[0] <= x0:
            i = 1
            while Xlist[i] <= x0:
                i += 1
        else:
            i = 1
            while Xlist[i] > x0:
                i += 1
        return ( Ylist[i-1] + ((x0 - Xlist[i-1])*((Ylist[i]-Ylist[i-1])/(Xlist[i]-Xlist[i-1]))) )
    
    #anything else
    else:
        print("invalid data to search!")
        return None

def eigen_rayleigh(A,x0,tolerance):
    # finding dominant eigen value and corresponding eigenvector for A matrix
    # x0 is the intial guess for eigenvector and the tolerance for this fnction is tolerance.
    # return eigen value and normalized eigenvector

    #initial run assumed not close enough
    Akx0 = (A*x0)
    print(Akx0)
    denom = (Akx0.dot(x0))
    Akx0 = (A*Akx0)
    num = (Akx0.dot(x0))
    k = 0
    lamk = num/denom
    
    # iterates till eigenvalue is close enough
    while k<50:
        denom = num
        evec = Akx0.copy()
        Akx0 = (A*Akx0)
        num = (Akx0.dot(x0))
        k += 1
        lamk2 = num / denom
        if abs(lamk - lamk2) < tolerance:
            evec *= (1/abs(evec))
            return lamk2,evec,k
        lamk = lamk2
    print('Eigenvalue not found within 50 steps')

def partialdiff_heatdiffusion(f0,xi,xf,ti,tf,nx,nt,times = []):
    # partial differential solver for heat diffusion equation.
    # f0 is the intial temperature function at time = 0
    # xi, xf and ti, tf are the initial and final x and t respectively.
    # nx is the number x divsions and nt is the number of t divisons
    # times is the list of time steps needed on the graphdata
    # returns the graphdata

    hx = (xf-xi)/nx
    ht = (tf-ti)/nt
    alpha = ht/(hx**2)

    # initial V0 data
    V0 = [f0(xi)]
    X = [xi]
    x = xi + hx
    for i in range(nx):
        X.append(x)
        V0.append(f0(x))
        x += hx
    Vlist = [V0[:]]
    # bulding V1 for each time step
    for i in range(1,nt+1):
        V1 = [(alpha*V0[1] + ((1-(2*alpha))*V0[0]))]
        for j in range(1,nx):
            V1.append((alpha*(V0[j+1]+V0[j-1]) + ((1-(2*alpha))*V0[j])))
        V1.append(alpha*V0[nx-1] + ((1-(2*alpha))*V0[nx]))

        # offset = V1[0] - V0[0]
        # for j in range(nx+1):
        #     V1[j] -= offset

        V0 = V1[:]

        if i in times:
            Vlist.append(V0[:])
            
    if bool(times) is False:
        return X,V0
    return X,Vlist

def radioactivity(Nlist,tlist,step,tf,seed1):
    # Does the action of radioactivity on N sets of particle whose initial numbers are in list Nlist
    # tlist is their given half lives of each set of particles

    #initialization
    n = len(Nlist)
    Nlist.append(0)
    result = Nlist[:]
    for i in range(len(result)):
        result[i] = [result[i]]
    for i in range(n):
        tlist[i] = ((m.log(2)/tlist[i]) * step) #converting half lives into probablities of decay for each
    ran = randgen(seed = seed1)
    t = step
    times = [0]
    decays = [[] for i in range(n)]

    #determining if decay occurs and add to data
    while t <= tf:
        for i in range(n):
            trials = Nlist[i]
            decay = 0
            for j in range(trials):
                if ran.gen() <= tlist[i]:
                    # print("decayed")
                    decay += 1
            Nlist[i] -= decay
            Nlist[i+1] += decay
            decays[i].append(decay)
        times.append(t)
        for i in range(len(result)):
            result[i].append(Nlist[i])
        t += step
    return times,result,decays

def invmat_LU(A):
    #inverse of matrix using LU
    # takes matrix in nested list form

    if len(A) !=  len(A[0]):
        print("Not square matrix, program exited!")
        return None

    n = len(A) #dimension

    # identity matrix column wise
    Icols = [0 for i in range(n)]
    for i in range(n):
        Icols[i] = [[0] for i in range(n)]
        Icols[i][i][0] = 1
    # print(Icols)

    results = [[] for i in range(n)]
    for i in range(n):
        # solves each column of identity matrix
        A1 = mat_copy(A)
        col = LU_solve(A1,Icols[i])
        if col is None:
            print("unexpected eror, inverese doesnt exist!")
            return None
        for j in range(n):
            results[j].append(col[j][0])
    
    return results

        
def wallhole(Nl,step,tf,seed1):
    # Does the action of radioactivity on N sets of particle whose initial numbers are in list Nlist
    # tlist is their given half lives of each set of particles

    #initialization\
    N0 = Nl
    Nlist = [Nl,0]
    result = Nlist[:]

    for i in range(len(result)):
        result[i] = [result[i]]
    
    ran = randgen(seed = seed1)
    t = step
    times = [0]

    #determining if decay occurs and add to data
    while t <= tf:
        trials = Nlist[0]
        prob = trials/N0
        # print(prob)
        if ran.gen() < prob:
            Nlist[0] -= 1
            Nlist[1] += 1
        else:
            Nlist[0] += 1
            Nlist[1] -= 1
        
        times.append(t)
        for i in range(len(result)):
            result[i].append(Nlist[i])
        t += step
    return times,result


    