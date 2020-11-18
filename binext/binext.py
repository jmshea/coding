#!/usr/bin/env python
# coding: utf-8
'''
John M. Shea
10/22/20
v 0.6

Opinionated Python module to implement Galois field arithmetic for binary extension fields GF(2^m)

Offers two classes:

ff: Represents a full finite field 
   
    Initializer takes a field specifier, which can be either:
	field size (integer, power of 2)
	primitive polynomial, given as list of exponents with nonzero coefficient

    Primary purpose is to set up arithmetic tables for the ffelt class described below

    Offers 1 method:

    print_minpoly_table():
	prints out the minimum polynomials associated with this field, along with the 
        minimum power of alpha that is a root of the minimum polynomial. Produces a
        table of the same form as in Lin and Costello Appendix B

ffelt: Represents an element in a finite field

    Initializer takes two arguments, a power of the primitive element and a field specifier (see ff above)
    Elements are entered and returned as powers of the primitive element.

    Offers methods for adding, multiplying, and raising elements to powers

    Also offers:

    vec(): return the vector representation of an element
    conjugates(): returns a list of conjugates for the element
    minpoly(): returns the minimum polynomial (in GF[2]) for the element
   
'''
import numpy as np
import copy
import pandas

class ff:
    '''
    John M. Shea
    10/22/20
    
    Generate lookup tables for operations in a binary extension field
    
    Must pass it a specifier for the field, which may be either:
    - the order of the field (uses a default primitive polynomial)
    - a primitive polynomial for the field
    '''    
    
    primpolys={4:[0,1,2],      # x^2+x+1
               8:[0,1,3],      # x^3+x+1
               16:[0,1,4],     # x^4+x+1
               32:[0,2,5],     # x^5+x^2+1
               64:[0,1,6],     # x^6+x+1 
               128:[0,3,7],    # x^7+x^3+1
               256:[0,2,3,4,8] # x^8+x^4+x^3+x^2+1
              } # To do: add more primitive polys  (from App. B of Lin and Costello)
    

    
    def __init__(self, gspec, debug=False):
        
        if type(gspec)==int:
            if gspec in self.primpolys.keys():
                glist=self.primpolys[gspec]
            else:
                raise "No default primitive polynomial for field order " +str(gspec)
        elif type(gspec)==list:
            glist=gspec
        else:
            raise "Must specify g as either field order or list of nonzero primitive poly coeffs"
            
        # Now convert g in list form to an integer
        self.g=0
        #print(glist)
        for i in glist:
            self.g+=2**i
            #print(i,g)
        self.m=int(np.log2(self.g))
        self.q=2**self.m
        self.fmt="{0:0"+str(self.m)+"b}"
        if debug:
            print(self.m,self.q)
        
        # Now find the vector representations and store in forward and reverse dicts
        gp=self.g%self.q
        a=1
        self.power_to_poly={None:0}
        self.poly_to_power={0:None}
        for n in range(self.q-1):
            self.power_to_poly[n]=a
            self.poly_to_power[a]=n
            if debug:
                print("{0:3} {1:3}  ".format(n,a), self.fmt.format(a))
            a=a<<1
            #print(a)

            if a&self.q>0:
                a=(a&(self.q-1)) ^ gp
                
    def print_minpoly_table(self):
        a=ffelt(1, self)
        seen=[]
        sumdeg=0
        printed=0
        #print("{0:3}: {1:30}".format(0, str((a**0).minpoly())),end="")
        #printed=1
        for i in range(1, self.q-2,2):
            mp=(a**i).minpoly(True)
            #print((a**i), mp)
            if mp not in seen:
                if printed%2==0:
                    print("{0:3}: {1:30}".format(i, str(mp)),end="")
                else:
                    print("{0:3}: {1:30}".format(i, str(mp)))
                seen+=[mp]
                printed+=1
                sumdeg+=(len(mp)-1)
                if sumdeg>=self.q-2:
                    break

    def add_table(self):
        num_additions = np.square(self.q-1)

        # generate basic element of the field
        a = ffelt(1, self, suppressField=True)

        # generate all unique elements
        elements = []
        for i in range(self.q-1):
            elements.append(a**i)

        #print(elements)

        # Now build the 2D array of additions:

        sum_array = np.zeros((self.q-1, self.q-1))

        sum_list = []
        for i, a_i in enumerate(elements):
            inner_list = []
            for j, a_j in enumerate(elements):

                #sum_array[i,j] = a_i + a_j
                inner_list.append(a_i + a_j)

            sum_list.append(inner_list)

        element_strs = []
        for e in elements:
            element_strs.append(str(e))

        GF_addition_table = pandas.DataFrame(sum_list, index=element_strs, 
                                            columns=element_strs)
    
        #print(GF_addition_table)

        return GF_addition_table
    

    def mul_table(self):
        num_additions = np.square(self.q-1)

        # generate basic element of the field
        a = ffelt(1,self,suppressField=True)

        # generate all unique elements
        elements = []
        for i in range(self.q-1):
            elements.append(a**i)

        #print(elements)

        # Now build the 2D array of products:

        prod_array = np.zeros((self.q-1, self.q-1))

        prod_list = []

        for i, a_i in enumerate(elements):
            inner_list = []
            for j, a_j in enumerate(elements):
                inner_list.append(a_i * a_j)

            prod_list.append(inner_list)


        # Now make the pretty table using pandas:
        # make string column/row labels
        element_strs = []
        for e in elements:
            element_strs.append(str(e))

        GF_product_table = pandas.DataFrame(prod_list, index=element_strs, 
                                            columns=element_strs)

        #print(GF_addition_table)

        return GF_product_table  

'''Helper function needed by ffelt class'''
def bin_array(num, m):
    """Convert a positive integer num into an m-bit bit vector"""
    return np.array(list(np.binary_repr(num).zfill(m))).astype(np.int8)



class ffelt:
    '''
    John M. Shea
    10/22/20

    arguments:
    elt: a power of alpha (optional, will default to 1 if not given)
    f: a field specifier (required). See below for details.

    You must either pass it a ff object OR a specifier of the field,
    which can either be the order of the field (to use the default poly)
    or a list of the nonzero coefficients of the primitive polynomial of the field
    '''

    def __init__(self, elt, f=-1, debug=False, suppressField=False):
        #print(f, type(f), type(f)==ff)

        if f==-1:
            f=elt
            elt=1

        if type(f)==ff:
            self.f=f
        else:
            self.f=ff(f)

        self.q=self.f.q
        if elt == None:
            self.elt = None
        elif type(elt) == int:
            self.elt=elt % (self.q-1)
        else:
            raise "elt must be int or None, not "+str(type(elt))

        self.debug=debug
        self.suppressField=suppressField

    def __add__(self, a):

        if a==None: 
            result=self.elt

        if a==1:
            result=self.f.power_to_poly[self.elt]^self.f.power_to_poly[0]
                    
            if self.debug:
                print(self.f.fmt.format(self.f.power_to_poly[self.elt]))        
                print(self.f.fmt.format(self.f.power_to_poly[a]))
                print(self.f.fmt.format(result))
                print(self.f.poly_to_power[result])
        elif type(a)==ffelt:
            if a.q==self.q:
                result=self.f.power_to_poly[self.elt]^a.f.power_to_poly[a.elt]
                if self.debug:
                    print(self.f.fmt.format(self.f.power_to_poly[self.elt]))        
                    print(self.f.fmt.format(self.f.power_to_poly[a.elt]))
                    print(self.f.fmt.format(result))
                    print(self.f.poly_to_power[result])

            else:
                print("element field size:", self.q)
                print("argument field size:", a.q)
                raise "Cannot add elements from different field sizes"
        else:
            raise "Cannot add ffelt with element of type " + type(a)

        #return ffelt(self.f.poly_to_power[result], self.f, debug=self.debug)
        ff_result=copy.deepcopy(self)
        ff_result.elt=self.f.poly_to_power[result]
        return ff_result

    def __sub__(self,a):
        return self.__add__(a)

    def __mul__(self,a):
        if a==None or self.elt==None:
            result=None
        elif a==0:
            result = None
        elif a==1:
            result = self.elt
        elif type(a)==ffelt:
            if a.q==self.q:
                result = (self.elt+a.elt) % (self.q-1)
            else:
                raise "Cannot multiply elements from different field sizes"
        else:
            raise "Cannot multipy ffelt with element of type " + type(a)

        #return ffelt(result, self.f, debug=self.debug)
        ff_result=copy.deepcopy(self)
        ff_result.elt=result
        return ff_result

    def __truediv__(self,a):
        if a==None:
            raise "Division by zero not defined"
        elif self.elt==None:
            result=None
        elif type(a)==int:
            result = (self.elt-a) % (self.q-1)
        elif type(a)==ffelt:
            if a.q==self.q:
                result = (self.elt-a.elt) % (self.q-1)
            else:
                raise "Cannot multiply elements from different field sizes"
        else:
            raise "Cannot multipy ffelt with element of type " + type(a)

        #return ffelt(result, self.f, debug=self.debug)
        ff_result=copy.deepcopy(self)
        ff_result.elt=result
        return ff_result

    def __pow__ (self, a):

        if self.elt==None: 
            result=None
        elif type(a)==int:
            result = (self.elt*a) % (self.q-1)
        else:
            raise "power must be an integer"

        # Old way
        # return ffelt(result, self.f, debug=self.debug)
        ff_result=copy.deepcopy(self)
        ff_result.elt=result
        return ff_result
    
    def vec(self):
        
        return bin_array(self.f.power_to_poly[self.elt],self.f.m).tolist()

    def conjugates(self):
        result=[]
        conj=self.elt
        while True:
            result+=[ffelt(conj,self.f)]
            conj=(conj+conj)%(self.q-1)
            if conj==self.elt:
                break
        return result

    def __repr__ (self):
        #print("__repr__", self.elt)

        fieldString=""
        if not self.suppressField:
            fieldString=" GF("+str(self.q)+")"



        if self.elt==None:
            return ("0"+fieldString)
        elif self.elt==0:
            return ("1"+fieldString)
        elif self.elt>0:
            return ("a^"+str(self.elt)+fieldString)
        else:
            raise "Something went wrong in __repr__, self.elt="+str(self.elt)
        
    def __eq__ (self, a):
        if a==0:
            return self.elt==None
        elif a==1:
            return self.elt==0
        elif type(a)==ffelt:
            if a.q==self.q:
                return a.elt==self.elt
            else:
                raise "Cannot add elements from different field sizes"
                
    def minpoly(self, coeffs=False, debug=False):
        ''' 
        Find the minimum polynomial (in the base field) that has the element as a root

        If coeffs=False (default) returns a vector of polynomial coefficients

        If coeffs=True, returns the positions of the nonzero coefficients 
        (as in App. B of Lin and Costello)
        '''

        if self.elt > 0:
            num_conjugates=len(self.conjugates())
            if debug:
                print("Num conjugates=", num_conjugates)
            for i in range(2**(num_conjugates-1)):
                middle=bin_array(i,num_conjugates-1)
                candidate=np.hstack(([1],middle,[1]))
                if debug:
                    print(candidate)

                result=ffelt(None, self.f)
                for power, coeff in enumerate(candidate):
                    if coeff>0:
                        result= result + self**(len(candidate)-power-1)
                        if debug:
                            print(power, result)


                if debug:
                    print()
                    print(result)
                    print("---")
                if result==0:
                    candidate= candidate.tolist()
                    break
        elif self.elt==0:
            candidate=[1,1]
        elif self.elt==None:
            candidate=[1,0]
        else:
            raise "Something when wrong, self.elt="+str(self.elt)
            
        if  coeffs:
            return np.where(np.flip(candidate)>0)[0].tolist()
        else:
            return candidate
