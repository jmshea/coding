# coding
Error-control coding related code

## binext.py

**a Python module for operating on binary extension fields, GF(2^m)**

Provides classes for creating and operating in finite fields:

*class* ff: general finite field, generates and stores tables for translating between power and vector representations in GF(2^m)
  
  Also provides method for printing out table of minimum polynomials, in same format as Lin and Costello book

*class* ffelt: general finite field element, can be used to carry out arithmetic on the field

  Also provides methods to find conjugates, minimum polynomials, and to output vector representations
  
