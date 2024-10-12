# RSA-impelmentation
Implementation of  the RSA algorithm from scratch in C++

This project consists of implementing a C++ RSA program from scratch. First, we create a class called "Bignum" that handles large numbers. Then, we implement all the necessary operations for this class to use them in the RSA algorithm.
The program functions well, but the execution time becomes lengthy when the key size exceeds 256 bits.
Below are the execution times I observed while testing this program:
Key size = 128 bits: Execution time is less than 2 seconds.
Key size = 256 bits: Execution time ranges between 4 and 8 seconds.
Key size = 512 bits: Execution time ranges between 5 and 11 minutes.
Key size = 1024 bits: I haven't mustered the courage to test it!
To compile the program: $ make
To run it: $ make run

Thanks for reading

BADIANE BASSIROU
