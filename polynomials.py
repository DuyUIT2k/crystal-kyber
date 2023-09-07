import random
from utils import *

class PolynomialRing:
    """
    Initialise the polynomial ring:
        
        R = GF(q) / (X^n + 1) 
    """
    def __init__(self, q, n, ntt_helper=None):
        self.q = q # The finite field size
        self.n = n # The degree of the polynomial
        self.element = PolynomialRing.Polynomial # A polynomial class instance
        self.ntt_helper = ntt_helper # An optional helper class for performing Number Theoretic Transforms (NTT)

    # Generates a polynomial with coefficients [0, 1]. 
    # The "is_ntt" parameter specifies whether the polynomial should be in NTT form.
    def gen(self, is_ntt=False):
        return self([0,1], is_ntt=is_ntt) # use __call__ method -> Polynomial

    # Generates a random polynomial with coefficients in the range [0, q-1]. 
    # The "is_ntt" parameter specifies whether the polynomial should be in NTT form.
    def random_element(self, is_ntt=False):
        coefficients = [random.randint(0, self.q - 1) for _ in range(self.n)]

        return self(coefficients, is_ntt=is_ntt) # use __call__ method -> Polynomial
        
    # Convert a byte stream to the NTT representation with q = 3329 (Parses a byte array into a polynomial)
    # It extracts coefficients from the byte array and creates a polynomial.
    def parse(self, input_bytes, is_ntt=False):
        """
        Algorithm 1 (Parse)
        https://pq-crystals.org/kyber/data/kyber-specification-round3-20210804.pdf
        
        Parse: B^* -> R
        """
        i, j = 0, 0
        coefficients = [0 for _ in range(self.n)]
        while j < self.n:
            d1 = input_bytes[i] + 256*(input_bytes[i+1] % 16)
            d2 = (input_bytes[i+1] // 16) + 16*input_bytes[i+2]
            
            if d1 < self.q:
                coefficients[j] = d1
                j = j + 1
            
            if d2 < self.q and j < self.n:
                coefficients[j] = d2
                j = j + 1
                
            i = i + 3

        return self(coefficients, is_ntt=is_ntt)      

    # Performs Centered Binomial Distribution (CBD) on a byte array and converts it into a polynomial. 
    # This is used in the Kyber algorithm to generate random polynomials.
    def cbd(self, input_bytes, eta, is_ntt=False):
        """
        Algorithm 2 (Centered Binomial Distribution)
        https://pq-crystals.org/kyber/data/kyber-specification-round3-20210804.pdf
        
        Expects a byte array of length (eta * polynomial_degree / 4) ~ eta * (polynomial_degree >> 2)
        For Kyber, this is 64 eta.
        """
        assert (self.n >> 2)*eta == len(input_bytes) # Ensure that input are 64 bytes long
        coefficients = [0 for _ in range(self.n)]
        list_of_bits = bytes_to_bits(input_bytes) # Convert 64 bytes to 512 bits
        for i in range(self.n):
            a = sum(list_of_bits[2*i*eta + j]       for j in range(eta))
            b = sum(list_of_bits[2*i*eta + eta + j] for j in range(eta))
            coefficients[i] = a-b

        return self(coefficients, is_ntt=is_ntt)
        
    # Decodes a byte array into a polynomial. 
    # This is used for decoding in Kyber.
    def decode(self, input_bytes, l=None, is_ntt=False):
        """
        Decode (Algorithm 3)
        
        decode: B^32l -> R_q
        """
        if l is None:
            l, check = divmod(8*len(input_bytes), self.n)
            if check != 0:
                raise ValueError("Input bytes must be a multiple of (polynomial degree) / 8")
        else:
            if self.n*l != len(input_bytes)*8:
                raise ValueError("Input bytes must be a multiple of (polynomial degree) / 8")
        coefficients = [0 for _ in range(self.n)]
        list_of_bits = bytes_to_bits(input_bytes)
        for i in range(self.n):
            coefficients[i] = sum(list_of_bits[i*l + j] << j for j in range(l))

        return self(coefficients, is_ntt=is_ntt)
            
    def __call__(self, coefficients, is_ntt=False):
        if isinstance(coefficients, int):
            return self.element(self, [coefficients], is_ntt)
        if not isinstance(coefficients, list):
            raise TypeError(f"Polynomials should be constructed from a list of integers, of length at most d = {self.n}")
        
        return self.element(self, coefficients, is_ntt) # Call the __init__ method of Polynomial

    # Define how instances of PolynomialRing are represented as strings when using the repr() function
    def __repr__(self):
        return f"Univariate Polynomial Ring in x over Finite Field of size {self.q} with modulus x^{self.n} + 1"

    class Polynomial:
        def __init__(self, parent, coefficients, is_ntt=False):
            self.parent = parent # An instance of PolynomialRing class
            self.coeffs = self.parse_coefficients(coefficients) # Coefficients in a polynomial
            self.is_ntt = is_ntt # Is polynomial in NTT form or not

        # Check whether that all coefficients are zero or not
        def is_zero(self):
            """
            Return if polynomial is zero: f = 0
            """
            return all(c == 0 for c in self.coeffs)

        # Check whether the polynomial is a 0 degree (There is just a constant) or not
        def is_constant(self):
            """
            Return if polynomial is constant: f = c
            """
            return all(c == 0 for c in self.coeffs[1:])
            
        # Add more "0" coefficients to the polynomial ~ Add "0" padding
        def parse_coefficients(self, coefficients):
            """
            Helper function which right pads with zeros
            to allow polynomial construction as 
            f = R([1,1,1])
            """
            l = len(coefficients)
            if l > self.parent.n:
                raise ValueError(f"Coefficients describe polynomial of degree greater than maximum degree {self.parent.n}")
            elif l < self.parent.n:
                coefficients = coefficients + [0 for _ in range (self.parent.n - l)]

            return coefficients
            
        # Reduce all coefficients value by mod q with each coefficient
        def reduce_coefficents(self):
            """
            Reduce all coefficents modulo q
            """
            self.coeffs = [c % self.parent.q for c in self.coeffs]

            return self
 
        # Encodes the polynomial as a byte array
        def encode(self, l=None):
            """
            Encode (Inverse of Algorithm 3)
            """
            if l is None:
                l = max(x.bit_length() for x in self.coeffs)

            bit_string = ''.join(format(c, f'0{l}b')[::-1] for c in self.coeffs)

            return bitstring_to_bytes(bit_string)

        # Compresses the polynomial coefficients using lossy compression
        def compress(self, d):
            """
            Compress the polynomial by compressing each coefficent
            NOTE: This is lossy compression
            """
            compress_mod   = 2**d
            compress_float = compress_mod / self.parent.q
            self.coeffs = [round_up(compress_float * c) % compress_mod for c in self.coeffs]

            return self
            
        # Decompresses the polynomial coefficients
        def decompress(self, d):
            """
            Decompress the polynomial by decompressing each coefficent
            NOTE: This as compression is lossy, we have
            x' = decompress(compress(x)), which x' != x, but is 
            close in magnitude.
            """
            decompress_float = self.parent.q / 2**d
            self.coeffs = [round_up(decompress_float * c) for c in self.coeffs]

            return self

        # Add two coefficents and then modulo q
        def add_mod_q(self, x, y):
            """
            add two coefficents modulo q
            """
            tmp = x + y
            if tmp >= self.parent.q:
                tmp -= self.parent.q

            return tmp

        # Subtract two coefficents and then modulo q
        def sub_mod_q(self, x, y):
            """
            sub two coefficents modulo q
            """
            tmp = x - y
            if tmp < 0:
                tmp += self.parent.q

            return tmp
            
        # An implementation of polynomial multiplication using the schoolbook multiplication algorithm. 
        # This algorithm is suitable for polynomial rings over finite fields, particularly for rings of the form R_q = F_1[X]/(X^n + 1), where F_1 represents a finite field.
        def schoolbook_multiplication(self, other):
            """
            Naive implementation of polynomial multiplication
            suitable for all R_q = F_1[X]/(X^n + 1)
            """
            n = self.parent.n
            a = self.coeffs
            b = other.coeffs
            new_coeffs = [0 for _ in range(n)]

            for i in range(n):
                for j in range(0, n-i):
                    new_coeffs[i+j] += (a[i] * b[j])

            for j in range(1, n):
                for i in range(n-j, n):
                    new_coeffs[i+j-n] -= (a[i] * b[j])

            return [c % self.parent.q for c in new_coeffs]
        
        """
        The next four `Polynomial` methods rely on the parent
        `PolynomialRing` having a `ntt_helper` from 
        ntt_helper.py and are used for NTT speediness.
        """
        # Converts a polynomial to Number Theoretic Transform (NTT) form
        def to_ntt(self):
            if self.parent.ntt_helper is None:
                raise ValueError("Can only perform NTT transform when parent element has an NTT Helper")
            
            return self.parent.ntt_helper.to_ntt(self)
        
        # Converts a polynomial from NTT form and performs a multiplication by a Montgomery factor
        def from_ntt(self):
            if self.parent.ntt_helper is None:
                raise ValueError("Can only perform NTT transform when parent element has an NTT Helper")
            
            return self.parent.ntt_helper.from_ntt(self)
            
        # Converts the coefficients of a polynomial "poly" to Montgomery form
        def to_montgomery(self):
            """
            Multiply every element by 2^16 mod q
            
            Only implemented (currently) for n = 256
            """
            if self.parent.ntt_helper is None:
                raise ValueError("Can only perform Mont. reduction when parent element has an NTT Helper")
            
            return self.parent.ntt_helper.to_montgomery(self)
        
        # Multiplies two sets of polynomial coefficients "f_coeffs" and "g_coeffs" using NTT-based multiplication
        def ntt_multiplication(self, other):
            """
            Number Theoretic Transform multiplication.
            Only implemented (currently) for n = 256
            """
            if self.parent.ntt_helper is None:
                raise ValueError("Can only perform ntt reduction when parent element has an NTT Helper")
            
            if not (self.is_ntt and other.is_ntt):
                raise ValueError("Can only multiply using NTT if both polynomials are in NTT form")
            
            # function in ntt_helper.py
            new_coeffs = self.parent.ntt_helper.ntt_coefficient_multiplication(self.coeffs, other.coeffs)

            return self.parent(new_coeffs, is_ntt=True)

        # Returns the negation of the polynomial. Behave like 0 - polynomial ~ -polynomial
        def __neg__(self):
            """
            Returns -f, by negating all coefficients
            """
            neg_coeffs = [(-x % self.parent.q) for x in self.coeffs]
            
            return self.parent(neg_coeffs, is_ntt=self.is_ntt)

        # Returns the sum of the self and other polynomials. Behave like self_polynomial + other_polynomial
        def __add__(self, other):
            if isinstance(other, PolynomialRing.Polynomial):
                if self.is_ntt ^ other.is_ntt:                    
                    raise ValueError(f"Both or neither polynomials must be in NTT form before multiplication")
                new_coeffs = [self.add_mod_q(x,y) for x,y in zip(self.coeffs, other.coeffs)]
            elif isinstance(other, int):
                new_coeffs = self.coeffs.copy()
                new_coeffs[0] = self.add_mod_q(new_coeffs[0], other)
            else:
                raise NotImplementedError(f"Polynomials can only be added to each other")
            
            return self.parent(new_coeffs, is_ntt=self.is_ntt)

        # Returns the sum of the self polynomial (in the right side of +) and another operand (in the left side of +). Behave like other_operand + self_polynomial
        def __radd__(self, other):
            return self.__add__(other)

        # Returns the sum of the self and other polynomials. Behave like self_polynomial += other_polynomial
        def __iadd__(self, other):
            self = self + other

            return self

        # Returns the difference of the self and other polynomials. Behave like self_polynomial - other_polynomial
        def __sub__(self, other):
            if isinstance(other, PolynomialRing.Polynomial):
                if self.is_ntt ^ other.is_ntt:
                    raise ValueError(f"Both or neither polynomials must be in NTT form before multiplication")
                new_coeffs = [self.sub_mod_q(x,y) for x,y in zip(self.coeffs, other.coeffs)]
            elif isinstance(other, int):
                new_coeffs = self.coeffs.copy()
                new_coeffs[0] = self.sub_mod_q(new_coeffs[0], other)
            else:
                raise NotImplementedError(f"Polynomials can only be subracted from each other")
            
            return self.parent(new_coeffs, is_ntt=self.is_ntt)

        # Returns the subtraction of the self polynomial (in the right side of -) and another operand (in the left side of -). Behave like other_operand - self_polynomial
        def __rsub__(self, other):
            return self.__sub__(other)

        # Returns the difference of the self and other polynomials. Behave like self_polynomial -= other_polynomial
        def __isub__(self, other):
            self = self - other

            return self

        # Returns the product of the self and other polynomials. Behave like self_polynomial * other_polynomial
        def __mul__(self, other):
            if isinstance(other, PolynomialRing.Polynomial):
                if self.is_ntt and other.is_ntt:
                    return self.ntt_multiplication(other)
                elif self.is_ntt ^ other.is_ntt:
                     raise ValueError(f"Both or neither polynomials must be in NTT form before multiplication")
                else:
                    new_coeffs = self.schoolbook_multiplication(other)
            elif isinstance(other, int):
                new_coeffs = [(c * other) % self.parent.q for c in self.coeffs]
            else:
                raise NotImplementedError(f"Polynomials can only be multiplied by each other, or scaled by integers")
            
            return self.parent(new_coeffs, is_ntt=self.is_ntt)

        # Returns the product of the self polynomial (in the right side of *) and another operand (in the left side of *). Behave like other_operand * self_polynomial
        def __rmul__(self, other):
            return self.__mul__(other)

        # Returns the product of the self and other polynomials. Behave like self_polynomial *= other_polynomial
        def __imul__(self, other):
            self = self * other
            return self

        # Define exponentiation behavior for a polynomial instance
        def __pow__(self, n):
            if not isinstance(n, int):
                raise TypeError(f"Exponentiation of a polynomial must be done using an integer.")

            # Deal with negative scalar multiplication
            if n < 0:
                raise ValueError(f"Negative powers are not supported for elements of a Polynomial Ring")
            
            f = self
            g = self.parent(1, is_ntt=self.is_ntt)
            while n > 0:
                if n % 2 == 1:
                    g = g * f
                f = f * f
                n = n // 2

            return g

        # Define the equality (==) comparison between polynomial instances
        def __eq__(self, other):
            if isinstance(other, PolynomialRing.Polynomial):
                return self.coeffs == other.coeffs and self.is_ntt == other.is_ntt
            elif isinstance(other, int):
                if self.is_constant() and (other % self.parent.q) == self.coeffs[0]:
                    return True
                
            return False

        # Allow you to access the coefficients of the polynomial using indexing (e.g., polynomial[0], polynomial[1], ...)
        def __getitem__(self, idx):
            return self.coeffs[idx]

        # Define how instances of Polynomial are represented as strings when using the repr() function
        def __repr__(self):
            ntt_info = ""

            if self.is_ntt:
                ntt_info = " (NTT form)"

            if self.is_zero():
                return "0" + ntt_info

            info = []
            for i,c in enumerate(self.coeffs):
                if c != 0:
                    if i == 0:
                        info.append(f"{c}")
                    elif i == 1:
                        if c == 1:
                            info.append("x")
                        else:
                            info.append(f"{c}*x")
                    else:
                        if c == 1:
                            info.append(f"x^{i}")
                        else:
                            info.append(f"{c}*x^{i}")
                            
            return " + ".join(info) + ntt_info

        # Define how instances of Polynomial are represented as strings when using the str() function
        def __str__(self):
            return self.__repr__()