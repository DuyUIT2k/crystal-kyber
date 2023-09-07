class Module:
    def __init__(self, ring):
        self.ring = ring # PolynomialRing
        
    # Decode byte sequences into module elements, which are matrices of elements
    def decode(self, input_bytes, m, n, l=None, is_ntt=False):
        if l is None:
            # Input length must be 32*l*m*n bytes long
            l, check = divmod(8*len(input_bytes), self.ring.n*m*n)
            if check != 0:
                raise ValueError("Input bytes must be a multiple of (polynomial degree) / 8")
        else:
            if self.ring.n*l*m*n > len(input_bytes)*8:
                raise ValueError("Byte length is too short for given l")
            
        chunk_length = 32*l
        byte_chunks = [input_bytes[i:i+chunk_length] for i in range(0, len(input_bytes), chunk_length)]
        matrix = [[0 for _ in range(n)] for _ in range(m)] # m rows and n cols
        for i in range(m):
            for j in range(n):
                mij = self.ring.decode(byte_chunks[n*i+j], l=l, is_ntt=is_ntt)
                matrix[i][j] = mij

        return self(matrix)

    # Define how instances of Module are represented as strings when using the repr() function
    def __repr__(self):
        return f"Module over the commutative ring: {self.ring}"

    # Define how instances of Module are represented as strings when using the str() function
    def __str__(self):
        return f"Module over the commutative ring: {self.ring}"

    def __call__(self, matrix_elements):
        if not isinstance(matrix_elements, list):
            raise TypeError(f"Elements of a module are matrices, with elements .")

        if isinstance(matrix_elements[0], list):
            for element_list in matrix_elements:
                if not all(isinstance(aij, self.ring.element) for aij in element_list):
                    raise TypeError(f"All elements of the matrix must be elements of the ring: {self.ring}")
                
            return Module.Matrix(self, matrix_elements)
        elif isinstance(matrix_elements[0], self.ring.element):
            if not all(isinstance(aij, self.ring.element) for aij in matrix_elements):
                raise TypeError(f"All elements of the matrix must be elements of the ring: {self.ring}")
            
            return Module.Matrix(self, [matrix_elements])
        else:
            raise TypeError(f"Elements of a module are matrices, built from elements of the base ring.")

    class Matrix:
        def __init__(self, parent, matrix_elements):
            self.parent = parent # An instance of Module class
            self.rows = matrix_elements # A 2D matrix of Polynomial class
            self.m = len(matrix_elements) # Number of rows
            self.n = len(matrix_elements[0]) # Number of cols
            if not self.check_dimensions():
                raise ValueError("Inconsistent row lengths in matrix")

        # Get dimension of the matrix (row, col)
        def get_dim(self):
            return self.m, self.n

        # Check all dimension of the matrix
        def check_dimensions(self):
            return all(len(row) == self.n for row in self.rows)

        # Transpose the self.rows matrix
        def transpose(self):
            new_rows = [list(item) for item in zip(*self.rows)]

            return self.parent(new_rows)

        # Transpose the self.rows matrix and assign again for self
        def transpose_self(self):
            self.m, self.n = self.n, self.m
            self.rows = [list(item) for item in zip(*self.rows)]

            return self
            
        # Reduce all coefficients value by mod q with each coefficient for polynomials in a row ~ polynomial ring
        def reduce_coefficents(self):
            for row in self.rows:
                for ele in row:
                    ele.reduce_coefficents()
            return self
            
        # Converts the coefficients of a polynomial "poly" to Montgomery form for polynomials in a row ~ polynomial ring
        def to_montgomery(self):
            for row in self.rows:
                for ele in row:
                    ele.to_montgomery()

            return self

        # Encodes the polynomial as a byte array for polynomials in a row ~ polynomial ring
        def encode(self, l=None):
            output = b""
            for row in self.rows:
                for j in range(self.n): # self.n: number of columns in a matrix
                    output += row[j].encode(l=l)

            return output
            
        # Compresses the polynomial coefficients using lossy compression for polynomials in a row ~ polynomial ring
        def compress(self, d):
            for row in self.rows:
                for ele in row:
                    ele.compress(d)

            return self
        
        # Decompresses the polynomial coefficients for polynomials in a row ~ polynomial ring
        def decompress(self, d):
            for row in self.rows:
                for ele in row:
                    ele.decompress(d)

            return self    
    
        # Converts a polynomial to Number Theoretic Transform (NTT) form for polynomials in a row ~ polynomial ring
        def to_ntt(self):
            for row in self.rows:
                for ele in row:
                    ele.to_ntt()

            return self
    
        # Converts a polynomial from NTT form and performs a multiplication by a Montgomery factor for polynomials in a row ~ polynomial ring
        def from_ntt(self):
            for row in self.rows:
                for ele in row:
                    ele.from_ntt()

            return self        
                    
        # Get a specific row in a matrix through index "i"
        def __getitem__(self, i):
            return self.rows[i]

        # Compare 2 matrices to see if they are equal or not
        def __eq__(self, other):
            return other.rows == self.rows

        # Returns the sum of the self and other matrices
        # Behave like self_matrix + other_matrix
        def __add__(self, other):
            if not isinstance(other, Module.Matrix):
                raise TypeError("Can only add matrcies to other matrices")
            
            if self.parent != other.parent:
                raise TypeError("Matricies must have the same base ring")
            
            if self.get_dim() != other.get_dim():
                raise ValueError("Matrices are not of the same dimensions")

            new_elements = []
            for i in range(self.m):
                new_elements.append([a+b for a,b in zip(self.rows[i], other.rows[i])])

            return self.parent(new_elements)

        # Returns the sum of the self matrix (in the right side of +) and another operand (in the left side of +)
        # Behave like other_operand + self_matrix
        def __radd__(self, other):
            return self.__add__(other)

        # Returns the sum of the self and other matrices 
        # Behave like self_matrix += other_matrix
        def __iadd__(self, other):
            self = self + other
            
            return self

        # Returns the difference of the self and other matrices 
        # Behave like self_matrix - other_matrix
        def __sub__(self, other):
            if not isinstance(other, Module.Matrix):
                raise TypeError("Can only subtract matrcies from other matrices")
            
            if self.parent != other.parent:
                raise TypeError("Matricies must have the same base ring")
            
            if self.get_dim() != other.get_dim():
                raise ValueError("Matrices are not of the same dimensions")

            new_elements = []
            for i in range(self.m):
                new_elements.append([a-b for a,b in zip(self.rows[i], other.rows[i])])

            return self.parent(new_elements)

        # Returns the difference of the self matrix (in the right side of -) and another operand (in the left side of -) 
        # Behave like other_operand - self_matrix
        def __rsub__(self, other):
            return self.__sub__(other)

        # Returns the difference of the self and other matrices. Behave like self_matrix -= other_matrix
        def __isub__(self, other):
            self = self - other

            return self

        # Returns the product of the self and other matrices. Behave like self_matrix @ other_matrix
        def __matmul__(self, other):
            """
            Denoted A @ B
            """
            if not isinstance(other, Module.Matrix):
                raise TypeError("Can only multiply matrcies with other matrices")
            
            if self.parent != other.parent:
                raise TypeError("Matricies must have the same base ring")
            
            if self.n != other.m:
                raise ValueError("Matrices are of incompatible dimensions")

            new_elements = [[sum(a*b for a,b in zip(A_row, B_col)) for B_col in other.transpose().rows] for A_row in self.rows]

            return self.parent(new_elements)

        # Define how instances of Matrix are represented as strings when using the repr() function
        def __repr__(self):
            if len(self.rows) == 1:
                return str(self.rows[0])
            
            max_col_width = []
            for n_col in range(self.n):
                max_col_width.append(max(len(str(row[n_col])) for row in self.rows))
                
            info = ']\n['.join([', '.join([f'{str(x):>{max_col_width[i]}}' for i,x in enumerate(r)]) for r in self.rows])

            return f"[{info}]"

