#BBFMM3D  



###1. INTRODUCTION

BBFMM3D is an open source package of the <a href="http://www.sciencedirect.com/science/article/pii/S0021999109004665">Black-box Fast Multipole Method</a> in 3 dimensions.   
The Black-box Fast Multipole Method is an O(N) fast multipole method, which is a technique to calculate sums of the form  

![](http://latex.codecogs.com/gif.latex?f%28x_i%29%20%3D%20%5Cdisplaystyle%20%5Csum_%7Bj%3D1%7D%5EN%20K%28x_i%2Cy_j%29%20%5Csigma_j%2C%20%5C%2C%5C%2C%5C%2C%20%5Cforall%20i%20%5Cin%5C%7B1%2C2%2C%5Cldots%2CN%5C%7D)

where ![](http://latex.codecogs.com/gif.latex?K%28x_i%2Cx_j%29) is kernel function, ![](http://latex.codecogs.com/gif.latex?x_i) are observation points, ![](http://latex.codecogs.com/gif.latex?y_j) are locations of sources, and ![](http://latex.codecogs.com/gif.latex?%5Csigma_i) are charges at corresponding locations.
BBFMM3D provides an O(N) solution to matrix-vector products of the type Ax. In that case the relation between A and K is:
![](http://latex.codecogs.com/gif.latex?A_%7Bij%7D%20%3D%20K%28x_i%2Cy_j%29)

This implementation of the FMM differs from other methods by the fact that it is applicable to all smooth kernels K. [Give examples of RBF kernels, 1/r, log r, Stokes, etc.].

The approximation scheme used in the FMM relies on Chebyshev interplation to construct low-rank approximations for well-separated clusters. In addition the use of Singular Value Decomposition ensures that the computational cost is minimal. In particular the rank is optimally chosen for a given error. 

Please cite the following paper if you use this code:

Fong, William, and Eric Darve. "The black-box fast multipole methodshod." Journal of Computational Physics 228, no. 23 (2009): 8712-8725. You can see details <a href="http://www.sciencedirect.com/science/article/pii/S0021999109004665">here</a>.

###2. DIRECTORIES AND FILES


	./examples/		:	Example input C++ codes; Needed to read input from user or from input file.  
	./src/			:	Source code in C++  
	./include/		:	Relevant header files  
	./exec/			:	Executables for BBFMM3D  
	./input/		:	The input file.  
	./README.md		:	This file  
	./License.md	:	License file  
	./Makefile		:	Makefile
	
###3. TUTORIAL
####3.1 To Get Started  
To check whether things are set up correctly, you can perform the following: Go to the directory where Makefile is in, then key in the following three commands in the terminal:

		make binary_file_mykernel
		cd exec/
		./binary_file_mykernel

####3.2 Basic usage

#####3.2.1 BBFMM3D with standard kernel

The basic usage of BBFMM3D with standard kernel is as follows: 

	#include"bbfmm3d.hpp"  
	...
	{
	double L;       // Length of simulation cell (assumed to be a cube)
    int n;          // Number of Chebyshev nodes per dimension
    doft dof;
    int Ns;         // Number of sources in simulation cell
    int Nf;         // Number of field points in simulation cell
    int m;
    int level;		// The number of levels in the hierarchy tree
    int use_chebyshev // label of whether the computation will use chebyshev interpolation formula or uniform interpolation formula 
    double eps;
    vector3 *source = new vector3[Ns];    // Position array for the source points
    vector3 *field = new vector3[Nf];     // Position array for the field points
    double *q =  new double[Ns*dof.s*m];  // Source array
    double *stress      =  new double[Nf*dof.f*m];// Field array (BBFMM calculation)
    …
	kernel_LaplacianForce Atree(L,level, n, eps, use_chebyshev);
    Atree.buildFMMTree();  // Build the fmm tree;
	
	/* The following can be repeated with different field, source, and q */
    
    H2_3D_Compute<kernel_LaplacianForce> compute(&Atree, field, source, Ns, Nf, q,m, stress);
    ...
    }
    
This example first build a FMM tree with these two lines:  

	kernel_LaplacianForce Atree(L,level, n,  eps, use_chebyshev);
    Atree.buildFMMTree();  
where kernel_LaplacianForce is a class of fmm tree using LaplacianForce kernel, the constructor takes 5 arguments:  

* L(double):   
	Length of simulation cell (assumed to be a cube).
* level(int):  
	The number of levels in the hierarchy tree
* n(int):  
	Number of Chebyshev nodes per dimension. 
	 
* eps(double):  
	Target accuracy; this is used to determine which
 singular values are kept after the SVD. cutoff.s and cutoff.f are computed
 using this epsilon.  
* use_chebyshev(int):  
	Label to indicate whether using chebyshev interpolation formula or uniform interpolation formula.  
	use_chebyshev = 1: chebyshev interplation formula;  
	use_chebyshev = 0: uniform interpolation formula(where FFT is used, now it just supports homogenous kernel, non-homogeneous will be added)
 
Once the tree is created, you can repeat matrix-vector product with different field, source and q(charge).(see **3.2.4**) The code shows an example using LapacianForce kernel:  

	H2_3D_Compute<kernel_LaplacianForce> compute(&Atree, field, source, Ns, Nf, q,m, stress);
The template class `H2_3D_Compute` is for computing matrix-vector product after the FMM tree is built. The constructor takes 8 arguments:  

* &Atree(T*):  
	A pointer to the FMM tree we just built. Here T is typename, FMM tree built with different kernels have different typenames.
* field(vector3*):  
	A pointer to position array of field points. vector3 is a struct storing x, y and z coordinates of a 3D point.  
* source(vector3*):  
	A pointer to position array of source points.  
* Ns(int):   
	Number of sources in simulation cell.
* Nf(int):  
	Number of fields in simulation cell  
* q(double*):  
	A pointer to the charges.
* m(int):   
	Number of sets of charges.  
* stress(double*):   
	A pointer to the result, and the result is stored column-wise in `stress`.

#####3.2.2 Options of provided kernels

Below are the details of the kernel functions K we have provided:  
( For all the kernel functions, we denote r to be Euclidean distance between x and y. )

Options of kernels:  

* LOGARITHM kernel:           
	usage: kernel_Logarithm  
	kernel function:  
    ![](http://latex.codecogs.com/gif.latex?K%28x%2Cy%29%20%3D%200.5%20%5Clog%28r%5E2%29%5C%2C%20%28r%5Cneq%200%29%3B%5C%2C%20K%28x%2Cy%29%3D%200%20%5C%2C%28r%3D0%29.)
	
	
* ONEOVERR2 kernel:  
	usage: kernel_OneOverR2  
	kernel function:  
    ![](http://latex.codecogs.com/gif.latex?K%28x%2Cy%29%20%3D%201%20/%20r%5E2%20%5C%2C%28r%20%5Cneq%200%29%3B%5C%2C%20K%28x%2Cy%29%3D%200%20%5C%2C%28r%3D0%29%24) 
	
* GAUSSIAN kernel:  
	usage: kernel_Gaussian  
	kernel function:  
	![](http://latex.codecogs.com/gif.latex?K%28x%2Cy%29%20%3D%20exp%28-r%5E2%29)   
	
* QUADRIC kernel:  
	usage: kernel_Quadric  
	kernel function:  
	![](http://latex.codecogs.com/gif.latex?K%28x%2Cy%29%20%3D%201%20&plus;%20r%5E2)
* INVERSEQUADRIC kernel:  
	usage: kernel_InverseQuadric  
	kernel function:  
	![](http://latex.codecogs.com/gif.latex?K%28x%2Cy%29%20%3D%201%20/%20%281&plus;r%5E2%29)	
* THINPLATESPLINE kernel:  
	usage:  kernel_ThinPlateSpline  
	kernel function:  
	![](http://latex.codecogs.com/gif.latex?K%28x%2Cy%29%20%3D%200.5%20r%5E2%20%5Clog%28r%5E2%20%29%5C%2C%20%28r%20%5Cneq%200%29%3B%5C%2C%20K%28x%2Cy%29%3D0%5C%2C%28r%3D0%29) 

* LAPLACIAN kernel:  
	usage: kernel_Lapacian  
	kernel function:  
	![](http://latex.codecogs.com/gif.latex?K%28x%2Cy%29%20%3D%201%20/%20r)
	
* ONEOVERR4 kernel:  
	usage: kernel_OneOverR4  
	kernel function:  
	![](http://latex.codecogs.com/gif.latex?K%28x%2Cy%29%20%3D%201%20/%20r%5E4)	
* LAPLACIANFORCE kernel:  
	usage: kernel_LaplacianForce  
* STOKES kernel (tensor kernel)  
	usage: kernel_Stokes  
	kernel function:  
	<img src="http://www.sciweavers.org/upload/Tex2Img_1441257351/render.png" border="0"/>
	i.e. <img src="http://www.sciweavers.org/upload/Tex2Img_1441257386/render.png" border="0"/> 
	where delta is an indicator function. vector r is the difference of x and y, see http://en.wikipedia.org/wiki/Stokes_flow#By_Green.27s_function:_the_Stokeslet for details.
	
	
    		
If you want to define your own kernel, please see **3.2.3**.

#####3.2.3 BBFMM3D with user defined kernel

The basic usage is almost the same as **3.2.1** except that you have to define your own routine of computing kernel. One example code is as follows:  

	#include "bbfmm3d.hpp"
	
	class myKernel: public H2_3D_Tree {
	public:
    myKernel(double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(L,level,n, epsilon, use_chebyshev){};
    virtual void setHomogen(string& kernelType,doft* dof) {
        homogen = -1;
        symmetry = 1;
        kernelType = "myKernel";
        dof->f = 1;
        dof->s = 1;
    }
    virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof) {
        vector3 diff;
        double rinv;
        
        diff.x = sourcepos.x - fieldpos.x;
        diff.y = sourcepos.y - fieldpos.y;
        diff.z = sourcepos.z - fieldpos.z;
        *K = 1. / sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
    }
	};

	

	...
	{
	…
	myKernel Atree(L,level, n, eps, use_chebyshev);
    Atree.buildFMMTree();  // Build the fmm tree;
	
	/* The following can be repeated with different field, source, and q */
    
    H2_3D_Compute<myKernel> compute(&Atree, field, source, Ns, Nf, q, m, stress);
    ...
    }
    
    
* dof(doft*):   
	A pointer to degree of freedom, which is the sizes of the small tensor matrix(tensor kernel). The FMM can also handle case of tensor kernel(not only scalar kernel). This struct stores information about the size of tensor kernel. 
	
You can define your own kernel inside `EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof)`, it takes field point, source point and degree of freedom(see **3.2.1**) as input, and pass the kernel value to K. 
                                                                
You also need to define information about kernel inside `setHomogen(string& kernelType,doft* dof)` 

* homogen:
	The homogeneous property of kernel.(The cost and memory requirements of the pre-computation step can be reduced for the case of homogeneous kernels.)
	* homogen = 0: if the kernel funtion is not homogeneous.
	* homogen = m: if the kernel funtion is homogeneous of degree m, i.e. ![](http://latex.codecogs.com/gif.latex?K%28%5Calpha%20x%2C%20%5Calpha%20y%29%20%3D%20%5Calpha%5Em%20K%28x%2Cy%29)
	* symmetry:   
	The symmetry property of the kernel.
	* symmetry = 0:  no symmetric property; 
	* symmetry = 1: symmetric kernel; 		[K(x,y) = K(y,x)]
	* symmetry = -1: anti-symmetric kernel.[K(x,y) = - K(y,x)]

*  kernelType:   
	The name of the kernel defined by user. The purpose of this is to label different files generated in the pre-computation step ( see **3.3** ) 
                      
                                
                           
#####3.2.4 Usage of multiple sources with same kernel

If you want to compute with different sources (source points, field points, charges) but with same kernel and number of Chebyshev nodes, you can do it in one file:
e.g.  

	class myKernel: public H2_3D_Tree {
	...
	}
	…
	
	{
	…
	/* Build FMM tree */
	myKernel Atree(L,level, n, eps, use_chebyshev);
    Atree.buildFMMTree(); 
	
	/* The following can be repeated with different field, source, and q */
    
    H2_3D_Compute<myKernel> compute1(&Atree, field1, source1, Ns1, Nf1, q1,m1, stress1);
	H2_3D_Compute<myKernel> compute2(&Atree, field2, source2, Ns2, Nf2, q2,m2, stress2);
    …
    }
    

The basic usage is already domonstrated in **3.2.1** and **3.2.3**. Once you have built the FMM tree, you can use different sources to compute the matrix-vector multiplication without rebuilding the tree. You can choose kernel type from standard kernels given by us ( see **3.2.2** ), or you can define your own kernel ( see **3.2.3** )

####3.3 Pre-computation

The power of this package is in the pre-computing part, which is much more computationally expensive than computing part. This package takes advantage of the fact that for a given kernel and number of chebyshev nodes, the precomputing part is the same, so for a fixed kernel and number of chebyshev nodes, it generates 3 files storing information of FMM tree in the folder /output. Everytime when we use the same kernel type and number of chebyshev nodes, we can directly read from the files, which would save a lot of time.

Note: it is true that sometimes with the pre-computation step, the code will be slower than direct calculation. But if the file already exists, then when doing more computations it will be faster than direct calculation. If you are using your own kernel, make sure to either change the kernelType or delete the existed file if you changed your kernel. 


####3.4 Test Interplation Error  
To give the user an idea of how large the interplation error is, we have provided a routine of computing the interplation error. If you want to test the interplation error between a cluster A (of size length) and a cluster B, where B is in the interplation list of A, you can do the following:  

	kernel_LaplacianForce testTree(&dof,1/pow(2,2),2, n, eps, use_chebyshev);
    double errtest = testInterplationErr(&testTree, 100, 100);
    
The first line constructs a FMM tree, you can see **3.2.1** for more details of the constructor. Here we more explanations of the second and third arguments:  

* second argument:  
  if you want to test the error with cluster of size `length`, then the second argument should be set to `4 x length`.
* third argument:  
  this argument should be set to 2.
  
The second line is the routine to compute interplation error.  

* &testTree (T* ):  
	A pointer to the tree you just created. 
* 100 (int):   
	Number of source points in the cluster.
* 100 (int):  
	Number of field points in the cluster.
	




###4. ROUTINES FOR INPUTING AND OUTPUTING DATA 
We have provided several routines for reading data from binary file, and writing data into binary file.	

####4.1 Reading meta data from text file
	
	void read_Metadata(const string& filenameMetadata,double& L, int& n, doft& dof, int& Ns, int& Nf, int& m, int& level);
The first argument, filenameMetadata is the filename for your meta data. L stores the length of simulation cell (assumed to be a cube); n stores the number of chebyshev nodes per dimension; dof stores (???); Ns stores the number of source points; Nf stores the number of field points; m stores the number of sets of charges; level stores the number of levels in the hierarchy tree;

**File format:**  
 
`L, n, dof.s, dof.f, Ns, Nf, m, level`

For example:

 	1,4,9,6,800,500,1,2
 	 	
####4.2 Reading from binary file  

	void read_Sources(const string& filenameField, vector3 *field, const int& Nf, const string& filenameSource, vector3 *source, const int& Ns, const string& filenameCharge, double *q, const int& m, const doft& dof);

The first argument filenameField, the forth argument filenameSource and the seventh argument filenameCharge are binary file names for field positions, source positions and charges respectively. Nf, Ns, m and dof are passed to this funtion. Nf is the number of field points, Ns is the number of source points, m is the number of sets of charges, and dof is the degree of freedom. The data of field is stored in `field`, the data of source is stored in `source`, and the data of charges is stored in `q` column-wise.  

**File format:** 
 
1. Binary file for source: 

	source positons are stored column-wise:   
	
		sourcepos0.x
		sourcepos1.x
		…
		sourcepos0.y
		sourcepos1.y
		…
		sourcepos0.z
		sourcepos1.z
		...

2. Binary file for field:

	field positions are stored column-wise:
		
		fieldpos0.x
		fieldpos1.x
		…
		fieldpos0.y
		fieldpos1.y
		…
		fieldpos0.z
		fieldpos1.z
		...
		
3. Binary file for charges:  
	charges are stored column-wise:  
	
		first set of charges
		second set of charges
		...
	
		
	
####4.3 Writing into binary file  
	
	void write_Into_Binary_File(const string& filename, double* outdata, int numOfElems);  
This first argument is the filename for your output data. The second argument is a pointer to the output data, and the last argument is the number of elements in the array of your output data.  


###5. EXAMPLES

We have provided several examples for BBFMM3D. Go to examples/, read through the files both must be self explanatory for the most part.
You can use our examples with your own input.
####5.1 Example file for generating binary/text files.
Both gen_binary_file.m and test.cpp shows how our test binary files are generated.
####5.2 Making changes to the examples for your own application

1. If you want to generate input through your own routine, and use the standard kernels:

    Go to `/examples/`, open `"get_input_through_routine_standard_kernel.cpp"`.        
    * To generate input through routine:   
    
        Change `SetMetaData()` and `SetSources()`  
    * To use standard kernels:   
      
      Choose the kernel type in `main()`, options of kernels are in **3.2.2**
  


3. If you want to generate input through your own routine, and use your own kernel:

    Go to `/examples/`, open `"get_input_through_routine_myKernel.cpp"`.    
    * To define your own kernel: 
     
      Modify `class myKernel`. 
    * To generate your input:  
    
      The same step as described in 1.

  	     
5. If you want to read input from binary file, and use standard kernel:
	
	Go to `/examples/`, open `"binary_file_standard_kernel.cpp"`. 			     
    * To change input filename:  
      
		`string filenameField  = "./../input/field_test.bin";`  
    	`string filenameSource = "./../input/source_test.bin";`
    	`string filenameCharge = "./../input/charge_test.bin";`  
    	`string filenameMetadata     =   "../input/metadata.txt";`  
    change them into your input filenames. 
    * To use standard kernels: 
     
  	    The same step as described in 1. ;  

6. If you want to read input from binary file, and use your own kernel:  

	Go to `/examples/`, open `"binary_file_mykernel.cpp"`.
	* To change the input filename:  
	
	  	The same step as described in 5.  
	* To define your own kernel:  
	
	   Modify `class myKernel`.
	 
When using our examples, make sure that the input file format are the same as described in  **4.**  	

####5.2 Run examples  

Here we give an example:  
If you want to use `"binary_file_standard_kernel.cpp"`  

1. As stated earlier to run the code, go to the appropriate directory and key in the following:

		make binary_file_standard_kernel

2. Make sure you have changed or cleaned the .o files from previous compilation. To clean the irrelevant files, key in:

		make clean

3. To tar the file, key in

		make tar

4. Read through the makefile for other options.

To run other .cpp files:  

1) `get_input_through_routine_myKernel.cpp`     
   key in: 

      	make get_input_through_routine_myKernel   
   
2) `get_input_through_routine_standard_kernel.cpp`     
   key in:      

   		make get_input_through_routine_standard_kernel
   
3) `binary_file_mykernel.cpp`  
   key in:
   
   		make binary_file_mykernel
4) `binary_file_standard_kernel.cpp`  
   key in:
   
   		make binary_file_standard_kernel

