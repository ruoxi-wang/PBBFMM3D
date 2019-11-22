# PBBFMM3D  

### 1. INTRODUCTION

PBBFMM3D is a parallel open source package of the Black-box Fast Multipole Method in 3 dimensions.   
The Black-box Fast Multipole Method is an O(N) fast multipole method, which is a technique to calculate sums of the form  

![](http://latex.codecogs.com/gif.latex?f%28x_i%29%20%3D%20%5Cdisplaystyle%20%5Csum_%7Bj%3D1%7D%5EN%20K%28x_i%2Cy_j%29%20%5Csigma_j%2C%20%5C%2C%5C%2C%5C%2C%20%5Cforall%20i%20%5Cin%5C%7B1%2C2%2C%5Cldots%2CN%5C%7D)

where ![](http://latex.codecogs.com/gif.latex?K%28x_i%2Cx_j%29) is kernel function, ![](http://latex.codecogs.com/gif.latex?x_i) are target points, ![](http://latex.codecogs.com/gif.latex?y_j) are srouce points, and ![](http://latex.codecogs.com/gif.latex?%5Csigma_i) are weights at corresponding locations.
PBBFMM3D provides an O(N) solution to matrix-vector products of the type Ax. In that case the relation between A and K is:
![](http://latex.codecogs.com/gif.latex?A_%7Bij%7D%20%3D%20K%28x_i%2Cy_j%29)

PBBFMM3D

*	applies to all non-oscillatory smooth kernels (e.g., RBF kernels, 1/r, log r, etc.) and only requires the kernel's numerical evaluation at data points.	
*	is parallelized with OpenMP for shared memory machines.
*	is efficient for both low and high accuracies (a choice of Chebyshev and uniform interpolations).
*  pre-computes and compresses far-field translation operators.
*  applies to multiple sets of weights in one pass.

Please cite the following papers if you use this code:

William Wang and Eric Darve. "The black-box fast multipole method." Journal of Computational Physics 228, no. 23 (2009): 8712-8725. <a href="http://www.sciencedirect.com/science/article/pii/S0021999109004665">link</a>

Ruoxi Wang, Chao Chen, Jonghyun Lee and Eric Darve. "PBBFMM3D: a parallel black-box method for kernel matrix-vector multiplication". <a href="https://arxiv.org/abs/1903.02153">arXiv</a>

### 2. DIRECTORIES AND FILES


	./examples/		:	Example input C++ codes
	./src/			:	Source code in C++  
	./include/		:	Relevant header files  
	./exec/			:	Executables for PBBFMM3D  
	./input/		:	The input file  
	./output/		:   Output M2L operators
	./README.md		:	This file  
	./License.md	:	License file  
	./Makefile		:	Makefile
	./python        :   Python interface for PBBFMM3D and REIG
	
### 3. TUTORIAL
#### 3.1 To Get Started  
To check whether things are set up correctly, you can perform the following: Go to the directory where Makefile is in, then key in the following three commands in the terminal:

		make binary_file_mykernel
		cd exec/
		./binary_file_mykernel

#### 3.2 Basic usage

##### 3.2.1 BBFMM3D with standard kernel

The basic usage of BBFMM3D with standard kernel is as follows: 

	#include"bbfmm3d.hpp"  
	...
	{
	    
    double L;                   // Length of simulation cell (assumed to be a cube)
    int interpolation_order;    // Number of interpolation nodes per dimension
    int Ns;                     // Number of sources in simulation cell
    int Nf;                     // Number of targes in simulation cell
    int nCols;                  // Number of columns for weights
    int tree_level;             // Tree level
    double eps = 1e-5 ;         // Target accuracy (SVD)
    int use_chebyshev = 1;      // 1: chebyshev interpolation; 0: uniform interpolation    
    
    std::vector<vector3> source(Ns); // Position array for the source points
    std::vector<vector3> target(Nf);  // Position array for the target points
    std::vector<double> weight(Ns*nCols); // Weight 
    std::vector<double> output(Nf*nCols);

       …
	kernel_Laplacian Atree(L, tree_level, interpolation_order, eps, use_chebyshev);
    Atree.buildFMMTree();  // Build the fmm tree;

	/* The following can be repeated with different field, source, and q */
    
    H2_3D_Compute<kernel_Laplacian> compute(Atree, target, source, weight, nCols, output);
    ...
    }
    
This example first build a FMM tree with these two lines:  

	kernel_Laplacian Atree(L, tree_level, interpolation_order, eps, use_chebyshev);
    Atree.buildFMMTree(); 
 
where kernel_Laplacian is a class of fmm tree using Laplacian kernel, the constructor takes 5 arguments:  

* L(double):   
	Length of simulation cell (assumed to be a cube).
* tree_level(int):  
	The number of levels in the hierarchy tree
* interpolation_order(int):  
	Number of interpolation nodes per dimension. 
* eps(double):  
	Target accuracy; this is used to determine which
 singular values are kept after the SVD. 
* use_chebyshev(int):  
	Label to indicate whether using Chebyshev interpolation formula or uniform interpolation formula.  
	use_chebyshev = 1: Chebyshev interplation formula;  
	use_chebyshev = 0: uniform interpolation formula (where FFT is used)
 
Once the tree is created, you can repeat matrix-vector product with different target, source and weight.(see **3.2.4**) The code shows an example using Lapacian kernel:  

	H2_3D_Compute<kernel_Laplacian> compute(Atree, target, source, weight, nCols, output);
The template class `H2_3D_Compute` is for computing matrix-vector product after the FMM tree is built. The constructor takes 8 arguments:  

* Atree (T):  
	The FMM tree we just built. Here T is typename, FMM tree built with different kernels have different typenames.
* target (std::vector<vector3>):  
	Position array of target points. vector3 is a struct storing x, y and z coordinates of a 3D point.  
* source (std::vector<vector3>):  
	Position array of source points.   
* weight (vector<double>):  
	Array of weights.
* nCols (int):   
	Number of sets of weights.  
* output (vector<double>):   
	The output, which is stored column-wise.

##### 3.2.2 Options of provided kernels

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
	
	
    		
If you want to define your own kernel, please see **3.2.3**.

##### 3.2.3 PBBFMM3D with customized kernel

The basic usage is almost the same as **3.2.1** except that you have to define your own routine of computing kernel. You need to fill in two pieces of information: the kernel definition and the kernel's homogeneous and symmetric properties. One example is as follows:  
	
	#include "bbfmm3d.hpp"


	class myKernel: public H2_3D_Tree {
	public:
	    myKernel(double L, int tree_level, int interpolation_order,  double epsilon, int use_chebyshev):H2_3D_Tree(L,tree_level,interpolation_order, epsilon, use_chebyshev){};
	    virtual void SetKernelProperty() {
	        homogen = -1;
	        symmetry = 1;
	        kernelType = "myKernel";
	    }
	    virtual double EvaluateKernel(vector3& targetpos, vector3& sourcepos) {
	        vector3 diff;        
	        diff.x = sourcepos.x - targetpos.x;
	        diff.y = sourcepos.y - targetpos.y;
	        diff.z = sourcepos.z - targetpos.z;
	        return 1. / sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
	    }
	};	

	...
	{
	…
	myKernel Atree(L, tree_level, interpolation_order, epsilon, use_chebyshev);
    Atree.buildFMMTree();  // Build the fmm tree;
	
	/* The following can be repeated with different field, source, and q */
    
    H2_3D_Compute<myKernel> compute(Atree, target, source, weight, nCols, output);
    ...
    }
 
	
You can define your own kernel inside `EvaluateKernel(vector3& targetpos, vector3& sourcepos)`, it takes target point and source point as input, and return the kernel value. 
                                                                
You also need to define information about kernel inside `SetKernelProperty()` 

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
                      
                                
                           
##### 3.2.4 Usage of multiple sources with same kernel

If you want to compute with different sources (source points, target points, weight) but with same kernel and interpolation order, you can do it in one file:
e.g.  

	class myKernel: public H2_3D_Tree {
	...
	}
	…
	
	{
	…
	/* Build FMM tree */
	myKernel Atree(L,tree_level,interpolation_order, epsilon, use_chebyshev);
    Atree.buildFMMTree(); 
	
	/* The following can be repeated with different field, source, and q */
    
    H2_3D_Compute<myKernel> compute1(Atree, target1, source1, weight1, nCol1, output1);
	H2_3D_Compute<myKernel> compute2(Atree, target2, source2, weight2, nCol2, output2);
    …
    }
    

The basic usage is already domonstrated in **3.2.1** and **3.2.3**. Once you have built the FMM tree, you can use different sources to compute the matrix-vector multiplication without rebuilding the tree. You can use either standard kernels provided ( see **3.2.2** ), or customized kernels ( see **3.2.3** )

#### 3.3 Pre-computation

The power of this package is in the pre-computing part, which pre-computes and compresses far-field translation operators for later calculations. These operators and tree-related information are stored in 3 files in the folder /output. Each time a same kenrel and interpolation order is used, it directly reads from the files. 

Note: If you are using your own kernel, make sure to either change the kernelType (used in filenames) or delete exsiting files if you changed your kernel. 

#### 3.4 Input parameter choice guidance
For interpolation order `p`, a higher p improves the accuracy but also increases the runtime and memory. For the interpolation scheme (indicated by `use_chebyshev`), when a low accuracy is requested, Chebyshev scheme is more efficient; when a high accuracy is requested, uniform scheme is more efficient. For the tree levels `tree_level`, it is often chosen such that the leaf has approximately 60 points. For the prescribed accuracy of SVD `eps`, it is chosen to be similar to the desired approximation accuracy.


#### 3.5 Test Interplation Error  
To give users an idea of how large the interplation error is, we have provided a routine of computing the interplation error. If you want to test the interplation error between a cluster A (of size length) and a cluster B, where B is in the interplation list of A, you can do the following:  

	kernel_Laplacian testTree(1/pow(2,2), 2, interpolation_order, eps, use_chebyshev);
    double errtest = testInterplationErr(&testTree, 100, 100);
    
The first line constructs a FMM tree, you can see **3.2.1** for more details of the constructor. Here we provide more explanations to the second and third arguments:  

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
	Number of target points in the cluster.
	




### 4. ROUTINES FOR INPUTING AND OUTPUTING DATA 
We have provided several routines for reading data from binary file, and writing data into binary file.	

#### 4.1 Reading meta data from text file
	
	void read_Metadata(const string& filenameMetadata, double& L, int& interpolation_order, int& Ns, int& Nf, int& nCols, int& tree_level);

The first argument, `filenameMetadata` is the filename for your meta data. `L` stores the length of simulation cell (assumed to be a cube); `interpolation_order ` stores the interpolation order; `Ns` stores the number of source points; `Nf` stores the number of target points; `nCols ` stores the number of sets of weights; `tree_level ` stores the number of levels in the hierarchy tree;

**File format:**  
 
`L, interpolation_order, Ns, Nf, nCols, tree_level `

For example:

 	1,4,800,800,1,2
 	 	
#### 4.2 Reading from binary file  

	void read_Sources(const string& filenameField, std::vector<vector3>& target, const int& Nf, const string& filenameSource, std::vector<vector3>& source, const int& Ns, const string& filenameCharge, std::vector<double>& weight, const int& nCols);

The arguments `filenameField`, `filenameSource` and `filenameCharge` are binary file names for target positions, source positions and weights respectively. `Nf`, `Ns` and `nCols` are the number of target points, the number of source points, and the number of sets of weights. The funtion outputs target points, source points and weights into `target`, `source`, and `weight` (column-wise), respectively.  

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
		
		targetpos0.x
		targetpos1.x
		…
		targetpos0.y
		targetpos1.y
		…
		targetpos0.z
		targetpos1.z
		...
		
3. Binary file for weights:  
	weights are stored column-wise:  
	
		first set of weights
		second set of weights
		...
	
		
	
####4.3 Writing into binary file  
	
	void write_Into_Binary_File(const string& filename, double* outdata, int numOfElems);  
This first argument is the filename for your output data. The second argument is a pointer to the output data, and the last argument is the number of elements in the array of your output data.  


### 5. EXAMPLES

We have provided several examples for PBBFMM3D. The files in `examples/` should be self-explanatory.
You can use our examples with your own input.
#### 5.1 Making changes to the examples for your own application

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
	
	  	The same step as described in 3.  
	* To define your own kernel:  
	
	   Modify `class myKernel`.
	 
When using our examples, make sure that the input file format are the same as described in  **4.**  	

#### 5.2 Run examples  

Here we give an example:  
If you want to use `"binary_file_standard_kernel.cpp"`  

1. As stated earlier to run the code, go to the appropriate directory and key in the following:

		make binary_file_standard_kernel

2. Make sure you have changed or cleaned the .o files from previous compilation. To clean the irrelevant files, key in:

		make clean

3. To tar the file, key in

		make tar

4. Read through the Makefile for other options.

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

