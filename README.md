# GenomeAssembly
Final project for Berkeley CS267 parallel computing.
  
1) Build binary:  
	```
	make
	```
  
2) 2 ways of execution  
	a. run with any *{filename}.txt* file within *hw3-datasets/*  
		```
		make DATA={filename} run
		```

	b. run with any *{filename}.txt* file within *hw3-datasets/smaller/*  
		```
		make DATA={filename} runsmall
		```  
  
3) Test the correction  
	a. if run with *{filename}.txt* within *hw3-datasets/*  
		```
		make DATA={filename} check
		```

	b. if run with *{filename}.txt* within *hw3-datasets/smaller*  
		```
		make DATA={filename} checksmall
		```


