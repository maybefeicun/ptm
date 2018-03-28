# PTM & Sparse PTM
A java implement of PseudoDocTM and SparsePseudoDocTM.
Usage:
========
```java -jar ptm.jar [sparse] P K iter alpha1 alpha2 beta doc_path```

* `[sparse]		:	if add 'sparse', then run with Sparse PTM, othereise run with PTM.");`
* `P 			:	number of pesudo document.");`
* `K			:	number of topic.");`
* `iter  		:	number of iteration times.");`
* `alpha1		:	prior parameter alpha1.");`
* `alpha2		:	prior parameter alpha2.");`
* `beta			:	prior parameter beta.");`
* `doc_path		:	path of train file.(must be absolute path!)");`


The assembled ptm.jar and sample train data are in "./target/" .
