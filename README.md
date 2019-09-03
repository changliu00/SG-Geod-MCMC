# [Stochastic Gradient Geodesic MCMC Methods](http://papers.nips.cc/paper/6281-stochastic-gradient-geodesic-mcmc-methods)

[Chang Liu][changliu] (<chang-li14@mails.tsinghua.edu.cn> (deprecated); <liuchangsmail@gmail.com>),
[Jun Zhu][junzhu], and [Yang Song][yangsong]. NIPS 2016.

[Paper](http://ml.cs.tsinghua.edu.cn/~changliu/sggmcmc-sam/sggmc_nips2016.pdf)
[Appendix](http://ml.cs.tsinghua.edu.cn/~changliu/sggmcmc-sam/sggmc_supp_nips2016.pdf)
[Slides](http://ml.cs.tsinghua.edu.cn/~changliu/sggmcmc-sam/sggmc_beamer_nips2016.pdf)
[Poster](http://ml.cs.tsinghua.edu.cn/~changliu/sggmcmc-sam/sggmc_poster_nips2016.pdf)

## Introduction

The repository implements the proposed methods, Stochastic Gradient Geodesic Monte Carlo (SGGMC)
and geodesic Stochastic Gradient Nose-Hoover Thermostats (gSGNHT), and their application on the
Spherical Admixture Model (SAM) ([Reisinger et al., 2010](https://icml.cc/Conferences/2010/papers/45.pdf)),
where the accuracy and scalability of SGGMC and gSGNHT can be demonstrated.
Implementations of baseline methods, including Geodesic Monte Carlo (GMC)
([Byrne & Girolami, 2013](https://onlinelibrary.wiley.com/doi/full/10.1111/sjos.12036)),
variational inference (VI) ([Reisinger et al., 2010](https://icml.cc/Conferences/2010/papers/45.pdf))
and its stochastic version (StoVI) ([Hoffman et al., 2013](http://jmlr.org/papers/v14/hoffman13a.html)),
are also provided.

* Codes for MCMC methods (SGGMC, gSGNHT, GMC) are created by Chang Liu
  (except for the [Eigen library](http://eigen.tuxfamily.org/)).
  They are implemented in C++ with [OpenMP](https://www.openmp.org/) paralellization.
* Codes for variational inference methods (VI, StoVI) are modified by Chang Liu based on
  the codes provided by Reisinger et al. (<joeraii@cs.utaxes.edu>).
  They are implemented in MATLAB.
* "20News-diff" dataset is available at <http://ml.cs.tsinghua.edu.cn/~changliu/sggmcmc-sam/>.
  It is generated from the Matlab/Octave version of
  the [20Newsgroups dataset](http://www.qwone.com/~jason/20Newsgroups/).
  Check the "README.md" file in the data package for more details.
* "Wikipedia-150K" dataset is available at <http://ml.cs.tsinghua.edu.cn/~changliu/sggmcmc-sam/>.
  It is generated from the processed dataset by
  [Aonan Zhang](http://ml.cs.tsinghua.edu.cn/~aonan/datasets/wikipedia/).
  Check the "README.md" file in the data package for more details.

## Instructions

### MCMC methods (SGGMC, gSGNHT, GMC)

Corresponds to the folder "MCMC/". Implemented by Chang Liu in C++.

* Introduction:
	- "MCMC/SAM_GMC/" contains implementations of GMC-apprMH and GMC-bGibbs.
	- "MCMC/SAM_SGGMC/" and "MCMC/SAM_GSGNHT/" contain implementations
	  of SGGMC and gSGNHT, respectively.
	- "MCMC/myinclude/" contains basic common utility files,
	  including the [Eigen library](http://eigen.tuxfamily.org/) for linear algebra.

* Compiling:
	- Please ensure that you are on a Linux machine that has installed g++
	  (or other C++ compilers supporting [OpenMP](http://openmp.org/) and
	  the "Makefile" files are modified for the available compiler).
	- Execute the following commands:
	  ```
		cd MCMC/myinclude/
		make
		cd ..
		make
		cd SAM_GMC/
		make
		cd ../SAM_SGGMC/
		make
		cd ../SAM_GSGNHT/
		make
	  ```

* Usage:
	- First download and unzip our provided data, to the current folder.
	- Set the number of OpenMP threads (we use 24) by the console command:
	  ```
		export OMP_NUM_THREADS=24
	  ```
	- To run inference methods on the 20Newsgroups-different dataset,  
	  -----------------------------------------------------------------------------------
		To run method:	| change directory to:	| execute:
	  ------------------+-----------------------+----------------------------------------
		GMC-apprMH		| "MCMC/SAM_GMC/"		| `./samgmc tr settings_diff_approxMH.txt`
		GMC-bGibbs		| "MCMC/SAM_GMC/"		| `./samgmc tr settings_diff_bGibbs.txt`
		SGGMC-batch		| "MCMC/SAM_SGGMC/"		| `./samsggmc tr settings_diff_batch.txt`
		SGGMC-full		| "MCMC/SAM_SGGMC/"		| `./samsggmc tr settings_diff_full.txt`
		gSGNHT-batch	| "MCMC/SAM_GSGNHT/"	| `./samgsgnht tr settings_diff_batch.txt`
		gSGNHT-full		| "MCMC/SAM_GSGNHT/"	| `./samgsgnht tr settings_diff_full.txt`
	  ------------------+-----------------------+----------------------------------------
	  To run on the 150K subset of Wikipedia, just replace "diff" with "wiki" in the above table.
	  
	  Above executions will reproduce our results. To change parameters or settings,
	  you can either edit the settings file, or provide a new parameter while executing:
	  e.g. run `./samgmc tr settings_diff_approxMH.txt btEps 1e-5 ts_tmIntv 300`.
	  The character "%" in the settings file will comment the rest of the line after it.
	  Detailed format instructions of settings file can be found in
	  "MCMC/myinclude/myutils_3.hpp", lines 33-41.

	- Output of the programs:  
		+ A folder with name concatenated by the provided prefix, the date and hour,
		  and a number distinguishing names within the hour.

		In the folder:  
		+ Settings for this run is written.
		+ Successive samples of model hidden variables.
		+ Time of checkpoints.
		+ A folder "logperp_[distinguishing_number]/", with a result file
		  and files recording the settings for evaluating the log-perplexity.
		  The four columns in the results file are: burn-in size, total sample numbers,
		  elapsed wall time, and log-perplexity.

	- Other modes of the programs:
		+ `./samgmc ts [modelDirname] [settingsFilename] ([var]) ([val]) (...)`  
		  Evaluate the log-perplexity of samples contained in the directory "modelDirname".

		+ `./samgmc tw [modelDirname] [settingsFilename] ([var]) ([val]) (...)`  
		  Get top-words of the averaged topic samples.

		+ `./samgmc re [tsDirname]`  
		  Go on training the model. "tsDirname" is of the form "modelDirname/logperp_n",
		  so the "modelDirname" is then known by the program, as well as evaluation settings
		  and previous evaluations.

### Variational methods (VI, StoVI)

Corresponds to the folder "Variational/".
Modified by Chang Liu based on the MATLAB codes by Reisinger et al. (<joeraii@cs.utaxes.edu>).

* Introduction:
	- "Variational/VI/" and "Variational/StoVI/" correspond to the vanilla mean-field variational
	  inference method ([Reisinger et al., 2010](https://icml.cc/Conferences/2010/papers/45.pdf)) and
	  its stochastic version based on [Hoffman et al. (2013)](http://jmlr.org/papers/v14/hoffman13a.html)
	  for scalability.
	  Usage is the same in both folders.
	- Files with name containing the substring "CompWithSmp" ("for comparison with sampling methods")
	  are our added files.
	- File "Variational/StoVI/VariationalStoEM_CompWithSmp.m" is modified according to the
	  Stochastic Variational Inference framework of [Hoffman et al. (2013)](http://jmlr.org/papers/v14/hoffman13a.html).

* Usage:
	- Same in both "Variational/VI/" and "Variational/StoVI/" folders.
	- "CompWithSmp_diff.m" and "CompWithSmp_wiki.m" are the entrance for running
	  variational inference on "20News-different" and "150K-Wikipedia" datasets, respectively.
	- Perplexity evaluation for trained models is done in "MCMC/SAM_GMC/", using the command:
	  ```
	    ./samgmc ts [model_filename] settings_[diff/wiki]_vi.txt [other_options]
	  ```

## Citation
```
	@incollection{liu2016stochastic,
	  title = {Stochastic Gradient Geodesic {MCMC} Methods},
	  author = {Liu, Chang and Zhu, Jun and Song, Yang},
	  booktitle = {Advances in Neural Information Processing Systems 29},
	  editor = {D. D. Lee and M. Sugiyama and U. V. Luxburg and I. Guyon and R. Garnett},
	  pages = {3009--3017},
	  year = {2016},
	  publisher = {Curran Associates, Inc.},
	  url = {http://papers.nips.cc/paper/6282-stochastic-gradient-geodesic-mcmc-methods.pdf},
	  organization={NIPS Foundation},
	  address={Barcelona, Spain}
	}
```

[changliu]: http://ml.cs.tsinghua.edu.cn/~changliu/index.html
[junzhu]: http://ml.cs.tsinghua.edu.cn/~jun/index.shtml
[yangsong]: https://yang-song.github.io/

