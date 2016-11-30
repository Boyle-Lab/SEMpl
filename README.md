# SEM_CPP
C++ implementation of the SEM algorithm

#Compilation
Copy entire repository to computer, then cd into lib/libBigWig-master/ and run "make"

cd back into directory where makefile is present, or from the previous location "cd .." twice

#Tasks

Determine how to interface Bio::DB::BigWig module in original code with the current code.

Continue C++ implementation of algorithm.
	May use a single header file or multiple, but plan to equate
	each .cpp file with each .pl file
	Use descriptive variable names
	ALWAYS INITIALIZE VARIABLES
		can be quite difficult to debug uninitialized variables
	Use compiler flags as provided in Makefile
	Will use struct as opposed to class for simplicity
	Will use public data members for simplicity of struct
	
#Notes

11-16-2016
In process of finding a C++/C interface for BigWig files/software, as opposed to attempting to rewrite perl module,

Below notes are now irrelevant

Old below,	
functions used within Bio:DB::BigWig

	example perl lines within accum*.pl
		my $wig = new(-bigwig=>$hfile);
		my @features = $wig->features(-seq_id=>$seqid, -type=>'bin'.$total_size,-start=>$upstart,-end=>$upend);

Should only need to "translate" the two above functions.

Only "translate" the single argument version of new.
	
#Resources

http://search.cpan.org/~lds/Bio-BigFile/lib/Bio/DB/BigWig.pm

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2922891/

For tfmpvalue binary source code, use below

tfmpvalue.cpp Bonsai

Possible BigWig file/software interface below

https://github.com/deltadev/bbi

github.com/jayhesselberth/libBigWig

#Notes regarding Bio::DB::BigWig replacement
functions used from Bio::DB::BigWig in original implementation of algorithm
	new(-bigwig=>file); equilvalent to new("-bigwig", $hfile);
		creates new object of Bio::DB::BigWig type, $hfile points to the indexed .bw file

	features("-seq_id", $seqid, "-type", 'bin:'.$total_size, "-start", $upstart, "-end", $upend);
		$seqid is chromosome or contig name defining the range of interest
		'bin:' + $total_size is the type of feature to retrieve
		$upstart is the start of the range of interest
		$upend is the end of the range of interest
	pseudo-result below
		fetch $total_size intervals across region $upstart to $upend on chromosome $seqid
