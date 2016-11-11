# SEM_CPP
C++ implementation of the SEM algorithm

#Tasks

Determine how to interface Bio::DB::BigWig module in original code with the current code.

Continue C++ implementation of algorithm.
	May use a single header file or multiple, but plan to equate
	each .cpp file with each .pl file
	Use descriptive variable names
	ALWAYS INITIALIZE VARIABLES
		can be quite difficult to debug uninitialized variables
	Use compiler flags as provided in Makefile
	
#Notes
	
functions used within Bio:DB::BigWig
	example perl lines within accum*.pl
		my $wig = new(-bigwig=>$hfile);
		my @features = $wig->features(-seq_id=>$seqid, -type=>'bin'.$total_size,-start=>$upstart,-end=>$upend);

Should only need to "translate" the two above functions.
Only "translate" the single argument version of new.
	
#Resources

http://search.cpan.org/~lds/Bio-BigFile/lib/Bio/DB/BigWig.pm
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2922891/

