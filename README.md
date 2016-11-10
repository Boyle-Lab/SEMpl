# SEM_CPP
C++ implementation of the SEM algorithm

Tasks:
	Determine how to interface Bio::DB::BigWig module in original code 
		with the current code.
	
	Continue C++ implementation of algorithm.
		May use a single header file or multiple, but plan to equate
		each .cpp file with each .pl file
	
Notes:
	
	functions used within Bio:DB::BigWig
		example perl lines
			my $wig = new(-bigwig=>$hfile);
			my @features = $wig->features(-seq_id=>$seqid, -type=>'bin'.$total_size,-start=>$upstart,-end=>$upend);

			
