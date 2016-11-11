
use strict;
use Data::Dumper;

print "arguments given below \n @ARGV\n";

test(-bigwig => "the path");

#new(@ARGV);

sub test{
	my $self = shift;
	my $path = shift;

	if($self->retTrue("var")){
		print "it's true\n";
	}
	else{
		print "not true\n";
	}

}

sub retTrue{
	my $var = shift;
	return $var;
}

sub new{

print "first argument given is $_[0]\n";

my $self = shift;
print "self is $self\n";
print '@_ is ';
print "\n@_\n";
my %args = $_[0] =~ /^-/ ? @_ : (-bigwig=>shift);

#print "$#_\n";
#for(my $i = 0; $i < @_; $i++){
#	print $_[$i];
#	print $#_;
#}
print Dumper(%args);

}
