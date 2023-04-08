#/usr/bin/perl -w

$matrix_path = $ARGV[0];

if(!(-f $matrix_path)){
	print STDERR "Usage: perl transpose.pl Path_of_input_file > Output_file\n";
	exit;
}

my @genes;
my @cells;
my @cell_gene;

open(MAT,"<$matrix_path") || die "Cannot read $matrix_path.\n";
$cells = <MAT>;
chomp($cells);
@cells = split(/\t/,$cells);

while(<MAT>){
	chomp;
	my @temp = split(/\t/);
	$genes[$#genes+1] = $temp[0];
	@{$cell_gene[$#cell_gene+1]} = @temp[1..$#temp];
}
close(MAT);

print join("\t",@genes)."\n";
for(my $i=0;$i<@cells;$i++){
	my $buffer = "$cells[$i]\t";
	for(my $j=0;$j<@genes;$j++){
		$buffer.="$cell_gene[$j][$i]\t";
	}
	chop($buffer);
	print $buffer."\n";
}
