#/usr/bin/perl -w

########Setting System Parameters###################

$Threads_num = 1;
$Mallet_path = "/Mallet/bin";
$Iterations_num = 500;

########Parsing User Parameters#####################

$MTX_path = "";
$Feature_path = "";
$Barcode_path = "";
$Out_prefix = "";
$Topics_num = 0;
while(my $parameter = shift){
	if($parameter eq "-input"){
		$parameter = shift;
		if($parameter=~/^[\p{L}\p{N}\p{P}\/ _]+$/ and $parameter!~/^-/){
			if(-f "$parameter/matrix.mtx.gz"){
				$MTX_path = "$parameter/matrix.mtx.gz";
			}elsif(-f "$parameter/matrix.mtx"){
				$MTX_path = "$parameter/matrix.mtx";
			}else{
				print STDERR "No expression matrix file found!\nThe file name should be 'matrix.mtx(.gz)' under the MEX directory.\n";
				last;
			}
			if(-f "$parameter/barcodes.tsv.gz"){
				$Barcode_path = "$parameter/barcodes.tsv.gz";
			}elsif(-f "$parameter/barcodes.tsv"){
				$Barcode_path = "$parameter/barcodes.tsv";
			}elsif(-f "$parameter/cells.tsv.gz"){
				$Barcode_path = "$parameter/cells.tsv.gz";
			}elsif(-f "$parameter/cells.tsv"){
				$Barcode_path = "$parameter/cells.tsv";
			}else{
				print STDERR "No barcode/cell list file found!\nThe file name should be 'barcodes.tsv(.gz)' under the MEX directory.\n";
				last;
			}
			if(-f "$parameter/features.tsv.gz"){
				$Feature_path = "$parameter/features.tsv.gz";
			}elsif(-f "$parameter/features.tsv"){
				$Feature_path = "$parameter/features.tsv";
			}elsif(-f "$parameter/genes.tsv.gz"){
				$Feature_path = "$parameter/genes.tsv.gz";
			}elsif(-f "$parameter/genes.tsv"){
				$Feature_path = "$parameter/genes.tsv";
			}else{
				print STDERR "No feature/gene list file found!\nThe file name should be 'genes.tsv(.gz)' under the MEX directory.\n";
				last;
			}
		}else{
			print STDERR "Invalid MEX path!\n";
			last;
		}
	}elsif($parameter eq "-output"){
		$parameter = shift;
		if($parameter=~/^[\p{L}\p{N}\p{P}\/ _]+$/ and $parameter!~/^-/){
			$Out_prefix = $parameter;
		}else{
			print STDERR "Invalid output prefix!\n";
			last;
		}
	}elsif($parameter eq "-topics"){
		$parameter = shift;
		$parameter =~ s/^0+//;
		if($parameter=~/^\p{N}+$/){
			$Topics_num = $parameter;
		}else{
			print STDERR "Invalid topic number!\n";
			last;
		}
	}elsif($parameter eq "-help" || $parameter eq "-h" || $parameter eq "--help"){
		last;
	}else{
		print STDERR "Illegle parameter \"$parameter\"!\n";
		last;
	}
}
if($MTX_path eq "" or $Barcode_path eq "" or $Feature_path eq "" or $Out_prefix eq "" or $Topics_num == 0){
	print STDERR 'Usage: perl run_lda.pl -input Path_of_MEX_directory -output Prefix_of_results -topics Number_of_topics

The program recieves 3 parameters: -input, -output, and -topics

-input DIRECTORY 
The Path_of_MEX_directory is the "CellRanger Count" output direcotry, must contain 3 files:
(1) features.tsv(.gz) or genes.tsv(.gz)
(2) barcodes.tsv(.gz) or cells.tsv(.gz)
(3) matrix.mtx(.gz)
REQUIRED

-output FILENAME
The results include 3 files:
(1) prefix.c2k.txt - Cell-PF(Document-Topic) matrix
(2) prefix.k2g.txt - PF-Gene(Topic-Term) matrix
(3) prefix.mallet.log - Log for the Mallet run
REQUIRED

-topics INTEGER
Topics number for LDA training
REQUIRED

------------------------------------------------------------------------
Several global parameters need to be reset before running this program for the first time on a new device.

Threads_num = INTEGER
How many threads can be occupied during the LDA training.
Default is 1.

Mallet_path = DIRECTORY
The pre-installed Mallet program location.
Default is /Mallet/bin.

Iterations_num = INTEGER
Number of iterations for Gibbs sampling.
Default is 500.
The global parameters can be found at the top of the PERL script.

';
	exit;
}

########Preparing Input Data for Mallet#############


print STDERR localtime()." | Loading barcode/cell and feature/gene information...";

if($Barcode_path=~/gz$/){
	open(CELL,"gzip -dc $Barcode_path |") || die "Cannot read $Barcode_path.\n";
}else{
	open(CELL,"<$Barcode_path") || die "Cannot read $Barcode_path.\n";
}
@cells = <CELL>;
close(CELL);

if($Feature_path=~/gz$/){
	open(GENE,"gzip -dc $Feature_path |") || die "Cannot read $Feature_path.\n";
}else{
	open(GENE,"<$Feature_path") || die "Cannot read $Feature_path.\n";
}
@genes = <GENE>;
close(GENE);

%genes = [];
for(my $i=0;$i<@genes;$i++){
	chomp($genes[$i]);
	$genes[$i] =~ s/\t.*$//;
	$genes{lc($genes[$i])} = $i;
}

print STDERR "Done.\n";

print STDERR localtime()." | Loading expression matrix.\n";

if($MTX_path=~/gz$/){
	open(EXPR,"gzip -dc $MTX_path |") || die "Cannot read $MTX_path.\n";
}else{
	open(EXPR,"<$MTX_path") || die "Cannot read $MTX_path.\n";
}
while(<EXPR>){
	next if(/^%/);
	chomp;
	my @temp = split(/ /);
	if($temp[0] != scalar(@genes) || $temp[1] != scalar(@cells)){
		print STDERR "The number of barcodes/cells or features/genes do not match the expression matrix!\n";
		exit;
	}
	last;
}

print STDERR localtime()." | ".scalar(@cells)." cells x ".scalar(@genes)." genes are going to be processed...";
my $current_cell = -1;
my $corpus_buffer = "";
open(CORPUS,">$Out_prefix.corpus.txt") || die "Cannot write to $Out_prefix.corpus.txt.\n";
while(<EXPR>){
	chomp;
	my @temp = split(/ /);
	if($temp[1] != $current_cell){
		if($current_cell != -1){
			chop($corpus_buffer);
			print CORPUS $corpus_buffer."\n";
		}
		$current_cell = $temp[1];
		$corpus_buffer = "";
	}
	$corpus_buffer.="$genes[$temp[0]-1] " x $temp[2];
}
if($current_cell != -1){
	chop($corpus_buffer);
	print CORPUS $corpus_buffer."\n";
}
close(EXPR);
close(CORPUS);
print STDERR "Done.\n";

########Running Mallet##############################

print STDERR localtime()." | Modelling with Mallet.\n";
if(-f -x "$Mallet_path/mallet"){
	if(-f "$Out_prefix.corpus.txt"){
		print STDERR localtime()." | Importing data...";
		`$Mallet_path/mallet import-file --input $Out_prefix.corpus.txt --output $Out_prefix.mallet --keep-sequence --token-regex '[\\p{L}\\p{N}\\p{P}_]+' 2>$Out_prefix.mallet.log`;
		print STDERR "Done.\n";
		`rm $Out_prefix.corpus.txt`;
	}else{
		print STDERR "Failed to import $Out_prefix.corpus.txt!\n";
	}

	if(-f "$Out_prefix.mallet"){
		if($Threads_num > 1){
			print STDERR localtime()." | Parallel Gibbs sampling with $Threads_num threads...";
		}else{
			print STDERR localtime()." | Single thread Gibbs sampling...";
		}
		`$Mallet_path/mallet train-topics --input $Out_prefix.mallet --num-topics $Topics_num --output-doc-topics $Out_prefix.c2k.temp --topic-word-weights-file $Out_prefix.k2g.temp --num-threads $Threads_num --num-iterations $Iterations_num 2>>$Out_prefix.mallet.log`;
		print STDERR "Done.\n";
		`rm $Out_prefix.mallet`;
		print STDERR localtime()." | Modelling Completed.\n";
	}else{
		print STDERR "Failed to generate $Out_prefix.mallet!\n";
	}
}else{
	print STDERR "Invalid Mallet path: $Mallet_path!\n";
}

########Post Processing#############################

if(-f "$Out_prefix.c2k.temp"){
	print STDERR localtime()." | Processing parameter phi...";
	open(C2KI,"<$Out_prefix.c2k.temp") || die "Cannot read $Out_prefix.c2k.temp!\n";
	open(C2KO,">$Out_prefix.c2k.txt") || die "Cannot write to $Out_prefix.c2k.txt!\n";
	while(<C2KI>){
		s/^.*?\t.*?\t//;
		print C2KO $_;
	}
	close(C2KI);
	close(C2KO);
	print STDERR "Done.\n";
	`rm $Out_prefix.c2k.temp`;
}

if(-f "$Out_prefix.k2g.temp"){
	my @topic_gene = [];
	for(my $i=0;$i<$Topics_num;$i++){
		$topic_gene[$i] = [];
		for(my $j=0;$j<@genes;$j++){
			$topic_gene[$i][$j] = 0;
		}
	}
	print STDERR localtime()." | Processing parameter beta...";
	open(K2GI,"<$Out_prefix.k2g.temp") || die "Cannot read $Out_prefix.k2g.temp!\n";
	my @melt = <K2GI>;
	close(K2GI);
	my $temp_topic = 0;
	my $temp_sum = 0;
	for(my $i=0;$i<@melt;$i++){
		chomp($melt[$i]);
		my @temp = split(/\t/,$melt[$i]);
		if($temp[0] ne $temp_topic){
			for(my $j=0;$j<@genes;$j++){
				$topic_gene[$temp_topic][$j] = $topic_gene[$temp_topic][$j]/$temp_sum;
			}
			$temp_sum = 0;
			$temp_topic = $temp[0];
		}
		$topic_gene[$temp[0]][$genes{lc($temp[1])}] = $temp[2];
		$temp_sum += $temp[2];
	}
	for(my $j=0;$j<@genes;$j++){
		$topic_gene[$temp_topic][$j] = $topic_gene[$temp_topic][$j]/$temp_sum;
	}
	undef @melt;

	open(K2GO,">$Out_prefix.k2g.txt") || die "Cannot write to $Out_prefix.k2g.txt!\n";
	for(my $i=0;$i<$Topics_num;$i++){
		print K2GO join("\t",@{$topic_gene[$i]})."\n";
	}
	close(K2GO);
	print STDERR "Done.\n";
	`rm $Out_prefix.k2g.temp`;
}

if(-f "$Out_prefix.k2g.txt" and -f "$Out_prefix.c2k.txt"){
	print STDERR localtime()." | All finished.\n";
}else{
	print STDERR localtime()." | Modelling aborted!\n";
}

__END__
