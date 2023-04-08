#/usr/bin/perl -w

########Setting System Parameters###################

$Threads_num = 1;
$Mallet_path = "/Mallet/bin";
$Iterations_num = 1000;

########Parsing User Parameters#####################

$Expr_path = "";
$Out_prefix = "";
$Topics_num = 0;
while(my $parameter = shift){
	if($parameter eq "-input"){
		$parameter = shift;
		if($parameter=~/^[\p{L}\p{N}\p{P}\/ _]+$/ and $parameter!~/^-/){
			$Expr_path = $parameter;
		}else{
			print STDERR "Invalid expression matrix path!\n";
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
                print STDERR 'Usage: perl run_lda.pl -input Path_of_expression_matrix -output Prefix_of_results -topics Number_of_topics

The program recieves 3 parameters: -input, -output, and -topics

-input FILE
Single input file: Cells as rows, genes as columns. The first row is gene ID, the first column is cell ID. 
REQUIRED

-output FILENAME
The results include 5 files:
(1) prefix.cell.txt - Cell list
(2) prefix.gene.txt - Gene list
(3) prefix.c2k.txt - Cell-LF(Document-Topic) matrix
(4) prefix.k2g.txt - LF-Gene(Topic-Term) matrix
(5) prefix.mallet.log - Log for the Mallet run
REQUIRED

-topics INTEGER
Topic number
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
Number of iterations for Gibbs sampling
Default is 1000.
The global parameters can be found at the top of the PERL script.

';
                exit;
	}else{
		print STDERR "Illegle parameter \"$parameter\"!\n";
		last;
	}
}
if($Expr_path eq "" or $Out_prefix eq "" or $Topics_num == 0){
	print STDERR 'Usage: perl run_lda.pl -input Path_of_expression_matrix -output Prefix_of_results -topics Number_of_topics

The program recieves 3 parameters: -input, -output, and -topics

-input FILE
Single input file: Cells as rows, genes as columns. The first row is gene ID, the first column is cell ID. 
REQUIRED

-output FILENAME
The results include 5 files:
(1) prefix.cell.txt - Cell list
(2) prefix.gene.txt - Gene list
(3) prefix.c2k.txt - Cell-LF(Document-Topic) matrix
(4) prefix.k2g.txt - LF-Gene(Topic-Term) matrix
(5) prefix.mallet.log - Log for the Mallet run
REQUIRED

-topics INTEGER
Topic number
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
Number of iterations for Gibbs sampling
Default is 1000.
The global parameters can be found at the top of the PERL script.

';
	exit;
}

########Preparing Input Data for Mallet#############

my @cells;
my @genes;
my %genes;

if(-f $Expr_path){
	open(EXPR,"<$Expr_path") || die "Cannot read $Expr_path.\n";

	print STDERR localtime()." | Generating gene lists...";
	my $genes = <EXPR>;
	chomp($genes);
	$genes =~ s/^\t//;
	@genes = split(/\t/,$genes);
	open(GENE,">$Out_prefix.gene.txt") || die "Cannot write to $Out_prefix.gene.txt.\n";
	print GENE join("\n",@genes)."\n";
	close(GENE);
	for(my $i=0;$i<@genes;$i++){
		$genes{lc($genes[$i])} = $i;
	}
	print STDERR "Done.\n";

	print STDERR localtime()." | Transforming expresssion matrix to corpus...";
	open(CORPUS,">$Out_prefix.corpus.txt") || die "Cannot write to $Out_prefix.corpus.txt.\n";
	while(<EXPR>){
		chomp;
		my @temp = split(/\t/);
		$cells[$#cells+1] = $temp[0];
		my $corpus_buffer = "";
		for(my $i=0;$i<@genes;$i++){
			$corpus_buffer .= "$genes[$i] " x $temp[$i+1];
		}
		chop($corpus_buffer);
		print CORPUS $corpus_buffer."\n";
	}
	close(CORPUS);
	print STDERR "Done.\n".localtime()." | ".scalar(@cells)." cells x ".scalar(@genes)." genes were loaded.\n";

	print STDERR localtime()." | Generating cell lists...";
	open(CELL,">$Out_prefix.cell.txt") || die "Cannot write to $Out_prefix.cell.txt\n";
	print CELL join("\n",@cells)."\n";
	close(CELL);
	print STDERR "Done.\n";

	close(EXPR);
}else{
	print STDERR "Invalid expression matrix!\n";
}

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
