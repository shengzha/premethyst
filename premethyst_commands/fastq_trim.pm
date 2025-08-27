package premethyst_commands::fastq_trim;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_trim");

sub fastq_trim {

# defaults
$min_RL = 30;
$threads = 1;
$r2_trim = 0;
$a1 = "CTATCTCTTATA";
$a2 = "AGATCGGAAGAGC";

getopts("O:1:2:m:t:e:ua:b:U:", \%opt);

$die = "

premethyst fastq-trim (options) -O OutputPrefix (-1 read1.fq.gz -2 read2.fq.gz) OR (-U read1.fq.gz)
       or  trim

Runs trim galore with settings specified for sciMETv2/3.
Uses sciMETv2/3 adapter sequences (hard coded) and many
of the trim galore defaults in paired mode.

Trimming will complete and produce a 'trim.complete'
file, at which point subsequent processing can happen.
The command will continue to run, producing additional
QC stats that are not required prior to alignment.

Can be run for sciMETv3 Illumina Mode (-1 and -2) or
with the Ultima sequencing sciMETv3 design (-U).

Options:
-O   [STR]   Output Prefix (req)
-1   [STR]   Read1 fastq (req if Illumina)
-2   [STR]   Read2 fastq (req if Illumina)
-U   [STR]   Ultima unpaired read (req if Ultima)
-m   [INT]   Min read length (def = $min_RL)
-t   [INT]   Threads to use (def = $threads)
-e   [INT]   Trim bases from the end of read 2 after adapter trim (def = $r2_trim)
               (set to 0 for SL preps since Hmer is not incorporated, 10 for LA)
-u           Retain unpaired reads (Illumina; def = discard)

-a   [STR]   Adapter 1 sequence (Illumina only; def = $a1)
-b   [STR]   Adapter 2 sequence (Both; def = $a2)

Executable Commands (from $DEFAULTS_FILE)
   trim_galore: $trim_galore
   zcat:        $zcat
   pigz:        $pigz

";

# universal options

if (!defined $opt{'O'}) {die "\nERROR: Specify output prefix as -O\n$die"};
if (defined $opt{'t'}) {$threads = $opt{'t'}};
if (defined $opt{'e'}) {$r2_trim = $opt{'e'}};
if (defined $opt{'a'}) {$a1 = $opt{'a'}};
if (defined $opt{'b'}) {$a2 = $opt{'b'}};
if (defined $opt{'U'} && (defined $opt{'1'} || defined $opt{'2'})) {
	die "\nEROR: Cannot specify Ultima reads (-U) and Illumina reads (-1 or -2) in the same command!\n$die";
}

if (!defined $opt{'U'}) { # Illumina mode, paired reads
	if (!defined $opt{'1'} || !defined $opt{'2'}) {die "\nERROR: Reads 1 and 2 must be specified.\n$die"};

	$a_opts = (defined $opt{'a'} && defined $opt{'b'}) ? "-a $a1 -a2 $a2" : "";

	$trim_command = "$trim_galore $a_opts -j $threads --paired";
	$trim_command .= " --three_prime_clip_R2 $r2_trim" if $r2_trim > 0;
	$trim_command .= " --retain_unpaired" if defined $opt{'u'};
	$trim_command .= " $opt{'1'} $opt{'2'} >> $opt{'O'}.trim.log 2>> $opt{'O'}.trim.log";
} else { # Ultima Mode, single-end
	$trim_command = "$trim_galore -a $a2 -j $threads $opt{'U'}";
	$trim_command .= " --three_prime_clip_R2 $r2_trim" if $r2_trim > 0;
	$trim_command .= " >> $opt{'O'}.trim.log 2>> $opt{'O'}.trim.log";
}

system("echo 'RUNNING: $trim_command' >> $opt{'O'}.trim.log");
system($trim_command);

# rename to output prefix names
# trim_galore output is written to the current directory, need to remove the path from input file names to get the basename for output files
use File::Basename;
if (defined $opt{'1'}) { # Illumina
	$out1 = basename($opt{'1'}, '.fq.gz'); $out1 .= "_val_1.fq.gz";
	$out2 = basename($opt{'2'}, '.fq.gz'); $out2 .= "_val_2.fq.gz";
	$trimmed_r1 = "$opt{'O'}.trimmed.paired.R1.fq.gz";
	$trimmed_r2 = "$opt{'O'}.trimmed.paired.R2.fq.gz";
	system("mv $out1 $trimmed_r1");
	system("mv $out2 $trimmed_r2");

	if (defined $opt{'u'}) {
		$out1u = basename($opt{'1'}, '.fq.gz'); $out1u .= "_unpaired_1.fq.gz";
		$out2u = basename($opt{'2'}, '.fq.gz'); $out2u .= "_unpaired_2.fq.gz";
		$trimmed_unpaired_r1 = "$opt{'O'}.trimmed.unpaired.R1.fq.gz";
		$trimmed_unpaired_r2 = "$opt{'O'}.trimmed.unpaired.R2.fq.gz";
		system("mv $out1u $trimmed_unpaired_r1");
		system("mv $out2u $trimmed_unpaired_r2");
	}
} else { # Ultima
	$out = basename($opt{'U'}, '.fq.gz'); $out .= "_val.fq.gz";
	$trimmed_r1 = "$opt{'O'}.trimmed.fq.gz";
	system("mv $out $trimmed_r1");
}

# write trimming complete file to note further processing can happen
system("date > $opt{'O'}.trim.complete");


# get stats on the trim in barcode-aware manner.
$zcat = "$pigz -dc -p $threads" if ($threads > 1);
$untrimmed_r1 = (defined $opt{'1'}) ? $opt{'1'} : $opt{'U'};

%BARC_IN_ct = (); %BARC_OUT_ct = (); %BARC_R1_up = (); %BARC_R2_up = ();

open IN, "$zcat $untrimmed_r1 |";
while ($tag = <IN>) {
	chomp $tag; $tag =~ s/^@//; $tag =~ s/:.+$//;
	$BARC_IN_ct{$tag}++;
	$BARC_OUT_ct{$tag} = 0;
	$null = <IN>; $null = <IN>; $null = <IN>;
} close IN;

open IN, "$zcat $trimmed_r1 |";
while ($tag = <IN>) {
	chomp $tag; $tag =~ s/^@//; $tag =~ s/:.+$//;
	$BARC_OUT_ct{$tag}++;
	$null = <IN>; $null = <IN>; $null = <IN>;
} close IN;

if (defined $trimmed_unpaired_r1) {
	open IN, "$zcat $trimmed_unpaired_r1 |";
	while ($tag = <IN>) {
		chomp $tag; $tag =~ s/^@//; $tag =~ s/:.+$//;
		$BARC_R1_up{$tag}++;
		$null = <IN>; $null = <IN>; $null = <IN>;
	} close IN;

	open IN, "$zcat $trimmed_unpaired_r2 |";
	while ($tag = <IN>) {
		chomp $tag; $tag =~ s/^@//; $tag =~ s/:.+$//;
		$BARC_R2_up{$tag}++;
		$null = <IN>; $null = <IN>; $null = <IN>;
	} close IN;
}

open RPT, ">$opt{'O'}.trimmed.stats.txt";
foreach $barc (keys %BARC_IN_ct) {
	$pct = sprintf("%.2f", ($BARC_OUT_ct{$barc}/$BARC_IN_ct{$barc})*100);
	if (!defined $BARC_R1_up{$barc}) {$BARC_R1_up{$barc} = 0};
	if (!defined $BARC_R2_up{$barc}) {$BARC_R2_up{$barc} = 0};
	$pct1u = sprintf("%.2f", ($BARC_R1_up{$barc}/$BARC_IN_ct{$barc})*100);
	$pct2u = sprintf("%.2f", ($BARC_R2_up{$barc}/$BARC_IN_ct{$barc})*100);
	if (defined $trimmed_unpaired_r1) {
		print RPT "$barc\t$BARC_IN_ct{$barc}\t$BARC_OUT_ct{$barc}\t$pct\t$BARC_R1_up{$barc}\t$pct1u\t$BARC_R2_up{$barc}\t$pct2u\n";
	} else {
		print RPT "$barc\t$BARC_IN_ct{$barc}\t$BARC_OUT_ct{$barc}\t$pct\n";
	}
} close RPT;

exit;

}

1;