package premethyst_commands::fastq_align;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_align");

sub fastq_align {

getopts("O:1:2:t:o:R:XyYM:", \%opt);

$threads = 1;
$o_threads = 1;
$sort_mem = "4G";

$die = "

premethyst fastq-align (options) -R [reference path] -O [output prefix] -1 [read1.trimmed.fq.gz] (-2 [read2.trimmed.fq.gz])
      or   align

Wrapper for BSBOLT to run alignment of sciMETv2/3 reads.
reads can be a list that is comma-separated.

If only one read is specified, it assumes it is Ultima
sequencing where the Tn5 side is read 1.

Will sort output bam by read name.

Options:

-R   [STR]   Reference path (required)
    Shortcuts:
$ref_shortcuts
			   
-O   [STR]   Output prefix (required)

-1   [STR]   Trimmed read 1 (req)
-2   [STR]   Trimmed read 2 (paired, req)

-t   [INT]   Number of threads for alignment.
-o   [INT]   Number of threads for output.
               (also thread count for sorting by name)
-M   [#G]    GB used per thread in sorting (def = $sort_mem)

-X           Retain coord sorted bam (def is only name sorted)

-y           Perform alignment only and skip sorting by name

-Y           Resume sorting by name after coord sorted bam is generated (will check if bam exists)

Executable Commands (from $DEFAULTS_FILE)
   bsbolt:   $bsbolt
   samtools: $samtools
	
";

$start_time = localtime(time);

if (!defined $opt{'R'}) {
	die "\nERROR: Provide a reference as -R\n$die";
} else {
	if (defined $REF{$opt{'R'}}) {$ref = $REF{$opt{'R'}}}
	else {$ref = $opt{'R'}};
}

if (!defined $opt{'1'}) {die "\nERROR: Read 1 MUST be specified!\n$die"};
if (!defined $opt{'O'}) {die "\nERROR: Specify output as -O\n$die"};
if (defined $opt{'t'}) {$threads = $opt{'t'}};
if (defined $opt{'o'}) {$o_threads = $opt{'o'}};
if (defined $opt{'M'}) {$sort_mem = $opt{'M'}};

open LOG, ">>$opt{'O'}.bsbolt.log";

if (defined $opt{'2'}) {
	$align_call = "$bsbolt Align -F1 $opt{'2'} -F2 $opt{'1'} -t $threads -OT $o_threads -O $opt{'O'} -DB $ref >> $opt{'O'}.bsbolt.log 2>> $opt{'O'}.bsbolt.log";
} else {
	$align_call = "$bsbolt Align -F1 $opt{'1'} -t $threads -OT $o_threads -O $opt{'O'} -DB $ref >> $opt{'O'}.bsbolt.log 2>> $opt{'O'}.bsbolt.log";
}


print LOG "\n=== premethyst fastq-align ===\n";
print LOG "Command: $align_call\n\n\n";
if (!defined $opt{'Y'}) {
	system("$align_call");
	print LOG "\n\n\nAlignment done.\n";
} else {
	print LOG "Alignment skipped.\n";
	# Check if the BAM file exists and is not empty
	if (-s "$opt{'O'}.bam") {
		print LOG "$opt{'O'}.bam exists! Resume sorting by name.\n";
	} else {
		print LOG "$opt{'O'}.bam does NOT exist! Create empty bam for stub run.\n";
		system("touch $opt{'O'}.bam");
	}
}

$ts = localtime(time);
print LOG "Alignment complete for $opt{'O'}\nStart time: $start_time\nEnd time: $ts\n\n";

if (defined $opt{'y'}) {
	print LOG "Skip sorting by read name. Exit.\n";
} else {
	print LOG "Sorting by read name.\n";
	
	$sort_call = "$samtools sort -@ $o_threads -n -m $sort_mem $opt{'O'}.bam > $opt{'O'}.nsrt.bam";
	print LOG "Command: $sort_call\n";
	system("$sort_call");

	$te = localtime(time);
	print LOG "Sorting complete for $opt{'O'}\nStart time: $ts\nEnd time: $te\n\n";
}

if (!defined $opt{'X'}) {
	print LOG "Deleting $opt{'O'}.bam\n";
	system("rm -f $opt{'O'}.bam");
}

exit;

}

1;