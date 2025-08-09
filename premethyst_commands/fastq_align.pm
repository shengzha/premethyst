package premethyst_commands::fastq_align;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_align");

sub fastq_align {

getopts("O:1:2:t:o:R:XyYM:", \%opt);

$a_threads = 1;
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

-Y           Skip alignment and resume sorting by name (will check if coord sorted bam exists)

Executable Commands (from $DEFAULTS_FILE)
   bsbolt:   $bsbolt
   samtools: $samtools
	
";

# Check if both flags are on
if ( defined $opt{'Y'} && defined $opt{'y'} ) {
    die "\nERROR: Should not turn on -Y and -y at the same time!\n";
}


if (!defined $opt{'O'}) {die "\nERROR: Specify output as -O\n$die"};
if (defined $opt{'o'}) {$o_threads = $opt{'o'}}; # Apply to both aln and srt

if (defined $opt{'Y'}) { # Skip alignment and resume sorting by name
	# Check if the BAM file exists and is not empty
	if (! -f "$opt{'O'}.bam") {
		die "$opt{'O'}.bam does note exist!\n";
	}
} else { # Will perform alignment 
	if (!defined $opt{'R'}) {
		die "\nERROR: Provide a reference as -R\n$die";
	} else {
		if (defined $REF{$opt{'R'}}) {$ref = $REF{$opt{'R'}}}
		else {$ref = $opt{'R'}};
	}
	if (!defined $opt{'1'}) {die "\nERROR: Read 1 MUST be specified!\n$die"};
	if (defined $opt{'t'}) {$a_threads = $opt{'t'}}; # Only apply to bsbolt alignment

	# Construct cmd
	if (defined $opt{'2'}) {
		$align_call = "$bsbolt Align -F1 $opt{'1'} -F2 $opt{'2'} -t $a_threads -OT $o_threads -O $opt{'O'} -DB $ref >> $opt{'O'}.bsbolt.log 2>> $opt{'O'}.bsbolt.log";
	} else {
		$align_call = "$bsbolt Align -F1 $opt{'1'} -t $a_threads -OT $o_threads -O $opt{'O'} -DB $ref >> $opt{'O'}.bsbolt.log 2>> $opt{'O'}.bsbolt.log";
	}
}

if (defined $opt{'y'}) { # Perform alignment only and skip sorting by name
} else { # Will perform sorting 
	if (defined $opt{'M'}) {$sort_mem = $opt{'M'}}; # Only apply to sorting

	# Construct cmd
	$sort_call = "$samtools sort -@ $o_threads -n -m $sort_mem $opt{'O'}.bam > $opt{'O'}.nsrt.bam";
}



open LOG, ">>$opt{'O'}.bsbolt.log";

print LOG "\n=== premethyst fastq-align ===\n";


if (defined $opt{'Y'}) { # Skip alignment and resume sorting by name
	print LOG "Alignment skipped.\n";
} else { # Will perform alignment 
	
	print LOG "Command: $align_call\n\n\n";

	$ts = localtime(time);
	system("$align_call");
	$te = localtime(time);

	print LOG "\n\n\nAlignment done.\n";
	print LOG "Alignment complete for $opt{'O'}\nStart time: $ts\nEnd time: $te\n\n";

}

if (defined $opt{'y'}) { # Perform alignment only and skip sorting by name
	print LOG "Skip sorting by read name. Exit.\n";
} else { # Will perform sorting 
	print LOG "Command: $sort_call\n";
	$ts = localtime(time);
	system("$sort_call");
	$te = localtime(time);

	print LOG "Sorting complete for $opt{'O'}\nStart time: $ts\nEnd time: $te\n\n";

	# Can only delete coord sorted bam when name sorting is performed
	if (!defined $opt{'X'}) {
		print LOG "Deleting $opt{'O'}.bam\n";
		system("rm -f $opt{'O'}.bam");
	}

}

close LOG;

exit;

}

1;