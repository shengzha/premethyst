package premethyst_commands::fastq_align;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_align");

sub fastq_align {

getopts("O:1:2:t:o:R:", \%opt);

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

Options:

-R   [STR]   Reference path (required)
    Shortcuts:
$ref_shortcuts
			   
-O   [STR]   Output prefix (required)

-1   [STR]   Trimmed read 1 (req)
-2   [STR]   Trimmed read 2 (paired, req)

-t   [INT]   Number of threads for alignment.
-o   [INT]   Number of threads for output.

Executable Commands (from $DEFAULTS_FILE)
   bsbolt:   $bsbolt
   samtools: $samtools
	
";


if (!defined $opt{'R'}) {
	die "\nERROR: Provide a reference as -R\n$die";
} else {
	if (defined $REF{$opt{'R'}}) {$ref = $REF{$opt{'R'}}}
	else {$ref = $opt{'R'}};
}

if (!defined $opt{'1'}) {die "\nERROR: Read 1 MUST be specified!\n$die"};
if (!defined $opt{'O'}) {die "\nERROR: Specify output as -O\n$die"};
if (defined $opt{'t'}) {$a_threads = $opt{'t'}};
if (defined $opt{'o'}) {$o_threads = $opt{'o'}};

# Construct cmd
$align_call = "$bsbolt Align -F1 $opt{'1'}";
$align_call .= " -F2 $opt{'2'}" if defined $opt{'2'};
$align_call .= " -t $a_threads -OT $o_threads -O $opt{'O'} -DB $ref >> $opt{'O'}.bsbolt.log 2>> $opt{'O'}.bsbolt.log";

open LOG, ">$opt{'O'}.bsbolt.log";
$ts = localtime(time);
print LOG "$ts\tAlignment called.\n";
print LOG "Command: $align_call\n";
system("$align_call");

$ts = localtime(time);
print LOG "$ts\tDone.\n";


close LOG;

exit;

}

1;