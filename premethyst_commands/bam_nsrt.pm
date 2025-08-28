package premethyst_commands::bam_nsrt;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_nsrt");


sub bam_nsrt {

# defaults
$sort_threads = 1;
$sort_mem = "4G";

$die = "

premethyst bam-nsrt (options) -O [out prefix (req)] [bam files]
      or   nsrt

Sorts bam by read name. Will merge the input bam files if there are multiple.

Options	   
-O   [STR]      Output prefix (required)
-t   [INT]      Number of threads for sorting.
-m   [#G]       GB used per thread in sorting (def = $sort_mem)
-T   [PREFIX]   Write temporary files to PREFIX.nnnn.bam

Executable Commands (from $DEFAULTS_FILE)
   samtools:   $samtools
   premethyst: $premethyst
   
";

# The Getopt::Std module in Perl, which provides the getopts() function, 
# handles unknown options by issuing a warning and returning false. 
# To abort the program upon detection of an unknown option, the return value of getopts() can be checked, and die() can be called if it is false.
unless (getopts("O:t:m:T:", \%opt)) {
    die "Unknown option or missing argument.\n$die";
}


if (!defined $ARGV[0]) {die "\nERROR: Specify input as argument\n$die"};
if (!defined $opt{'O'}) {die "\nERROR: Specify output as -O\n$die"};
if (defined $opt{'t'}) {$sort_threads = $opt{'t'}};
if (defined $opt{'m'}) {$sort_mem = $opt{'m'}};


$sort_options = "$samtools sort -n -o $opt{'O'}.nsrt.bam -@ $sort_threads -m $sort_mem";
if (defined $opt{'T'}) {
	$sort_options .= " -T $opt{'T'}";
}

if (@ARGV > 1) { # multiple bam files 
	$sort_call = "$samtools cat @ARGV | $sort_options";
} else {
	$sort_call = "$sort_options $ARGV[0]";
}


print STDERR "Command: $sort_call\n";
system("$sort_call");

$ts = localtime(time);
print STDERR "$ts\tDone.\n";


}

1;