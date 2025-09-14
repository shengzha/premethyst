package premethyst_commands::bam_bss;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_bss");


sub bam_bss {

# defaults
$sort_threads = 1;
$sort_mem = "4G";

$die = "

premethyst bam-bss (options) [bam directory] 
      or   bss

Sort barcode-splitted bam by coord.

Options    
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
unless (getopts("t:m:T:", \%opt)) {
    die "Unknown option or missing argument.\n$die";
}



if (!defined $ARGV[0]) {die "\nERROR: Specify bam directory as argument\n$die"};
if (defined $opt{'t'}) {$sort_threads = $opt{'t'}};
if (defined $opt{'m'}) {$sort_mem = $opt{'m'}};

$sort_options = "$samtools sort -@ $sort_threads -m $sort_mem";
if (defined $opt{'T'}) {
	$sort_options .= " -T $opt{'T'}";
}

#-------------------------
# Configuration
#-------------------------
my $data_dir   = $ARGV[0];

#-------------------------
# Sort and move BAMs
#-------------------------
opendir(my $dh, $data_dir) or die "Can't open $data_dir: $!";
while (my $file = readdir($dh)) {
	next unless $file =~ /^(.+)\.bam$/;
	my $cellID = $1;

	my $unsorted_bam = "$data_dir/$file";
	my $sorted_bam   = "$data_dir/${cellID}.csrt.bam";

	print STDERR "Sorting BAM for $cellID...\n";
	system("$sort_options -o $sorted_bam $unsorted_bam") == 0
		or die "Failed to sort $unsorted_bam";

	unlink $unsorted_bam;  # Remove temporary file

	# # Optionally index the BAM
	# system("$samtools index $sorted_bam") == 0
	#     or warn "Failed to index $sorted_bam";

}
closedir($dh);

print STDERR "All done.\n";


}
1;
