package premethyst_commands::bam_bss;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_bss");


sub bam_bss {

# defaults
$threads = 1;
$sort_mem = "4G";

$die = "

premethyst bam-bss (options) [bam directory] 
      or   bss

Sort name-sorted barcode-split bam by coord, and remove duplicates (optional).

Options    
-t   [INT]      Number of threads to use (def = $threads)
-m   [#G]       GB used per thread in sorting (def = $sort_mem)
-T   [PREFIX]   Write temporary files to PREFIX.nnnn.bam
-c              Perform sorting only
-I              Do NOT index the sorted bam


Executable Commands (from $DEFAULTS_FILE)
   samtools:   $samtools
   premethyst: $premethyst
   
";

# The Getopt::Std module in Perl, which provides the getopts() function, 
# handles unknown options by issuing a warning and returning false. 
# To abort the program upon detection of an unknown option, the return value of getopts() can be checked, and die() can be called if it is false.
unless (getopts("t:m:T:cI", \%opt)) {
    die "Unknown option or missing argument.\n$die";
}


if (!defined $ARGV[0]) {die "\nERROR: Specify bam directory as argument\n$die"};
if (defined $opt{'t'}) {$threads = $opt{'t'}};
if (defined $opt{'m'}) {$sort_mem = $opt{'m'}};


$fixm_options = "$samtools fixmate -@ $threads -p -m";	
$sort_options = "$samtools sort -@ $threads -m $sort_mem";
if (defined $opt{'T'}) {
	$sort_options .= " -T $opt{'T'}";
}
$md_options = "$samtools markdup -@ $threads -r -S";
if (defined $opt{'T'}) {
	$md_options .= " -T $opt{'T'}";
}

$idx_options = "$samtools index -@ $threads";

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

	my $nsrt_bam = "$data_dir/$file";
	my $csrt_bam   = "$data_dir/${cellID}.csrt.bam";

	print STDERR "Process BAM for $cellID...\n";

	
	my $call = "$fixm_options $nsrt_bam - | $sort_options - -o $csrt_bam";
	print STDERR "$call\n";
	system("$call") == 0
		or die "Failed to sort $nsrt_bam";

	unlink $nsrt_bam;  # Remove temporary file

	my $rd_bam = "$data_dir/${cellID}.rd.bam";
	unless (defined $opt{'c'}) {
		# remove duplicate reads
		my $stats = "$data_dir/${cellID}.rd.stats";
		my $call = "$md_options -f $stats $csrt_bam $rd_bam";
		print STDERR "$call\n";
		system("$call") == 0
			or die "Failed to markdup $csrt_bam";
	}
	
    if (defined $opt{'I'}) {next;} # skip index

	my $res_bam = (defined $opt{'c'}) ? $csrt_bam : $rd_bam;
	my $call = "$idx_options $res_bam";
	print STDERR "$call\n";
	system("$call") == 0
		or die "Failed to index $res_bam";


}
closedir($dh);

print STDERR "All done.\n";


}
1;
