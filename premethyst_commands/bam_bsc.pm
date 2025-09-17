package premethyst_commands::bam_bsc;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_bsc");


sub bam_bsc {

# defaults
$threads = 1;

$die = "

premethyst bam-bsc (options) [bam directory] 
      or   bsc

Call methytaion from duplicates-removed barcode-split bam.

Options    
-t   [INT]      Number of threads to use (def = $threads)
-d   [DIR]      BSBolt DB directory
-o   [DIR]      Output directory
-C              NOT only output CpG sites in CGmap file (def = CpG only)


Executable Commands (from $DEFAULTS_FILE)
   bsbolt:   $BSBolt
   premethyst: $premethyst
   
";

unless (getopts("t:d:o:C", \%opt)) {
    die "Unknown option or missing argument.\n$die";
}


if (!defined $ARGV[0]) {die "\nERROR: Specify bam directory as argument\n$die"};
if (!defined $opt{'d'}) {die "\nERROR: Specify BSBolt DB directory as -d\n$die"};
if (!defined $opt{'o'}) {die "\nERROR: Specify output directory as -o\n$die"};
if (defined $opt{'t'}) {$threads = $opt{'t'}};


$bsb_option = "$bsbolt CallMethylation -t $threads -min 1 -ignore-ov -verbose -DB $opt{'d'}";	
unless (defined $opt{'C'}) {$bsb_option .= " -CG";}

#-------------------------
# Configuration
#-------------------------
my $data_dir   = $ARGV[0];
my $out_dir  = $opt{'o'};

#-------------------------
# Prepare output dirs
#-------------------------
use File::Path qw(make_path);
make_path($out_dir) unless -d $out_dir;

#-------------------------
# Sort and move BAMs
#-------------------------
opendir(my $dh, $data_dir) or die "Can't open $data_dir: $!";
while (my $file = readdir($dh)) {
	next unless $file =~ /^(.+)\.rd.bam$/;
	my $cellID = $1;

	print STDERR "Call meth for $cellID...\n";
	my $rd_bam = "${data_dir}/${file}";
	my $out_prefix = "${out_dir}/${cellID}";
	my $call = "$bsb_option -I $rd_bam -O $out_prefix >> ${out_prefix}.bsbolt.call.log 2>> ${out_prefix}.bsbolt.call.log";
	print STDERR "$call\n";
	system("$call") == 0
		or die "Failed to CallMethylation for $cellID";
}
closedir($dh);

print STDERR "All done.\n";


}
1;
