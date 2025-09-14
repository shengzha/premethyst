package premethyst_commands::bam_split;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_split");


sub bam_split {

$in_threads = 1;

$die = "

premethyst bam-split (options) -w [barcode whitelist] -o [output directory] [input bam]
      or   split-bam

Split name sorted bam by barcodes.

Options    
-w   [STR]      Barcode whitelist file, one barcode per line (required)
-o   [DIR]      Output directory 
-T   [INT]      Threads for the bam file read-in. (def = $in_threads)

Executable Commands (from $DEFAULTS_FILE)
   samtools:   $samtools
   premethyst: $premethyst
   
";

# The Getopt::Std module in Perl, which provides the getopts() function, 
# handles unknown options by issuing a warning and returning false. 
# To abort the program upon detection of an unknown option, the return value of getopts() can be checked, and die() can be called if it is false.
unless (getopts("w:o:T:", \%opt)) {
	die "Unknown option or missing argument.\n$die";
}


if (!defined $ARGV[0]) {die "\nERROR: Specify input as argument\n$die"};
if (!defined $opt{'w'}) {die "\nERROR: Specify barcode whitelist as -w\n$die"};
if (!defined $opt{'o'}) {$opt{'o'} = $ARGV[0]; $opt{'o'} =~ s/\.nsrt\.bam$/\.split/}; # parse from input bam
if (defined $opt{'T'}) {$in_threads = $opt{'T'}};

#-------------------------
# Configuration
#-------------------------
my $in_bam   = $ARGV[0];
my $out_dir  = $opt{'o'};    # Change this as needed
my $barcode_file = $opt{'w'};  # List of valid barcodes

#-------------------------
# Prepare output dirs
#-------------------------
use File::Path qw(make_path);
make_path($out_dir) unless -d $out_dir;

#-------------------------
# Read valid cell barcodes
#-------------------------
my %id_lookup;
open my $list_fh, '<', $barcode_file or die "Can't open barcode list: $!";
while (<$list_fh>) {
	chomp;
	$id_lookup{$_} = 1;
}
close $list_fh;

#-------------------------
# Extract SAM header
#-------------------------
my $header = `$samtools view -H $in_bam`;

#-------------------------
# Read BAM line by line
#-------------------------
open(my $IN, "-|", "$samtools view -@ $in_threads $in_bam 2>/dev/null") or die "Can't read BAM: $!";

my $prev_cellID = '';
my $cellID;
my $OUT;
while (my $line = <$IN>) {
	chomp $line;
	my @fields = split(/\t/, $line);

	# Extract cell ID (prefix before first colon in read name)
	($cellID, undef) = split(/:/, $fields[0], 2);


	# If switching to a new cell ID
	if ($cellID ne $prev_cellID) {

		# Skip if not in valid barcode list
		unless (exists $id_lookup{$cellID}) {
			print STDERR "$cellID is NOT in the list.\n";
			$prev_cellID = $cellID;  # still update so we don't repeat
			next;
		}

		# Close previous output
		if (defined $OUT) {
			close $OUT;
			print STDERR "Closed output for $prev_cellID\n";
		}

		# Open new output
		my $cell_bam = "$out_dir/${cellID}.bam";
		open($OUT, "|-", "$samtools view -bS - > $cell_bam") or die "Can't open $cell_bam: $!";
		print $OUT $header;

		print STDERR "Opened output for $cellID\n";
		$prev_cellID = $cellID;
	}

	# Write line if we’re writing to an output
	if (defined $OUT) {
		print $OUT "$line\n";
	}
}
close $IN;
if (defined $OUT) {
	close $OUT;
	print STDERR "Closed output for $cellID\n";
}

print STDERR "All done.\n";


}
1;
