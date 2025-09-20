package premethyst_commands::bam_split;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_split");


sub bam_split {

$threads = 1;
$n_buckets = 1;

$die = "

premethyst bam-split (options) -w [barcode whitelist] -o [output directory] [input bam]
      or   split

Split name sorted bam by barcodes.

Options    
-w   [STR]      Barcode whitelist file, one barcode per line (required)
-o   [DIR]      Output directory 
-t   [INT]      Threads for the bam file read-in. (def = $threads)
-b   [INT]      Distribute the bams to a number of buckets (when > 1) (def = $n_buckets)

Executable Commands (from $DEFAULTS_FILE)
   samtools:   $samtools
   premethyst: $premethyst
   
";

unless (getopts("w:o:t:b:", \%opt)) {
	die "Unknown option or missing argument.\n$die";
}


if (!defined $ARGV[0]) {die "\nERROR: Specify input as argument\n$die"};
if (!defined $opt{'w'}) {die "\nERROR: Specify barcode whitelist as -w\n$die"};
if (!defined $opt{'o'}) {$opt{'o'} = $ARGV[0]; $opt{'o'} =~ s/\.nsrt\.bam$/\.split/}; # parse from input bam
if (defined $opt{'t'}) {$threads = $opt{'t'}};
if (defined $opt{'b'}) {$n_buckets = $opt{'b'}};

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
open(my $IN, "-|", "$samtools view -@ $threads $in_bam 2>/dev/null") or die "Can't read BAM: $!";

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

		# 1. Close previous output
		if (defined $OUT) {
			close $OUT;
			print STDERR "Closed output for $prev_cellID\n";
		}

		# 2. Skip if not in valid barcode list
		unless (exists $id_lookup{$cellID}) {
			print STDERR "$cellID is NOT in the list.\n";
			$OUT = undef;
			$prev_cellID = $cellID;
			next;
		}

		# 3. Open new output
		my $cell_bam = "$out_dir/${cellID}.bam";
		open($OUT, "|-", "$samtools view -bS - > $cell_bam") 
			or die "Can't open $cell_bam: $!";
		print $OUT $header;
		print STDERR "Opened output for $cellID\n";

		$prev_cellID = $cellID;
	}

	# Write line if we have a valid output
	print $OUT "$line\n" if defined $OUT;
}
close $IN;
if (defined $OUT) {
	close $OUT;
	print STDERR "Closed output for $cellID\n";
}

print STDERR "Done splitting.\n";

unless ($n_buckets > 1) {return;}
print STDERR "\nRedistribute bams to $n_buckets buckets...\n";

use File::Basename;
use Digest::MD5 qw(md5_hex);
use File::Copy;


#------------------------
# Read BAMs
#------------------------
opendir(my $dh, $out_dir) or die "Cannot open directory '$out_dir': $!";
my @bam_files = grep { /\.bam$/ && -f "$out_dir/$_" } readdir($dh);
closedir($dh);


#------------------------
# Process each BAM file
#------------------------
foreach my $bam_file (@bam_files) {
	my ($cellID) = $bam_file =~ /^(.+)\.bam$/;

	my $subfolder;

	# Use hash bucket
	my $hash = md5_hex($cellID);
	my $bucket = hex(substr($hash, 0, 2)) % $n_buckets;
	$subfolder = "$out_dir/bucket_$bucket";
	

	make_path($subfolder) unless -d $subfolder;

	my $src = "$out_dir/$bam_file";
	my $dst = "$subfolder/$bam_file";

	print "Moving $bam_file → $dst\n";
	move($src, $dst) or warn "Failed to move $bam_file: $!";
}

print STDERR "All done.\n";

}
1;
