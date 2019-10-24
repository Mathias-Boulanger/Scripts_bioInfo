#!/usr/bin/perl -w

use strict;

my $to_remove_from_file_title = "/home/chipseq/Bureau/Formation/Rnaseq/Count/"; # titre contient  
my $to_remove_from_file_title1 = "_hg19.count"; # titre contient  
my $to_remove_from_file_title2 = "_ERCC92.count"; # titre contient  
my @header; $header[0] = "Gene_id";



my %output_table;
for (my $i = 0; $i <= $#ARGV; $i++) {
    if (-e($ARGV[$i])) {
	open (F, $ARGV[$i]) or die "Cannot open file: ".$ARGV[$i]."\n";
	while(<F>) {
	    $_ =~ /^([^\t]+)\t([0-9]+)/;
	    $output_table{$1}[$i] = $2; # gene_id file_name = value
	}
	close F;
	my $title = $ARGV[$i];
	$title =~ s/$to_remove_from_file_title//;
	$title =~ s/$to_remove_from_file_title1//;
	$title =~ s/$to_remove_from_file_title2//;
	push @header, $title;
    } else {
	print STDERR "invalid parameter: ".$ARGV[$i].""; exit;
    }
}

print join("\t", @header)."\n";
foreach (sort {$a cmp $b;} keys %output_table) {
    print $_."\t";
    print join("\t", @{$output_table{$_}})."\n";
}
# paste <file1> <file2> <file3> <file4>| awk '{print $1,$2,$4,$8}'
