#! /bin/perl

# Creates a hierarchy file for McClintock from a GFF created by RepeatMasker and
# the fasta TE database used to run RepeatMasker.

# Open the GFF file output from a RepeatMasker run 
open (GFFIN, $ARGV[0]);
# Open the TE Fasta consensus sequence database used with RepeatMasker
open (FASTAIN, $ARGV[1]);
# Open a Gff to be the output GFF file with corrected IDs
open (GFFOUT, ">$ARGV[0]_ID");
# Open a file that will contain a ready made hierachy file for McClintock
open (HIERARCHY, ">$ARGV[2]");
# Open a Fasta to be the output Fasta file with corrected IDs
open (FASTAOUT, ">$ARGV[1]_ID");

# For each line in the fasta file
while (my $TE = <FASTAIN>)
{
	# Reset to the top of the GFF input
	seek (GFFIN, 0,0);
	# If this line is the header of a TE sequence
	if ($TE =~ m/^>/)
	{
		# Save the whole ID, split on the RepeatMasker dividing character
		my @teID = split (/[>#]/, $TE);
		my $family = $teID[1];
		
		print FASTAOUT ">$family";

		my $family_count = 1;
		# For every GFF record in the input
		while (my $GFF_record = <GFFIN>)
		{
			# If the record is not header line
			if ($GFF_record !~ m/^#/)
			{
				my @GFF_fields = split (/\t/, $GFF_record);
				my @exact_family = split (/[:"]/, $GFF_fields[8]);
				chomp($family);
				# If the search family matches the family of the GFF record
				if ($exact_family[2] eq  $family)
				{
					# Print the GFF fields to the new GFF
					print GFFOUT join("\t", @GFF_fields[0 .. 1]), "\t";
					print GFFOUT "$family\t";
					print GFFOUT join("\t", @GFF_fields[3 .. 7]), "\t";
					# Print the family name as the ID with a count appended to make it unique
					print GFFOUT "ID=$family","_","$family_count\n";
					# Print the TE ID and family name to the hierarchy file
					print HIERARCHY "$family","_","$family_count\t$family\n";
					$family_count++;
				} 
			}
		}
	}
	else
	{
		print FASTAOUT "$TE";
	}
}

