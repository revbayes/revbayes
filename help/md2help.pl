use strict;
use warnings;
use List::Util qw[min max];


sub find_block_indent {
	my ($arg) = @_;
        my @lines = split(/\n/,$arg);

        my $indent;
        foreach my $line1 (@lines)
        {
                my $line2 = $line1;
                $line2 =~ s/^ +//;
                my $whitewidth = length($line1) - length($line2);
                next if ($line2 eq "");

                if (defined($indent)) {
                        $indent = min($indent,$whitewidth);
                } else {
                        $indent = $whitewidth;
                }
        }
        return $indent;
}

sub parse_entry {
	my ($key, $value) = @_;

	my $ret;

	# return array reference for see_also and authors entries
	if ( $key eq "see_also" || $key eq "authors") {
		@{$ret} = ();
		my @items = split '\n', $value;
		foreach my $item (@items) {
			# remove leading and trailing spaces
			$item =~ s/(^\s+|\s+$)//g;
			push @{$ret}, $item;
		}
	# parse yaml entry for each reference
	# fill return array reference with a hash for each ref
	} elsif ( $key eq "references" ) {
		@{$ret} = ();
		$value =~ s/^\s*- //;
		my @items = split '- ', $value;
		foreach my $item (@items) {
			
			my %reference;
			$reference{'citation'} = "";
			$reference{'doi'} = "";
			$reference{'url'} = "";

			$item =~ s/^\s*//mg;
			$item =~ s/\\?\n\s*\\?/ /sg;
			$item =~ s/ (?=(citation|doi|url))/\n/g;
			while($item =~ /(citation|doi|url): (.*)$/mg) {
				$reference{$1} = $2;
				$reference{$1} =~ s/null\s*//;
			}

			push @{$ret}, \%reference;
		}
	# return all other fields as simple strings 
	} else {
		$ret = $value;

                # Replace tabs with 8 spaces
                $ret =~ s/\t/        /g;

		# remove block indentation from examples
		if ( $key eq "example")
                {
                        my $indent = find_block_indent($ret);

                        if (defined($indent))
                        {
                            my @lines = split(/\n/,$ret);
                            foreach my $line (@lines)
                            {
                                $line = substr($line,$indent);
                                $line =~ s/\s*$//s;
                            }
                            $ret = join("\n",@lines);
                        }
		}

		# remove trailing spaces
		$ret =~ s/\s*$//s;
		#$ret =~ s/\n/\\n/sg;
	}
	return $ret;
}

# print .cpp preample
print "/**\n";
print " * This file was generated automatically.\n";
print " * It is not intended to be human readable.\n";
print " * See help/README.md for details.\n";
print " */\n\n";
print "#include \"RbHelpDatabase.h\"\n\n";
print "using namespace std;\n\n";
print "RevBayesCore::RbHelpDatabase::RbHelpDatabase()\n{\n";

# go through each markdown file
foreach my $file (sort @ARGV) {
	$file =~ /([^\/]+)$/;
	my $name = $1;
	$name =~ s/\.md$//;

	open FILE, $file;
	my $key;
	my $value = "";
	my %entry;
	while (<FILE>) {
		# find each markdown header
		if ($_ =~ /^\#\# (.*)$/) {
			if ( $key ) {
				$entry{$key} = parse_entry($key, $value);
			}

			$key = $1;
			$value = "";
		} else {
			$value .= $_;
		}
	}
	close FILE;

	$entry{$key} = parse_entry($key, $value);

	# print .cpp map entries for each markdown field
	foreach my $key (sort keys %entry) {
		if ( ref($entry{$key}) eq 'ARRAY' ) {
			my @cpp;
			foreach my $item ( @{$entry{$key}} ) {
				if ( ref($item) eq 'HASH' ) {
					my %hash = %{$item};
					print "\thelp_references[string(\"$name\")].push_back(RbHelpReference(R\"($hash{'citation'})\",R\"($hash{'doi'})\",R\"($hash{'url'})\"));\n";
				} else {
					print "\thelp_arrays[string(\"$name\")][string(\"$key\")].push_back(string(R\"($item)\"));\n";
				}
			}
		} else {
			$entry{$key} =~ s/\s+$//g;
			if ( $entry{$key} ne "" )
			{
				print "\thelp_strings[string(\"$name\")][string(\"$key\")] = string(R\"($entry{$key})\");\n";
			}
		}
	}
}

print "}\n";
