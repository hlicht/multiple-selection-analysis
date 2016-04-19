#!/usr/bin/perl -w
#
#	This script takes two sequence files in tab format and align the orthologus
#	sequences and estimate the pairwise dN/dS ratio.
#	
#	Input file needs to be in fasta-tab format, and all letters in Upper-case!
#	
#	For full usage info type:
#	perl AlignAnddNdSAnalysis.pl
#
#	This script was written by Henrik H. De Fine Licht (hhdefinelicht@plen.ku.dk)
#	No copyright or license
#
################################################################
use strict;


my $bbhfile;
my $emfile;
my $drfile;
my $framefile;



my $emfastaid;
my $drfastaid;
my %bbhhash;
my %emseqhash;
my %drseqhash;
my $bbhid;
my $emid;
my $drid;
my $emseq;
my $drseq;
my $emcontigid;
my $drcontigid;
my $emid_frame;
my $drid_frame;
my $frameinfo;
my %framehash;

my $clseq;
my $lgcontig;
my $clcontig;
my $lgseq;
my $length;
my $clustalline;
my @filearray;
my $newlength;
my $last;
my $first;
my $seclast;
my $codemlline;
my $alnlength;
my $dNdS;
my $dN;
my $dS;
my $newlast;
my $newseclast;
my $line;
my @seqarray;
my @lgseqarray;
my $newlgseq;
my @clseqarray;
my $newclseq;
my $newlength2;
my $t;
my $N;
my $S;


# First deal with passed parameters
#########################################################################################
if ($#ARGV == -1) {
    &usage;
}
$bbhfile = "";
$emfile = "";
$drfile = "";
$framefile = "";

my %my_args = @ARGV;
for my $i (sort keys %my_args) {
    if ($i eq "-b") {
	$bbhfile = $my_args{$i};
    }
    elsif ($i eq "-p") {
	$emfile = $my_args{$i};
    }
    elsif ($i eq "-s") {
	$drfile = $my_args{$i};
    }
    elsif ($i eq "-f") {
	$framefile = $my_args{$i};
    }
    else {
	print "\nUnrecognized argument: $i\n\n";
	&usage;
    }
}

open (BBHFILE,'<', "$bbhfile") or die "$!";

# Open the sorted tab files and store them in a hash with contig as key and the
# two seqs as values
open (EMFILE,'<', "$emfile") or die "$!";
open (DRFILE,'<', "$drfile") or die "$!";

#Open special tab file with frame info with Eminfo first
open (FRAMEFILE,'<', "$framefile") or die "$!"; 



# Start parsing the files
########################################################################################




while ( my $bbhline = <BBHFILE> ) {
    chomp $bbhline;
    ($bbhid, $emid, $drid) = split ' ', $bbhline;
    push @{ $bbhhash{$bbhid} },($emid,$drid);
}

while ( my $emline = <EMFILE> ) {
    chomp $emline;
    ($emfastaid, $emseq) = split ' ', $emline;
    $emseqhash{$emfastaid} =  $emseq;
}

while ( my $drline = <DRFILE> ) {
    chomp $drline;
    ($drfastaid, $drseq) = split ' ', $drline;
    $drseqhash{$drfastaid} = $drseq;

}

while ( my $frameline = <FRAMEFILE> ) {
    chomp $frameline;
    ($emid_frame, $frameinfo, $drid_frame) = split ' ', $frameline;
    $framehash{$emid_frame} = $frameinfo;
}


close EMFILE;
close DRFILE;
close BBHFILE;

# create a temporary fasta file of each sequence pair and send it to PRANK
########################################################################################

foreach my $hashkey ( sort keys %bbhhash ) {
    ($emcontigid, $drcontigid) = @{ $bbhhash{$hashkey} };
    next if ( $emcontigid eq '' | $drcontigid eq '');
    open (TEMPFASTA,'>tempcontig.fasta') or die "$!";
    if ( $framehash{$emcontigid} =~ /-1/ ) {
	#reverse the DNA sequence
	my $revcomp = reverse($drseqhash{$drcontigid});
	#complement the DNA sequence
	$revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
	print TEMPFASTA ">Em$emcontigid\n$emseqhash{$emcontigid}\n>Dr$drcontigid\n$revcomp\n";
	close TEMPFASTA;
    } 
    
    elsif ( $framehash{$emcontigid} =~ /-2/ ) {
	#reverse the DNA sequence
	my $revcomp = reverse($drseqhash{$drcontigid});
	#complement the DNA sequence
	$revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
	#remove the first character form the string to get the correct frame
	substr($revcomp, 0, 1) = ""; 
	print TEMPFASTA ">Em$emcontigid\n$emseqhash{$emcontigid}\n>Dr$drcontigid\n$revcomp\n";
	close TEMPFASTA;
    } 
    
    elsif ( $framehash{$emcontigid} =~ /-3/ ) {
	#reverse the DNA sequence
	my $revcomp = reverse($drseqhash{$drcontigid});
	#complement the DNA sequence
	$revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
	#remove the first two characters form the string to get the correct frame
	substr($revcomp, 0, 2) = "";
 	print TEMPFASTA ">Em$emcontigid\n$emseqhash{$emcontigid}\n>Dr$drcontigid\n$revcomp\n";
	close TEMPFASTA;
    } 
    
    elsif ( $framehash{$emcontigid} =~ /\+1/ ) {
	print TEMPFASTA ">Em$emcontigid\n$emseqhash{$emcontigid}\n>Dr$drcontigid\n$drseqhash{$drcontigid}\n";
	close TEMPFASTA;
    } 
    
    elsif ( $framehash{$emcontigid} =~ /\+2/ ) {
	substr($drseqhash{$drcontigid}, 0, 1) = "";
	print TEMPFASTA ">Em$emcontigid\n$emseqhash{$emcontigid}\n>Dr$drcontigid\n$drseqhash{$drcontigid}\n";
	close TEMPFASTA;
    } 
    
    elsif  ( $framehash{$emcontigid} =~ /\+3/ ) {
	substr($drseqhash{$drcontigid}, 0, 2) = "";
	print TEMPFASTA ">Em$emcontigid\n$emseqhash{$emcontigid}\n>Dr$drcontigid\n$drseqhash{$drcontigid}\n";
	close TEMPFASTA;
    } 
    
    else {
	# Print the sequence as is if no frame info, should not be the case because all BBH 
	# pairs has a hit!!!!
	print TEMPFASTA ">Em$emcontigid\n$emseqhash{$emcontigid}\n>Dr$drcontigid\n$drseqhash{$drcontigid}\n";
	close TEMPFASTA;
    }

# Get sequence on a single line in the fasta file
#########################################################################################
    open (FASTASEQS,'<',"tempcontig.fasta") or die "$!";
    open (FASTAONELINE,'>',"tempcontigoneline.fasta") or die "$!";


    @filearray = <FASTASEQS> ;
    foreach $line ( @filearray ) {
	unless ( $line =~ /^>/ ) {
	    $line =~ s/\s//g;
	    chomp $line;
	}
	
	$line =~ s/>Dr/\n>Dr/; # Start of header for sequence 2 IMPORTANT!!!!
    }
    
    print FASTAONELINE @filearray;
    close FASTASEQS;
    close FASTAONELINE;

   
    if ( system( "rm tempcontig.fasta" ) != 0 ) {
	die "error deleting tempcontig.fasta!!!\n";
    }


# remove trailing basepairs so the two seqs can be divided by three
    open (FASTAONELINE,'<',"tempcontigoneline.fasta") or die "$!";
    open (FASTAMODIF,'>',"tempcontigmodif.fasta") or die "$!";

    @filearray =  <FASTAONELINE>;
    print "@filearray\n"; # testcode to print temporary files for troubleshooting

    chomp $filearray[1];
    $length = length $filearray[1];
    if ( $length % 3 == 0 ) {
	$filearray[1] = $filearray[1] . "\n";
#	print FASTAMODIF @filearray; # Unmask for troubleshooting
	system( "cp tempcontigmodif.fasta tempcontigmodif_inframe1.fasta" );
    } 

    elsif ( $length % 3 == 1 ) {
	chop $filearray[1];
	$filearray[1] = $filearray[1] . "\n";
#	print FASTAMODIF @filearray; # Unmask for troubleshooting
	system( "cp tempcontigmodif.fasta tempcontigmodif_inframe2.fasta" );
    } 
    
    elsif ( $length % 3 == 2 ) {
	chop $filearray[1];
	chop $filearray[1];
	$filearray[1] = $filearray[1] . "\n";
#	print FASTAMODIF @filearray; # Unmask for troubleshooting
	system( "cp tempcontigmodif.fasta tempcontigmodif_inframe3.fasta" );
    } 

    else {
	die "problem parsing temcontigmodif.fasta: $!\n";
    }
    
    chomp $filearray[3];
    $length = length $filearray[3];
    if ( $length % 3 == 0 ) {
	$filearray[3] = $filearray[3] . "\n";
	print FASTAMODIF @filearray;
	system( "cp tempcontigmodif.fasta tempcontigmodif_inframe1.fasta" );
    } 

    elsif ( $length % 3 == 1 ) {
	chop $filearray[3];
	$filearray[3] = $filearray[3] . "\n";
	print FASTAMODIF @filearray;
	system( "cp tempcontigmodif.fasta tempcontigmodif_inframe2.fasta" );
    } 

    elsif ( $length % 3 == 2 ) {
	chop $filearray[3];
	chop $filearray[3];
	$filearray[3] = $filearray[3] . "\n";
	print FASTAMODIF @filearray;
	system( "cp tempcontigmodif.fasta tempcontigmodif_inframe3.fasta" );
    } 

    else {
	die "problem parsing temcontigmodif.fasta: $!\n";
    }


   
    if ( system( "rm tempcontigoneline.fasta" ) != 0 ) {
	die "error deleting tempcontigoneline.fasta!!!\n";
    }

# Align sequences using PRANK
#
# It may be necessary to give the full path to the prank program. That is change "prank" 
# to the full path, ex. "/home/programs/prank"
#
########################################################################################

    if ( system( "prank -d=tempcontigmodif.fasta -o=prankout -f=paml -codon -once" ) != 0 ) {
	die "error executing prank command!!!\n";
    }

 
    system( "cp prankout.best.phy prankout.best_check.phy");# testcode to print temporary
    # files for troubleshooting
    system( "mv prankout.best.phy clustalout.phy");
    
    open (CLUSTALALN,'<',"clustalout.phy") or die "$!";
    open (CLUSTALFORMAT,'>',"clustalformat.phy") or die "$!";

    # Get sequence on a single line in the alignment file
    @filearray = <CLUSTALALN> ;
    foreach $line ( @filearray ) {
	unless ( $line =~ /\d+/ ) {
	    $line =~ s/\s//g;
	    chomp $line;
	}
	
	$line =~ s/Dr/\nDr/; # Start of header for sequence 2 IMPORTANT!!!!
    
    }
    
    print CLUSTALFORMAT @filearray;
    close CLUSTALALN;
    close CLUSTALFORMAT;
    system( "cp clustalformat.phy clustalformat_check.phy");# testcode to print temporary
    # files for troubleshooting
    system( "mv clustalformat.phy clustalout.phy");

    open (CLUSTALALN2,'<',"clustalout.phy") or die "$!";
    open (CLUSTALFORMAT2,'>',"clustalformat.phy") or die "$!";

    # if the primary sequence file is aligned with the first positions as gaps '-'
    # these positions should be removed to maintain the alignment in frame
    my $nucl = 0;
    @filearray =  <CLUSTALALN2> ;
#    print "$filearray[2]\n"; # Unmask for troubleshooting
    if ( $filearray[2] =~ /^-/ ) {
#	print "NOTICE THIS NOTICE THIS\n"; # Unmask for troubleshooting
	@lgseqarray = split '', $filearray[2];
#	print "@lgseqarray"; # Unmask for troubleshooting
	foreach $line ( @lgseqarray ) {
	    if ( $line =~ /-/ ) {
#		print "$line\n"; # Unmask for troubleshooting
		$nucl = $nucl + 1;
	    } else {
		last;
	    }
	}


# Now remove the positions from the alignment file
	splice @lgseqarray, 0, $nucl;
	$newlgseq = join '', @lgseqarray;

	@clseqarray = split '', $filearray[4];
	splice @clseqarray, 0, $nucl;
	$newclseq = join '', @clseqarray;

# And change the header of the phylip format alignment file
	$filearray[0] =~ /2\s+(\d+)/;
	$length = $1;
	$newlength2 = $length - $nucl;
	$filearray[0] =~ s/(2\s+)\d+/$1$newlength2/;

	
	print CLUSTALFORMAT2 "$filearray[0]$filearray[1]$newlgseq$filearray[3]$newclseq";
	close CLUSTALFORMAT2;
	system( "cp clustalformat.phy clustal_fixlen"); # testcode to print temporary
    # files for troubleshooting
	system( "mv clustalformat.phy clustalout.phy");
    }
       
    close CLUSTALALN2;


    
# remove trailing basepairs so the alignment can be divided by three
    open (CLUSTALALN3,'<',"clustalout.phy") or die "$!";
    open (CLUSTALFORMAT3,'>',"clustalformat.phy") or die "$!";

    @filearray =  <CLUSTALALN3>;
    if ( $filearray[0] =~ /2\s+\d+/ ) {
	$filearray[0] =~ /2\s+(\d+)/;
	$length = $1;
	if ( $length % 3 == 0 ) {
	    print CLUSTALFORMAT3 @filearray;
	    system( "cp clustalformat.phy clustalformat_inframe.phy" );
	} 
	
	elsif ( $length % 3 == 1 ) {
	    $newlength = $length - 1;
#		@filearray = <CLUSTALALN> ;
	    $last = pop (@filearray);
	    $seclast = $filearray[2];
	    chomp $last;
	    chop $last;
	    $newlast = $last . "\n";
	    chomp $seclast;
	    chop $seclast;
	    $newseclast = $seclast . "\n";
	    $filearray[2] = $newseclast;
	    push @filearray, $newlast;
		
	    $first = shift (@filearray);
	    $first =~ s/(2\s+)\d+/$1$newlength/;
	    unshift @filearray, $first;
	
	    print CLUSTALFORMAT3 @filearray;
	    system( "cp clustalformat.phy clustalformat_minus1.phy" );
	} 
	
	elsif ( $length % 3 == 2 ) {
	    $newlength = $length - 2;
#		@filearray = <CLUSTALALN> ;
	    $last = pop (@filearray);
	    $seclast = $filearray[2];
	    chomp $last;
	    chop $last;
	    chop $last;
	    $newlast = $last . "\n";
	    chomp $seclast;
	    chop $seclast;
	    chop $seclast;
	    $newseclast = $seclast . "\n";
	    $filearray[2] = $newseclast;
	    push @filearray, $newlast;
		
	    $first = shift (@filearray);
	    $first =~ s/(2\s+)\d+/$1$newlength/;
	    unshift @filearray, $first;
		
	    print CLUSTALFORMAT3 @filearray;
	    system( "cp clustalformat.phy clustalformat_minus2.phy" );
	} 

	else {
	    die "problem parsing clustal phylip alignment: $!\n";
	}
    }

	# replace the original phylip alignment with the formatted clustal alignment 
	# and send it to codeml.
    system( "mv clustalformat.phy clustalout.phy");
 
    system( "/home/henrik/Programs/paml4.8/bin/codeml" ); # != 0 );


    system( "cat output_pairwise >> output_pairwise_concat" ); # testcode to print temporary
    # files for troubleshooting
    system( "cat clustalout.phy >> clustalout_phy_concat" ); # testcode to print temporary
    # files for troubleshooting
    
    

    # Parse the codeml output and write output file
    
    open (CLUSTALALN4,'<',"clustalout.phy") or die "$!";
    @filearray =  <CLUSTALALN4>;
    if ( $filearray[0] =~ /2\s+\d+/ ) {
	$filearray[0] =~ /2\s+(\d+)/;
	$alnlength = $1;
    }

    open (CODEMLOUTPUT,'<', "output_pairwise") or die "Problem reading codeml output file: $!";
    open (RESULTFILE,'>>', "dNdS_analysis_output.tab");
	
    print RESULTFILE "$hashkey\t$emcontigid\t$drcontigid\t";
    while ( $codemlline = <CODEMLOUTPUT> ) {
	if ( $codemlline =~ /^After\sdeleting\sgaps/ ) {
	    $codemlline =~ /(\d+)\ssites/;
	    $alnlength = $1;
	} 

	elsif ( $codemlline =~ /^t=\s/ ) {
	    chomp $codemlline;
	    print RESULTFILE "$codemlline\t";
	}	
    
    
    }
    print RESULTFILE "length= $alnlength\n";

    if ( system( "rm clustalout.phy" ) != 0 ) {
	die "error deleting clustalout.phy!!!\n";
    }

}

    



#-----------------------------------------------------------------------
sub usage {
    print "\nUsage: ./AlignAnddNdSAnalysis.pl -p primary seq file -s secondary seq file -b bbh file -f framefile\n\n";
    print "Parameters:\n";
    print "-p primary sequence file in tab format and in frame. Uppercase letters.\n\n";
    print "-s secondary sequence file in tab format. Uppercase letters.\n\n";
    print "-b bbh input file\n\n";
    print "-f frame file\n\n";
    print "bbh input file: File with orthologous sequence information (for example from Bidirectional Best BLAST hits). This is a tab file with three columns separated by tabs with a unique number in the first column followed by the sequence names of the orthologous sequence pair to be tested on each line with the primary sequence name in the second column and secondary sequence in the third column.\n\n";
    print "frame file: Tab file with frame info from the bidirectional best blast hits. First column is sequence name of primary sequence, second column is the frame info (ex. '+1') and third column is the secondary sequence name.\nPossible frames allowed: +1, +2, +3, -1, -2, -3\n\n";
    print "Output will be written to the file dnds_analysis_output.tab.\n\n";
    print "Columns in output file:\n";
    print "1. Orthologous pair number from bbh file\n";
    print "2. Primary sequence name\n";
    print "3. Secondary sequence name\n";
    print "4. t\n5. S\n6. N\n7. dNdS\n8. dN\n9. dS\n10. Alignment length\n\n";
    print "Dependencies:\n";
    print "Executable program files of prank and codeml, including the codeml.ctl file needs to be in the same library!!!!!.\n\n";
    print "A number of temporary files are written that can be deleted after a succesful run.\n";
    print "Especially the files 'output_pairwise_concat' and 'clustalout_phy_concat' are good to double-check results afterwards, as they contain the codeml output and alignments, respectively\n\n";

    exit;
}
#-----------------------------------------------------------------------





