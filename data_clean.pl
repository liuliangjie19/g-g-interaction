#!/usr/bin/perl -w
use strict;

sub clear {
  my @samstat = @_;
  my $i=0;
  foreach(@samstat){
    $i=$_;
    if ($i==-1 || $i==-2){
      return 0;
    }
  }
  return 1;
}

print "GGMU!\n";
my $gene1file = $ARGV[0];
my $gene2file = $ARGV[1];
my %gene1sam;
my %gene2sam;
my $line;
my $i=0;
my @linelist;
my @gene1list;
my @gene2list;

print "The data is about the gene $gene1file and the gene $gene2file!\n";
print "reading the infomation about the gene1\n";

open(GENE1, "<", $gene1file) || die "1 can not open the gene1 file\n";

while(<GENE1>){
  chomp;
  #print $_;
  if (/^Sample ID/){
    next;
  }
  $line = $_;
  @linelist = split('\t', $line);
  #print @linelist[2..$#linelist],$#linelist;
  #last;
  $gene1sam{$linelist[0]} = clear(@linelist[2..$#linelist]);
}

close(GENE1);
print "reading the infomation about the gene2\n";

open(GENE2, "<", $gene2file) || die "2 can not open the gene2 file\n";

while(<GENE2>){
  chomp;
  if (/^Sample ID/){
    next;
  }
  $line = $_;
  @linelist = split('\t', $line);
  $gene2sam{$linelist[0]} = clear(@linelist[2..$#linelist]);
}

close(GENE2);
open(GENE1, "<", $gene1file) || die "3 can not open the gene1 file\n";
open(GENE2, "<", $gene2file) || die "4 can not open the gene2 file\n";
@gene1list = <GENE1>;
@gene2list = <GENE2>;
close(GENE1);
close(GENE2);
open(GENE1, ">", "clear_".$gene1file) || die "5 can not open the gene1 file\n";
open(GENE2, ">", "clear_".$gene2file) || die "6 can not open the gene2 file\n";
my $sampleID;
for($i=0;$i<=$#gene2list;$i++){
  if ($i==0){
    #print $gene1list[$i],$gene2list[$i];
    next;
  }
  chomp($gene1list[$i]);
  chomp($gene2list[$i]);
  @linelist = split('\t', $gene1list[$i]);
  $sampleID = $linelist[0];
  if ($gene2sam{$sampleID}+$gene1sam{$sampleID}==2){
    print "$sampleID work out\n";
    print GENE1 "$gene1list[$i]\n";
    print GENE2 "$gene2list[$i]\n";
  }
}
close(GENE1);
close(GENE2);
